"""
brillouin.jl — Spontaneous Brillouin scattering noise.

Included into the top-level `BIFROST` module by `BIFROST.jl` (not its own
module). Reuses `_simpson` / `_trapz` and `SPEED_OF_LIGHT_M_PER_S` from earlier
includes.

Backward (counter-propagating) Brillouin only — forward (GAWBS) is ~30 dB
weaker and is not modelled here.

# Fiber-aware top-level API

```julia
spbs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump; T_K=fiber.T_ref_K, ...)
brillouin_threshold(fiber::Fiber, λ_pump; T_K=fiber.T_ref_K)
```

# References
- Agrawal, *Nonlinear Fiber Optics*, 6th ed. (2019), Ch. 9.
- Smith, R.G., *Appl. Opt.* 11, 2489 (1972).
- Nikles, M. et al., *J. Lightwave Technol.* 15, 1842 (1997).
- Kobyakov, A. et al., *Adv. Opt. Photon.* 2, 1 (2010).
"""

# ── Physical constants (private) ──────────────────────────────────────────────
# `SPEED_OF_LIGHT_M_PER_S` is in scope from material-properties.jl.
const _HBAR_BRIL = 1.054571817e-34
const _KB_BRIL   = 1.380649e-23

# ── Material constants (tabulated bulk values, no fits) ──────────────────────
# All values are tabulated for pure SiO₂ or pure GeO₂ glass at 297 K.
# References:
#   ρ_SiO₂   = 2200 kg/m³  — CRC Handbook of Chemistry & Physics
#   ρ_GeO₂   = 3650 kg/m³  — CRC Handbook
#   C₁₁_SiO₂ = 78.5 GPa    — Boyd, Nonlinear Optics 4e Table 9.2
#   C₁₁_GeO₂ = 60.0 GPa    — Kobyakov, Sauer & Chowdhury, Adv. Opt. Photon. 2, 1 (2010), §2.3
#   M_SiO₂   = 60.08 g/mol — molar mass
#   M_GeO₂   = 104.61 g/mol
const _RHO_SIO2 = 2200.0      # kg/m³
const _RHO_GEO2 = 3650.0      # kg/m³
const _C11_SIO2 = 78.5e9      # Pa  — longitudinal elastic stiffness
const _C11_GEO2 = 60.0e9      # Pa
const _M_SIO2   = 60.08       # g/mol
const _M_GEO2   = 104.61      # g/mol

# Measured Brillouin linewidth for fused silica at 297 K and ~1550 nm pump.
# This is a measured-quantity material datum, not a fit. Source: Agrawal NLFO 6e Ch. 9.
const BRIL_DFREQ_HZ  = 25.0e6
const BRIL_GAMMA_RAD = 2π * BRIL_DFREQ_HZ

# Smith (1972) SBS threshold convention. Derived from solving the coupled SBS
# equations with quantum-noise boundary conditions; not a fit.
const BRIL_G_THRESH  = 21.0

# ── Binary-glass mixing rules (textbook linear mass-weighted mixing) ─────────

"Mass fraction of GeO₂ for a given molar fraction `x` of GeO₂."
mass_fraction_geo2(x::Real) = x * _M_GEO2 / ((1 - x) * _M_SIO2 + x * _M_GEO2)

"""
    glass_density(x::Real) → kg/m³

Density of a SiO₂–GeO₂ binary glass at molar fraction `x` of GeO₂.
Linear mass-weighted mixing between bulk endpoint values.
"""
function glass_density(x::Real)
    y = mass_fraction_geo2(x)
    return (1 - y) * _RHO_SIO2 + y * _RHO_GEO2
end

"""
    glass_c11(x::Real) → Pa

Longitudinal elastic stiffness `C₁₁` of a SiO₂–GeO₂ binary glass at molar
fraction `x` of GeO₂. Linear mass-weighted mixing between bulk endpoint values.
"""
function glass_c11(x::Real)
    y = mass_fraction_geo2(x)
    return (1 - y) * _C11_SIO2 + y * _C11_GEO2
end

"""
    acoustic_velocity_bulk(x::Real = 0.0) → m/s

Bulk longitudinal acoustic velocity in a SiO₂–GeO₂ binary glass from
`v_A = √(C₁₁ / ρ)`. Returns the *bulk* material value; for guided-fibre
applications the effective acoustic velocity is lower (~10 % for SMF-28) due
to confinement — implement Kobyakov 2010 §3 to get the guided correction.
"""
acoustic_velocity_bulk(x::Real=0.0) = sqrt(glass_c11(x) / glass_density(x))

# Backward-compatible alias for the previous constant — now derived, not tabulated.
const BRIL_VA_SIO2 = acoustic_velocity_bulk(0.0)    # ≈ 5973 m/s

# ── Brillouin peak gain — anchored at a measured material datum ─────────────
#
# The "first-principles" Agrawal NLFO 6e Eq. (9.1.7) formula
#   g_B = (prefactor) · γ_e² / (c · λ_p² · ρ · v_A · Δν_FWHM)
# has at least three published prefactor conventions and a γ_e definition
# (γ_e = ρ ∂ε/∂ρ vs n^4·p₁₂) that differ by O(1) factors. Picking any one and
# computing absolute values gives g_B for pure silica in the range
# 2.5×10⁻¹¹ – 2×10⁻¹⁰ m/W — none of which match the universally-cited measured
# value of 5×10⁻¹¹ m/W.
#
# Even worse, Eq. 9.1.7 with bulk mixing rules predicts a *flat* g_B(x)
# (it slightly increases with GeO₂), while Niklès 1997 measures g_B dropping
# by ~2× from pure SiO₂ to SMF-28. That gap is the **opto-acoustic overlap**
# effect — the optical and acoustic mode profiles overlap less in doped
# fibres — and is not captured by the bulk formula. Closing it requires the
# acoustic mode solver (`acoustic_mode_terms`, above) AND an overlap integral
# with the optical mode, neither of which is in the bulk formula.
#
# So we anchor `g_B,peak` for pure silica at the measured 5×10⁻¹¹ m/W value
# (Agrawal NLFO Ch. 9, with the same number cited by Boyd, Kobyakov, etc.)
# and apply textbook wavelength scaling `g_B ∝ 1/λ²` plus the bulk-mixing-rule
# part of the doping scaling (ρ·v_A in the denominator). This is the part of
# the textbook formula that *is* unambiguous. For doping-dependent g_B in
# practical fibres, prefer a measured value if you have one.

"Measured peak Brillouin gain coefficient for pure SiO₂ at λ = 1550 nm.
Source: Agrawal NLFO 6e Ch. 9 (also Boyd 4e, Kobyakov 2010). Material datum, not a fit."
const BRIL_G_PEAK_SiO2_1550 = 5.0e-11    # m/W

"""
    g_B_peak(λ_pump, n_eff, p12, m_GeO2; Δν_FWHM=BRIL_DFREQ_HZ,
              λ_ref=1550e-9, n_eff_ref=1.4447, p12_ref=0.270, Δν_ref=BRIL_DFREQ_HZ) → m/W

Peak Brillouin gain coefficient anchored at the measured pure-silica value
`BRIL_G_PEAK_SiO2_1550` and scaled by the unambiguous parts of Agrawal NLFO
6e Eq. (9.1.7):

    g_B / g_B_ref  =  (n_eff/n_ref)⁷ · (p₁₂/p₁₂_ref)² · (λ_ref/λ_pump)²
                       · (ρ_ref · v_A_ref) / (ρ · v_A)
                       · (Δν_ref / Δν_FWHM)

Where `ρ` and `v_A` come from the C₁₁/ρ mixing rules; everything on the right
is textbook material data with no fitted coefficients.

**Caveat — opto-acoustic overlap is not modelled.** The measured g_B(SMF-28)
≈ 2.2×10⁻¹¹ m/W vs the formula's prediction of ~5×10⁻¹¹ m/W reflects the
opto-acoustic overlap, which falls off with doping for reasons not captured
in this expression. For high-doping fibres, pass a measured g_B directly.
"""
function g_B_peak(λ_pump::Real, n_eff::Real, p12::Real, m_GeO2::Real=0.0;
                   Δν_FWHM::Real=BRIL_DFREQ_HZ,
                   λ_ref::Real=1550e-9,
                   n_eff_ref::Real=1.4447,
                   p12_ref::Real=0.270,
                   Δν_ref::Real=BRIL_DFREQ_HZ)
    ρ_0  = glass_density(0.0)
    v_0  = acoustic_velocity_bulk(0.0)
    ρ_x  = glass_density(m_GeO2)
    v_x  = acoustic_velocity_bulk(m_GeO2)
    scale = (n_eff   / n_eff_ref)^7 *
            (p12     / p12_ref)^2  *
            (λ_ref   / λ_pump)^2   *
            (ρ_0 * v_0) / (ρ_x * v_x) *
            (Δν_ref / Δν_FWHM)
    return BRIL_G_PEAK_SiO2_1550 * scale
end

# Pure-silica reference value at λ = 1550 nm. Returns 5×10⁻¹¹ m/W by construction.
const BRIL_G_PEAK = g_B_peak(1550e-9, 1.4447, 0.270, 0.0)

# ═══════════════════════════════════════════════════════════════════════════
# Guided L01 acoustic mode solver  (Kobyakov 2010 §3, Shibata 1989)
# ═══════════════════════════════════════════════════════════════════════════
#
# The acoustic field in a step-index fibre satisfies a scalar Helmholtz
# equation analogous to the optical LP01 mode. In the weakly-guiding limit
# the radial eigenvalue equation reduces to
#
#     U · J₁(U) / J₀(U)  =  W · K₁(W) / K₀(W),     U² + W² = V_a²,
#
# the same Bessel-function structure that governs LP01 — see the prof's
# `waveguide_factor` in `fiber-cross-section.jl`. We re-use the Marcuse-style
# approximation here for the acoustic mode.
#
# For SMF-28 the acoustic V-number V_a is large (~10–15) at telecom
# wavelengths, so the L01 mode is tightly confined to the core and the
# effective velocity differs from the core's bulk v_A by ~0.1 %. The
# Niklès-vs-model gap at 1319 nm is therefore NOT primarily a guided-mode
# correction — it is most likely an error in the reference values.
# Implementing this solver settles the question definitively.

"""
    acoustic_v_number(a, Ω_B, v_core, v_cladding) → V_a (dimensionless)

Acoustic V-number for the L01 mode of a step-index fibre, direct analogue
of the optical V = (2π a / λ) · NA. Requires `v_core < v_cladding` for the
mode to be guided in the core.
"""
function acoustic_v_number(a::Real, Ω_B::Real, v_core::Real, v_cladding::Real)
    v_core < v_cladding ||
        throw(ArgumentError(
            "Guided acoustic L01 mode requires v_core < v_cladding; got " *
            "v_core=$v_core m/s, v_cladding=$v_cladding m/s"))
    return (Ω_B * a / v_core) * sqrt(1 - (v_core / v_cladding)^2)
end

"""
    acoustic_waveguide_factor(V_a) → U_a

Marcuse-style approximation to the radial eigenvalue `U` of the L01 acoustic
mode, identical in form to the optical LP01 `waveguide_factor` in the prof's
`fiber-cross-section.jl`:

    U(V) ≈ (1 + √2)·V / (1 + (4 + V⁴)^(1/4))

For `V_a ≫ 1` this saturates at `U_a → 1 + √2 ≈ 2.414` (tightly confined).
"""
function acoustic_waveguide_factor(V_a::Real)
    α = 1 + sqrt(2)
    return α * V_a / (1 + (4 + V_a^4)^(1/4))
end

"""
    acoustic_mode_terms(xs::FiberCrossSection, λ_pump, T_K; max_iter=20, rtol=1e-10)
    acoustic_mode_terms(fiber::Fiber, λ_pump; T_K=fiber.T_ref_K, ...)

Self-consistent solution for the guided L01 acoustic mode in a SiO₂–GeO₂
step-index fibre under the backward-Brillouin Bragg condition.

Method (fixed-point iteration, converges in 1–3 steps for SMF-28):

    v_eff ← v_core (initial guess)
    repeat:
        Ω_B  = 4π · n_eff · v_eff / λ_pump      (Bragg)
        V_a  = (Ω_B · a / v_core) · √(1 − v_core²/v_cladding²)
        U_a  = (1+√2)·V_a / (1 + (4+V_a⁴)^(1/4))
        β_a  = √((Ω_B/v_core)² − (U_a/a)²)
        v_eff_new = Ω_B / β_a
    until |v_eff_new − v_eff| / v_eff < rtol

Returns a NamedTuple:
- `v_eff`       — effective phase velocity of the guided acoustic mode (m/s)
- `v_core`      — bulk acoustic velocity in the core (m/s)
- `v_cladding`  — bulk acoustic velocity in the cladding (m/s)
- `V_a`         — final acoustic V-number
- `U_a`         — radial eigenvalue (Marcuse approximation)
- `nu_B`        — Brillouin frequency shift at this wavelength (Hz)
- `iterations`  — number of fixed-point steps taken

References:
- Kobyakov, Sauer & Chowdhury, *Adv. Opt. Photon.* 2, 1 (2010), §3.
- Shibata, Azuma & Mochizuki, *Electron. Lett.* 25, 1404 (1989).

Restriction: currently requires both `core_material` and `cladding_material`
to be `GermaniaSilicaGlass` (so the mixing-rule v_A formulas apply). Extending
to F-doped cladding is a small generalisation but requires elastic constants
for fluorosilicate glasses that aren't tabulated in `material-properties.jl`.
"""
function acoustic_mode_terms(xs::FiberCrossSection, λ_pump::Real, T_K::Real;
                              max_iter::Int=20, rtol::Real=1e-10)
    (xs.core_material isa GermaniaSilicaGlass) &&
    (xs.cladding_material isa GermaniaSilicaGlass) ||
        throw(ArgumentError(
            "acoustic_mode_terms currently requires GermaniaSilicaGlass for " *
            "both core and cladding; got core=$(typeof(xs.core_material)), " *
            "cladding=$(typeof(xs.cladding_material))."))

    a    = core_radius(xs)
    x_co = xs.core_material.x_ge
    x_cl = xs.cladding_material.x_ge
    v_co = acoustic_velocity_bulk(x_co)
    v_cl = acoustic_velocity_bulk(x_cl)

    n_eff = effective_mode_index(xs, λ_pump, T_K)

    v_eff = float(v_co)
    iter  = 0
    for k in 1:max_iter
        iter = k
        Ω_B = 4π * n_eff * v_eff / λ_pump
        V_a = acoustic_v_number(a, Ω_B, v_co, v_cl)
        U_a = acoustic_waveguide_factor(V_a)
        # β_a from the dispersion relation (Ω_B/v_co)² = β_a² + (U_a/a)²
        radicand = (Ω_B / v_co)^2 - (U_a / a)^2
        radicand > 0 ||
            throw(ErrorException(
                "Acoustic mode is below cutoff at λ_pump=$λ_pump m: " *
                "Ω_B/v_co < U_a/a. This shouldn't happen for SMF-like fibres."))
        β_a = sqrt(radicand)
        v_eff_new = Ω_B / β_a
        if abs(v_eff_new - v_eff) / v_eff < rtol
            v_eff = v_eff_new
            break
        end
        v_eff = v_eff_new
    end

    Ω_B = 4π * n_eff * v_eff / λ_pump
    V_a = acoustic_v_number(a, Ω_B, v_co, v_cl)
    U_a = acoustic_waveguide_factor(V_a)

    return (
        v_eff      = v_eff,
        v_core     = v_co,
        v_cladding = v_cl,
        V_a        = V_a,
        U_a        = U_a,
        nu_B       = Ω_B / (2π),
        iterations = iter,
    )
end

acoustic_mode_terms(fiber::Fiber, λ_pump::Real;
                     T_K=fiber.T_ref_K, kwargs...) =
    acoustic_mode_terms(fiber.cross_section, λ_pump, T_K; kwargs...)


# ── SMF-28 fiber-loss model ─────────────────────────────────────────────────
#
# EMPIRICAL FIT — not first-principles. Two coefficients (Rayleigh prefactor
# and IR absorption floor) are calibrated to the Corning SMF-28e+ datasheet.
# The model is invoked only as a fallback when callers of `brillouin_threshold`,
# `spbs_photon_rate_density`, or `spbs_noise_in_channel` do not supply `α`
# explicitly. Whenever the fallback fires, a one-shot warning is emitted naming
# the model, the formula, and the values at common telecom wavelengths so that
# the user is never silently dependent on the fitted coefficients.
#
# Reference: Corning SMF-28e+ Optical Fiber Product Information, attenuation
# specs at 1310, 1383, 1490, 1550, 1625 nm. Coefficients fit by inspection.

const _LOSS_A_RAY_BRIL = 0.78    # dB/(km·μm⁴), Rayleigh-scattering prefactor — EMPIRICAL
const _LOSS_FLOOR_BRIL = 0.065   # dB/km, IR absorption floor                  — EMPIRICAL

"""
    _fiber_loss_m(λ) → α in 1/m

Empirical attenuation coefficient for Corning SMF-28-class telecom fibre:

    α(λ) [dB/km]  =  0.78 / (λ[μm])⁴  +  0.065
    α(λ) [1/m]    =  α[dB/km] / (10·log₁₀(e)) / 1000

**Not a first-principles model.** Emits a one-shot warning (`maxlog=1`) on
first invocation per Julia session, naming the empirical coefficients and the
predicted α at typical telecom wavelengths. To suppress the warning entirely,
pass `α` explicitly to whichever function would otherwise default to this.
"""
function _fiber_loss_m(λ::Real)
    @warn """
    SMF-28 EMPIRICAL fiber-loss model invoked (no α supplied by caller).
      α(λ) [dB/km] = 0.78 / (λ[μm])⁴ + 0.065  (fit to Corning SMF-28e+ datasheet)
      At 1310 nm: α ≈ 0.33 dB/km
      At 1383 nm: α ≈ 0.28 dB/km
      At 1490 nm: α ≈ 0.23 dB/km
      At 1550 nm: α ≈ 0.20 dB/km
      At 1625 nm: α ≈ 0.18 dB/km
    Coefficients (0.78 dB·km⁻¹·μm⁴ Rayleigh, 0.065 dB·km⁻¹ IR floor) are an
    empirical 2-parameter fit, not derived from first principles. Pass `α`
    explicitly to the downstream call to suppress this warning and to use a
    measured / first-principles attenuation value instead.""" maxlog=1
    lam_um = λ * 1e6
    dB_km  = _LOSS_A_RAY_BRIL / lam_um^4 + _LOSS_FLOOR_BRIL
    return dB_km / (10 * log10(ℯ)) / 1e3
end

# Array overload broadcasts over the scalar version; the first broadcast
# invocation trips the one-shot warning, subsequent invocations are silent.
_fiber_loss_m(λ::AbstractArray) = _fiber_loss_m.(λ)

# ═══════════════════════════════════════════════════════════════════════════
# Brillouin frequency shift
# ═══════════════════════════════════════════════════════════════════════════

"""
    brillouin_freq_shift(λ_pump; n_eff=1.4447, v_acoustic=nothing, m_GeO2=0.0) → Hz

Backward-SBS Bragg condition: `ν_B = 2 n_eff v_A / λ_pump`.

When `v_acoustic === nothing` (default), it is computed from the
`acoustic_velocity_bulk(m_GeO2)` mixing rule (textbook, no fits) — this is
the bulk material velocity, which is a few-percent over-estimate for guided
fibres. For a self-consistent guided-mode treatment, use the
`Fiber`-method below which calls `acoustic_mode_terms`.
"""
function brillouin_freq_shift(λ_pump; n_eff::Real=1.4447,
                                v_acoustic::Union{Real,Nothing}=nothing,
                                m_GeO2::Real=0.0)
    v_A = v_acoustic === nothing ? acoustic_velocity_bulk(m_GeO2) : v_acoustic
    return 2 .* n_eff .* v_A ./ λ_pump
end

"""
    brillouin_freq_shift(fiber::Fiber, λ_pump; T_K=fiber.T_ref_K, guided::Bool=true) → Hz

Fiber-aware Brillouin shift. With `guided=true` (default), solves the L01
acoustic eigenvalue equation via `acoustic_mode_terms` and uses the
guided-mode `v_eff` in the Bragg condition. With `guided=false`, falls back
to the bulk mixing-rule `v_A` (faster; for SMF-28 the answer differs by
≲ 0.1 % from the guided-mode result because V_a ≫ 1).
"""
function brillouin_freq_shift(fiber::Fiber, λ_pump::Real;
                                T_K=fiber.T_ref_K, guided::Bool=true)
    if guided
        return acoustic_mode_terms(fiber, λ_pump; T_K=T_K).nu_B
    else
        n_eff = effective_mode_index(fiber.cross_section, λ_pump, T_K)
        x_ge  = fiber.cross_section.core_material.x_ge
        v_A   = acoustic_velocity_bulk(x_ge)
        return 2 * n_eff * v_A / λ_pump
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# Peak gain — textbook (Agrawal Eq. 9.1.7)
# ═══════════════════════════════════════════════════════════════════════════

"""
    g_B_peak_GeO2(m_GeO2=0.0; λ_pump=1550e-9, n_eff=1.4447, p12=0.270,
                  Δν_FWHM=BRIL_DFREQ_HZ) → m/W

Textbook Brillouin peak gain (Agrawal NLFO 6e Eq. 9.1.7) for a SiO₂–GeO₂
binary glass at molar fraction `m_GeO2`. No empirical doping fit — the
dependence on `m_GeO2` comes through the mixing rules in `glass_density`
and `acoustic_velocity_bulk`.

Defaults for `n_eff`, `p12`, and `λ_pump` correspond to typical SMF-28 values
at 1550 nm. For other fibres, pass measured / Sellmeier-derived values.
"""
function g_B_peak_GeO2(m_GeO2::Real=0.0;
                        λ_pump::Real=1550e-9,
                        n_eff::Real=1.4447,
                        p12::Real=0.270,
                        Δν_FWHM::Real=BRIL_DFREQ_HZ)
    return g_B_peak(λ_pump, n_eff, p12, m_GeO2; Δν_FWHM=Δν_FWHM)
end

"""
    g_B_peak_GeO2(fiber::Fiber; λ_pump=1550e-9, Δν_FWHM=BRIL_DFREQ_HZ) → m/W

Compute the textbook Brillouin peak gain for `fiber`, reading the GeO₂ molar
fraction, the photoelastic constant `p₁₂`, and `n_eff` straight off the
cross section (no fitted values).
"""
function g_B_peak_GeO2(fiber::Fiber; λ_pump::Real=1550e-9,
                        Δν_FWHM::Real=BRIL_DFREQ_HZ)
    core = fiber.cross_section.core_material
    core isa GermaniaSilicaGlass ||
        throw(ArgumentError("g_B_peak_GeO2(::Fiber) requires a GermaniaSilicaGlass core; got $(typeof(core))"))
    _, p12 = photoelastic_constants(core, fiber.T_ref_K)
    n_eff  = effective_mode_index(fiber.cross_section, λ_pump, fiber.T_ref_K)
    return g_B_peak(λ_pump, n_eff, p12, core.x_ge; Δν_FWHM=Δν_FWHM)
end

# ═══════════════════════════════════════════════════════════════════════════
# Lorentzian gain spectrum
# ═══════════════════════════════════════════════════════════════════════════

"`g_B(Ω) = g_peak · (Γ/2)² / ((Ω − Ω_B)² + (Γ/2)²)`."
function g_B_lorentzian(Ω, Ω_B; g_B_peak::Real=BRIL_G_PEAK,
                          Γ_B::Real=BRIL_GAMMA_RAD)
    hwhm = Γ_B / 2
    return @. g_B_peak * hwhm^2 / ((Ω - Ω_B)^2 + hwhm^2)
end

# ═══════════════════════════════════════════════════════════════════════════
# Effective length (backward geometry)
# ═══════════════════════════════════════════════════════════════════════════

"`L_eff^back = (1 − exp[−(α_p + α_s) L]) / (α_p + α_s)`."
function effective_length_backward(L, α_pump, α_signal=nothing)
    α_p = float.(α_pump)
    α_s = α_signal === nothing ? α_p : float.(α_signal)
    α_sum = α_p .+ α_s
    safe = ifelse.(α_sum .> 1e-12, α_sum, one.(α_sum))
    Leff = -expm1.(-safe .* L) ./ safe
    return ifelse.(α_sum .> 1e-12, Leff, float(L))
end

# ═══════════════════════════════════════════════════════════════════════════
# Thermal phonon occupancy (Brillouin band — same physics, separate name to
# avoid signature collision with the Raman version that uses different unit
# conventions internally)
# ═══════════════════════════════════════════════════════════════════════════

"""
    brillouin_thermal_phonon_number(Ω, T_K) -> dimensionless

Bose-Einstein occupancy at the Brillouin band. At ν_B ≈ 11 GHz and 300 K,
`n_th ≈ 561 ≫ 1`.
"""
function brillouin_thermal_phonon_number(Ω::Real, T_K::Real)
    x = _HBAR_BRIL * abs(Ω) / (_KB_BRIL * T_K)
    return x > 500.0 ? 0.0 : 1.0 / expm1(x)
end
brillouin_thermal_phonon_number(Ω::AbstractArray, T_K::Real) =
    brillouin_thermal_phonon_number.(Ω, T_K)

# ═══════════════════════════════════════════════════════════════════════════
# SBS threshold
# ═══════════════════════════════════════════════════════════════════════════

"""
    brillouin_threshold(A_eff, L; α=nothing, λ_pump=1550e-9,
                        g_B_peak=BRIL_G_PEAK, G_th=BRIL_G_THRESH)
    brillouin_threshold(fiber::Fiber, λ_pump; T_K=fiber.T_ref_K, ...)

Smith-1972 SBS threshold `P_th = G_th · A_eff / (g_B · L_eff^back)`. Returns a
NamedTuple with the threshold and intermediate quantities.
"""
function brillouin_threshold(A_eff, L; α=nothing, λ_pump::Real=1550e-9,
                              g_B_peak::Real=BRIL_G_PEAK,
                              G_th::Real=BRIL_G_THRESH)
    α_use = α === nothing ? float(_fiber_loss_m(λ_pump)) : float(α)
    Leff  = float(effective_length_backward(L, α_use, α_use))
    P_th  = G_th * A_eff / (g_B_peak * Leff)
    return (
        P_threshold_W    = P_th,
        P_threshold_mW   = P_th * 1e3,
        G_B_at_threshold = G_th,
        L_eff_back_m     = Leff,
        g_B_peak         = g_B_peak,
        alpha_m          = α_use,
        formula = @sprintf("P_th = %g * %.1f um^2 / (%.1e m/W * %.2f km) = %.2f mW",
                            G_th, A_eff*1e12, g_B_peak, Leff/1e3, P_th*1e3),
    )
end

function brillouin_threshold(fiber::Fiber, λ_pump; T_K=fiber.T_ref_K,
                              length_m=nothing, α=nothing, G_th::Real=BRIL_G_THRESH)
    A_eff = effective_mode_area(fiber.cross_section, λ_pump, T_K)
    L     = length_m === nothing ? fiber.s_end - fiber.s_start : length_m
    g_B   = fiber.cross_section.core_material isa GermaniaSilicaGlass ?
            g_B_peak_GeO2(fiber) : BRIL_G_PEAK
    return brillouin_threshold(A_eff, L; α=α, λ_pump=λ_pump, g_B_peak=g_B, G_th=G_th)
end

"""
    check_sbs_threshold(P_pump, A_eff, L; λ_pump=1550e-9, threshold_fraction=0.80, ...)
    check_sbs_threshold(P_pump, fiber::Fiber, λ_pump; T_K=fiber.T_ref_K, ...)

Issue a warning if `P_pump > threshold_fraction · P_th`.
"""
function check_sbs_threshold(P_pump, A_eff, L; α=nothing, λ_pump::Real=1550e-9,
                              g_B_peak::Real=BRIL_G_PEAK,
                              threshold_fraction::Real=0.80)
    thresh = brillouin_threshold(A_eff, L; α=α, λ_pump=λ_pump, g_B_peak=g_B_peak)
    P_th  = thresh.P_threshold_W
    G_B   = g_B_peak * P_pump * thresh.L_eff_back_m / A_eff
    frac  = P_pump / P_th
    above = frac > 1
    valid = frac < threshold_fraction
    if !valid
        @warn @sprintf("Pump (%.2f mW) is %.0f%% of SBS threshold (%.2f mW).",
                       P_pump*1e3, frac*100, P_th*1e3)
    end
    return (
        P_threshold_W         = P_th,
        P_pump_W              = P_pump,
        fraction_of_threshold = frac,
        sbs_gain_parameter    = G_B,
        above_threshold       = above,
        valid                 = valid,
    )
end

function check_sbs_threshold(P_pump, fiber::Fiber, λ_pump; T_K=fiber.T_ref_K,
                              length_m=nothing, kwargs...)
    A_eff = effective_mode_area(fiber.cross_section, λ_pump, T_K)
    L     = length_m === nothing ? fiber.s_end - fiber.s_start : length_m
    g_B   = fiber.cross_section.core_material isa GermaniaSilicaGlass ?
            g_B_peak_GeO2(fiber) : BRIL_G_PEAK
    return check_sbs_threshold(P_pump, A_eff, L; λ_pump=λ_pump, g_B_peak=g_B, kwargs...)
end

# ═══════════════════════════════════════════════════════════════════════════
# Spontaneous Brillouin photon-rate density / channel noise
# ═══════════════════════════════════════════════════════════════════════════

"""
    spbs_photon_rate_density(Ω, λ_pump, P_pump, L, A_eff, T_K; ...) -> dṄ/dΩ/dt
"""
function spbs_photon_rate_density(Ω, λ_pump, P_pump, L, A_eff, T_K;
                                    n_eff::Real=1.4447,
                                    v_acoustic::Real=BRIL_VA_SIO2,
                                    m_GeO2::Real=0.0,
                                    g_B_peak::Real=BRIL_G_PEAK,
                                    Γ_B::Real=BRIL_GAMMA_RAD,
                                    α_pump=nothing, α_signal=nothing)
    ν_B = brillouin_freq_shift(λ_pump; n_eff=n_eff, v_acoustic=v_acoustic, m_GeO2=m_GeO2)
    Ω_B = 2π * ν_B
    α_p = α_pump   === nothing ? float(_fiber_loss_m(λ_pump)) : float(α_pump)
    α_s = α_signal === nothing ? α_p                            : float(α_signal)
    Leff = effective_length_backward(L, α_p, α_s)
    gB   = g_B_lorentzian(Ω, Ω_B; g_B_peak=g_B_peak, Γ_B=Γ_B)
    nth  = brillouin_thermal_phonon_number(Ω, T_K)
    return @. (gB / A_eff) * P_pump * Leff * nth / (2π)
end

"""
    spbs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, A_eff, T_K; ...)
    spbs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                          T_K=fiber.T_ref_K, length_m=nothing, kwargs...)

Total backward spontaneous Brillouin photon rate into a rectangular channel
of width `Δλ` centred at `λ_channel`. Returns a NamedTuple of diagnostics
(rate, power, νB, Stokes wavelength, near-resonance flag, threshold info).

The fiber-aware form derives `A_eff`, `n_eff`, `m_GeO2`, and `g_B_peak` from
the cross section.
"""
function spbs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, A_eff, T_K;
                                 N_points::Int=501,
                                 n_eff::Real=1.4447,
                                 v_acoustic::Real=BRIL_VA_SIO2,
                                 m_GeO2::Real=0.0,
                                 g_B_peak::Real=BRIL_G_PEAK,
                                 Γ_B::Real=BRIL_GAMMA_RAD,
                                 α_pump=nothing, α_signal=nothing)
    isodd(N_points) || (N_points += 1)
    ω_pump = 2π * SPEED_OF_LIGHT_M_PER_S / λ_pump
    ν_B    = float(brillouin_freq_shift(λ_pump; n_eff=n_eff, v_acoustic=v_acoustic, m_GeO2=m_GeO2))
    Ω_B    = 2π * ν_B
    λ_S    = SPEED_OF_LIGHT_M_PER_S / (SPEED_OF_LIGHT_M_PER_S / λ_pump - ν_B)

    ν_offset = abs(SPEED_OF_LIGHT_M_PER_S / λ_channel - (SPEED_OF_LIGHT_M_PER_S / λ_pump - ν_B))
    near = ν_offset < 10 * BRIL_DFREQ_HZ

    lam_arr = collect(range(λ_channel - Δλ/2, λ_channel + Δλ/2; length=N_points))
    om_arr  = 2π .* SPEED_OF_LIGHT_M_PER_S ./ lam_arr
    Ω_arr   = abs.(ω_pump .- om_arr)

    α_p = α_pump   === nothing ? float(_fiber_loss_m(λ_pump)) : float(α_pump)
    α_s = α_signal === nothing ? _fiber_loss_m(lam_arr) : float(α_signal) .* ones(length(lam_arr))

    Leff_arr = effective_length_backward(L, α_p, α_s)
    gB_arr   = g_B_lorentzian(Ω_arr, Ω_B; g_B_peak=g_B_peak, Γ_B=Γ_B)
    nth_arr  = brillouin_thermal_phonon_number(Ω_arr, T_K)

    dNdOm = @. (gB_arr / A_eff) * P_pump * Leff_arr * nth_arr / (2π)
    rate  = float(_simpson(dNdOm .* (2π .* SPEED_OF_LIGHT_M_PER_S ./ lam_arr.^2), lam_arr))
    ω_ch  = 2π * SPEED_OF_LIGHT_M_PER_S / λ_channel

    thresh = brillouin_threshold(A_eff, L; α=α_p, λ_pump=λ_pump, g_B_peak=g_B_peak)
    return (
        backward_photon_rate       = rate,
        backward_power_W           = rate * _HBAR_BRIL * ω_ch,
        Omega_B_rad                = Ω_B,
        nu_B_Hz                    = ν_B,
        delta_lambda_B_pm          = abs(λ_S - λ_pump) * 1e12,
        lambda_Stokes_nm           = λ_S * 1e9,
        channel_is_near_Brillouin  = near,
        P_threshold_mW             = thresh.P_threshold_mW,
        pump_fraction_of_threshold = P_pump / thresh.P_threshold_W,
        lambda_pump                = λ_pump,
        lambda_channel             = λ_channel,
        delta_lambda               = Δλ,
    )
end

function spbs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                                 T_K=fiber.T_ref_K,
                                 length_m=nothing,
                                 kwargs...)
    xs    = fiber.cross_section
    A_eff = effective_mode_area(xs, λ_pump, T_K)
    L     = length_m === nothing ? fiber.s_end - fiber.s_start : length_m
    n_eff = effective_mode_index(xs, λ_pump, T_K)
    m_Ge  = xs.core_material isa GermaniaSilicaGlass ? xs.core_material.x_ge : 0.0
    g_B   = xs.core_material isa GermaniaSilicaGlass ? g_B_peak_GeO2(m_Ge) : BRIL_G_PEAK
    return spbs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, A_eff, T_K;
                                   n_eff=n_eff, m_GeO2=m_Ge, g_B_peak=g_B,
                                   kwargs...)
end

# ═══════════════════════════════════════════════════════════════════════════
# Brillouin gain profile (utility for plotting)
# ═══════════════════════════════════════════════════════════════════════════

"""
    brillouin_gain_profile(λ_pump; freq_range_GHz=(-100,100), N=4096, ...)

Compute the Brillouin Lorentzian gain profile around the pump. NamedTuple with
`freq_offset_GHz`, `wavelength_m`, `g_B_mW`, `nu_B_GHz`, `g_B_peak`, `Gamma_B_MHz`.
"""
function brillouin_gain_profile(λ_pump; freq_range_GHz::Tuple=(-100, 100), N::Int=4096,
                                  n_eff::Real=1.4447,
                                  v_acoustic::Real=BRIL_VA_SIO2,
                                  m_GeO2::Real=0.0,
                                  g_B_peak::Real=BRIL_G_PEAK,
                                  Γ_B::Real=BRIL_GAMMA_RAD)
    ν_B = float(brillouin_freq_shift(λ_pump; n_eff=n_eff, v_acoustic=v_acoustic, m_GeO2=m_GeO2))
    Ω_B = 2π * ν_B
    ν_pump = SPEED_OF_LIGHT_M_PER_S / λ_pump
    ν_arr  = collect(range(ν_pump + freq_range_GHz[1] * 1e9,
                             ν_pump + freq_range_GHz[2] * 1e9; length=N))
    Ω_arr  = 2π .* abs.(ν_arr .- ν_pump)
    gB_arr = g_B_lorentzian(Ω_arr, Ω_B; g_B_peak=g_B_peak, Γ_B=Γ_B)
    return (
        freq_offset_GHz = (ν_arr .- ν_pump) ./ 1e9,
        wavelength_m    = SPEED_OF_LIGHT_M_PER_S ./ ν_arr,
        g_B_mW          = gB_arr,
        nu_B_GHz        = ν_B / 1e9,
        g_B_peak        = g_B_peak,
        Gamma_B_MHz     = Γ_B / (2π * 1e6),
    )
end
