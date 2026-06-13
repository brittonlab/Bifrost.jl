"""
silica_raman.jl — Spontaneous Raman scattering noise in **fused silica fibre**.

This file is included into the `Nonlinear` submodule of `Bifrost` by
`Bifrost.jl`. The submodule's `using ..MaterialProperties`, `..FiberCS`, and
`..FiberPath` lines bring the types and helpers this file relies on
(`Fiber`, `FiberCrossSection`, `nonlinear_coefficient`, `SPEED_OF_LIGHT_M_PER_S`)
into scope.

# Material scope

All Raman parameters here describe pure fused silica: the Blow-Wood τ₁, τ₂
time constants; the Lin & Agrawal (2006) 50-point tabulated spectrum; the
Hollenbeck-Cantrell (2002) 13-mode Gaussian-damped oscillator parameters;
and `fR = 0.18`. GeO₂ doping is not modelled — for telecom-grade core
concentrations this is typically a small correction at the 13 THz peak but
can be significant at high doping. 

# Configuration

Model choice and fractional Raman contribution are bundled into a single
[`RamanConfig`](@ref) value, passed per call. There is **no module-level
default** — pick the model explicitly at each call site so the choice is
visible in the code that uses it.

  `:bw`         Blow-Wood (1989) two-parameter damped-oscillator model.
                Fast, analytic closed form. Good near 13 THz; overestimates
                the tail beyond ~20 THz and misses the 15 THz shoulder.

  `:tabulated`  Lin & Agrawal (2006) 50-point KK-normalised cubic spline.
                Accurate over 0-25 THz; correct tail fall-off.

  `:hc`         Hollenbeck & Cantrell (2002) 13-mode Gaussian-damped oscillator
                model. Captures the 15.2 THz shoulder, ~160 fs phase reversal,
                and structure up to ~1 ps. Precomputed at module load via a
                vectorised sine transform; subsequent calls use a cubic-spline
                interpolator.

# Fiber-aware top-level API

```julia
sprs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                      T_K=fiber.T_ref_K, config=RamanConfig(), ...)
```

derives `γ = (2π/λ)·n₂/A_eff` from `fiber.cross_section` and the
fiber length from `fiber.s_end - fiber.s_start`, then forwards to the raw
`(L, γ, T_K)` form, which is also kept available for reuse outside a `Fiber`
context.

# References
- K. J. Blow & D. Wood, IEEE J. Quantum Electron. 25, 2665 (1989).
- G. P. Agrawal, *Nonlinear Fiber Optics*, 6th ed. (2019).
- Q. Lin & G. P. Agrawal, *Opt. Lett.* 31, 3086 (2006).
- D. Hollenbeck & C. D. Cantrell, *J. Opt. Soc. Am. B* 19, 2886 (2002).
"""

using DataInterpolations
using FFTW
using Printf

# ── Physical constants (private) ──────────────────────────────────────────────
# `SPEED_OF_LIGHT_M_PER_S` comes from material-properties.jl and is in scope.
const _HBAR_RAMAN = 1.0545718e-34   # J·s
const _KB_RAMAN   = 1.380649e-23    # J/K
const _C_CM       = SPEED_OF_LIGHT_M_PER_S * 100   # cm/s

# ── Blow-Wood constants (silica) ──────────────────────────────────────────────
const RAMAN_TAU1_SI = 12.2e-15      # s
const RAMAN_TAU2_SI = 32.0e-15      # s
const RAMAN_FR_SI   = 0.18          # fractional Raman contribution (silica)
const RAMAN_OMEGA_R = 2π * 13.2e12  # rad/s, Raman peak

# ── Model + fR config (replaces the previous module-level Ref) ────────────────
const VALID_RAMAN_MODELS = (:bw, :tabulated, :hc)

"""
    RamanConfig(; model=:bw, fR=RAMAN_FR_SI)
    RamanConfig(model::Symbol, fR::Real)

Bundle of the two parameters that specify a silica Raman response: the model
symbol and the fractional Raman contribution `fR`. The inner constructor
validates both — `model ∈ (:bw, :tabulated, :hc)` and `0 ≤ fR ≤ 1` — so
downstream functions can trust the values without re-checking.

# Examples
```julia
g_R(Ω, γ; config=RamanConfig())                       # defaults: :bw, fR=0.18
g_R(Ω, γ; config=RamanConfig(model=:hc))              # HC model, default fR
g_R(Ω, γ; config=RamanConfig(:tabulated, 0.20))       # positional form
```
"""
struct RamanConfig
    model::Symbol
    fR::Real
    function RamanConfig(model::Symbol, fR::Real)
        model ∈ VALID_RAMAN_MODELS ||
            throw(ArgumentError(
                "model must be one of $VALID_RAMAN_MODELS, got $(repr(model))"))
        0 ≤ fR ≤ 1 ||
            throw(ArgumentError(
                "fR (fractional Raman contribution) must be in [0, 1], got $fR"))
        return new(model, fR)
    end
end
RamanConfig(; model::Symbol=:bw, fR::Real=RAMAN_FR_SI) = RamanConfig(model, fR)

# ── Numerical helpers (private; reused by silica_brillouin.jl) ────────────────
"""Trapezoidal integration on a possibly non-uniform 1-D grid."""
function _trapz(y::AbstractVector, x::AbstractVector)
    @assert length(y) == length(x)
    s = zero(promote_type(eltype(y), eltype(x)))
    @inbounds for i in 1:length(y)-1
        s += (y[i] + y[i+1]) * (x[i+1] - x[i]) / 2
    end
    return s
end

"""Composite Simpson's rule on a possibly non-uniform 1-D grid (scipy-compatible)."""
function _simpson(y::AbstractVector, x::AbstractVector)
    n = length(y)
    @assert n == length(x)
    n < 2 && return zero(promote_type(eltype(y), eltype(x)))
    n == 2 && return (x[2] - x[1]) * (y[1] + y[2]) / 2

    s = zero(promote_type(eltype(y), eltype(x)))
    i = 1
    @inbounds while i + 2 <= n
        h0 = x[i+1] - x[i]
        h1 = x[i+2] - x[i+1]
        h  = h0 + h1
        s += h/6 * ((2 - h1/h0) * y[i] +
                    (h^2 / (h0*h1)) * y[i+1] +
                    (2 - h0/h1) * y[i+2])
        i += 2
    end
    if i < n
        s += (x[n] - x[n-1]) * (y[n-1] + y[n]) / 2
    end
    return s
end

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Blow-Wood model
# ═══════════════════════════════════════════════════════════════════════════════

"Blow-Wood impulse response hR(t). Causal, normalised: ∫₀^∞ hR(t) dt = 1."
function h_R_time(t::Real)
    τ1, τ2 = RAMAN_TAU1_SI, RAMAN_TAU2_SI
    t > 0 || return 0.0
    prefac = (τ1^2 + τ2^2) / (τ1 * τ2^2)
    return prefac * exp(-t/τ2) * sin(t/τ1)
end
h_R_time(t::AbstractArray) = h_R_time.(t)

"Analytical Fourier transform of `h_R_time`. Returns complex h̃R(Ω)."
function h_R_freq(Ω::Real)
    τ1, τ2 = RAMAN_TAU1_SI, RAMAN_TAU2_SI
    num = (τ1^2 + τ2^2) / (τ1^2 * τ2^2)
    den = (1/τ2 - im*Ω)^2 + (1/τ1)^2
    return num / den
end
h_R_freq(Ω::AbstractArray) = h_R_freq.(Ω)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Hollenbeck-Cantrell (HC) 13-mode model
# ═══════════════════════════════════════════════════════════════════════════════

const _HC_NU_CM = [
    56.25, 100.00, 231.25, 362.50, 463.00, 497.00, 611.50,
    691.67, 793.67, 835.50, 930.00, 1080.00, 1215.00,
]
const _HC_A = [
    1.00, 11.40, 36.67, 67.67, 74.00, 4.50, 6.80,
    4.60, 4.20, 4.50, 2.70, 3.10, 3.00,
]
const _HC_G_FWHM_CM = [
    52.10, 110.42, 175.00, 162.50, 135.33, 24.50, 41.50,
    155.00, 59.50, 64.30, 150.00, 91.00, 160.00,
]
const _HC_L_FWHM_CM = [
    17.37, 38.81, 58.33, 54.17, 45.11, 8.17, 13.83,
    51.67, 19.83, 21.43, 50.00, 30.33, 53.33,
]

const _HC_OMEGA_V = 2π .* _C_CM .* _HC_NU_CM
const _HC_GAMMA   = π  .* _C_CM .* _HC_G_FWHM_CM
const _HC_GAMMA_L = π  .* _C_CM .* _HC_L_FWHM_CM

function _hc_h_R_unnorm(t::AbstractVector)
    tf = float.(t)
    h  = zeros(Float64, length(tf))
    @inbounds for k in eachindex(tf)
        tk = tf[k]
        tk > 0 || continue
        acc = 0.0
        for i in eachindex(_HC_A)
            acc += _HC_A[i] *
                   exp(-_HC_GAMMA_L[i] * tk) *
                   exp(-_HC_GAMMA[i]^2 * tk^2 / 4) *
                   sin(_HC_OMEGA_V[i] * tk)
        end
        h[k] = acc
    end
    return h
end
_hc_h_R_unnorm(t::Real) = _hc_h_R_unnorm([float(t)])[1]

# Precompute HC spline at module load
const _HC_T_MAX = 3.0e-12
const _HC_DT    = 0.5e-15
const _N_T_HC   = Int(round(_HC_T_MAX / _HC_DT))
const _t_hc     = collect(range(0.0; step=_HC_DT, length=_N_T_HC))

const _h_hc_raw  = _hc_h_R_unnorm(_t_hc)
const _HC_NORM   = _trapz(_h_hc_raw, _t_hc)
const _h_hc_norm = _h_hc_raw ./ _HC_NORM

const _HC_OMEGA_MAX  = 2π * 40.0e12
const _HC_N_OMEGA    = 801
const _omega_hc_grid = collect(range(0.0, _HC_OMEGA_MAX; length=_HC_N_OMEGA))

const _im_hR_hc_pos = let
    out = zeros(Float64, _HC_N_OMEGA)
    @inbounds for k in 1:_HC_N_OMEGA
        ω = _omega_hc_grid[k]
        s = 0.0
        prev = _h_hc_norm[1] * sin(ω * _t_hc[1])
        for j in 2:length(_t_hc)
            tj  = _t_hc[j]
            cur = _h_hc_norm[j] * sin(ω * tj)
            s  += (prev + cur) * (tj - _t_hc[j-1]) / 2
            prev = cur
        end
        out[k] = s
    end
    out
end

const _HC_SPLINE = DataInterpolations.CubicSpline(_im_hR_hc_pos, _omega_hc_grid)

"Im[h̃R(Ω)] for the Hollenbeck-Cantrell 13-mode model. Antisymmetric in Ω."
function im_h_R_hc(Ω::Real)
    Ωf    = float(Ω)
    abs_Ω = abs(Ωf)
    abs_Ω > _HC_OMEGA_MAX && return 0.0
    val = _HC_SPLINE(abs_Ω)
    isfinite(val) || (val = 0.0)
    val = max(val, 0.0)
    return Ωf >= 0.0 ? val : -val
end
im_h_R_hc(Ω::AbstractArray) = im_h_R_hc.(Ω)

"HC normalised impulse response hR(t)."
h_R_time_hc(t::Real) = _hc_h_R_unnorm(float(t)) / _HC_NORM
h_R_time_hc(t::AbstractArray) = _hc_h_R_unnorm(t) ./ _HC_NORM

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Lin-Agrawal (2006) tabulated model
# ═══════════════════════════════════════════════════════════════════════════════

const _FREQ_THZ = [
    0.00,  0.50,  1.00,  1.50,  2.00,  2.50,  3.00,  3.50,
    4.00,  4.50,  5.00,  5.50,  6.00,  6.50,  7.00,  7.50,
    8.00,  8.50,  9.00,  9.50, 10.00, 10.50, 11.00, 11.50,
   12.00, 12.50, 13.00, 13.20, 13.50, 14.00, 14.50, 15.00,
   15.20, 15.50, 16.00, 16.50, 17.00, 17.50, 18.00, 18.50,
   19.00, 19.50, 20.00, 21.00, 22.00, 23.00, 24.00, 25.00,
   27.00, 30.00,
]

const _IM_HR_RAW = [
    0.000, 0.005, 0.015, 0.032, 0.055, 0.080, 0.105, 0.138,
    0.175, 0.215, 0.260, 0.305, 0.355, 0.405, 0.455, 0.508,
    0.562, 0.617, 0.672, 0.724, 0.772, 0.818, 0.858, 0.893,
    0.926, 0.960, 0.988, 1.000, 0.990, 0.962, 0.928, 0.895,
    0.900, 0.888, 0.850, 0.790, 0.710, 0.610, 0.495, 0.380,
    0.275, 0.185, 0.118, 0.042, 0.013, 0.004, 0.001, 0.000,
    0.000, 0.000,
]

@assert length(_FREQ_THZ) == length(_IM_HR_RAW) "Tabulated Raman: length mismatch"

const _FREQ_RAD = _FREQ_THZ .* 1e12 .* 2π

const _NORM_C = let
    mask    = _FREQ_RAD .> 0
    om      = _FREQ_RAD[mask]
    integ   = _IM_HR_RAW[mask] ./ om
    kk_int  = (2/π) * _simpson(integ, om)
    1.0 / kk_int
end

const _IM_HR_NORM = _IM_HR_RAW .* _NORM_C
const _TAB_SPLINE = DataInterpolations.CubicSpline(_IM_HR_NORM, _FREQ_RAD)

# Largest |Ω| at which the silica Raman gain is experimentally constrained.
# Pulled from the Lin & Agrawal (2006) tabulation; beyond this point the
# tabulated and HC models return 0 by construction and the BW analytic form
# is unvalidated. Used by `_validate_raman_channel` to warn callers when
# their integration window falls outside the model's domain.
const _RAMAN_OMEGA_MAX = _FREQ_RAD[end]

"Im[h̃R(Ω)] using the Lin & Agrawal (2006) tabulated spectrum."
function im_h_R_tabulated(Ω::Real)
    Ωf      = float(Ω)
    abs_Ω   = abs(Ωf)
    max_rad = _FREQ_RAD[end]
    (abs_Ω > max_rad || abs_Ω == 0.0) && return 0.0
    val = _TAB_SPLINE(abs_Ω)
    isfinite(val) || return 0.0
    val = max(val, 0.0)
    return Ωf >= 0.0 ? val : -val
end
im_h_R_tabulated(Ω::AbstractArray) = im_h_R_tabulated.(Ω)

"""
    h_R_time_tabulated(t; n_omega=2^17)

Time-domain Raman response from the tabulated h̃R(Ω) via inverse Fourier transform.
Slow on first call (large IFFT). For fast eval use `h_R_time` (Blow-Wood).
"""
function h_R_time_tabulated(t::AbstractArray; n_omega::Int=2^17)
    tf = float.(t)
    Om_max  = 40e12 * 2π
    Om_grid = collect(range(0.0, Om_max; length=n_omega ÷ 2 + 1))
    dOm     = Om_grid[2] - Om_grid[1]
    Im_h    = im_h_R_tabulated(Om_grid)

    Om2  = Om_grid .^ 2
    Re_h = zeros(Float64, length(Im_h))
    integrand = similar(Im_h)
    @inbounds for k in eachindex(Om_grid)
        ω2k = Om2[k]
        for j in eachindex(Om_grid)
            integrand[j] = (j == k) ? 0.0 : Om_grid[j] * Im_h[j] / (Om2[j] - ω2k)
        end
        Re_h[k] = (2/π) * _trapz(integrand, Om_grid)
    end
    Re_h[1] = 1.0

    h_tilde = Re_h .+ im .* Im_h
    h_full  = vcat(h_tilde, conj.(reverse(h_tilde[2:end-1])))

    n_full  = length(h_full)
    h_t_fft = real(ifft(ifftshift(h_full))) .* (Om_max / π)
    dt    = 2π / (n_full * dOm)
    t_fft = ifftshift(((0:n_full-1) .- n_full ÷ 2) .* dt)

    perm      = sortperm(t_fft)
    t_fft_s   = t_fft[perm]
    h_t_fft_s = h_t_fft[perm]

    h_out = zeros(Float64, length(tf))
    @inbounds for i in eachindex(tf)
        ti = tf[i]
        (ti <= 0 || ti < t_fft_s[1] || ti > t_fft_s[end]) && continue
        j = searchsortedfirst(t_fft_s, ti)
        if j == 1
            h_out[i] = h_t_fft_s[1]
        elseif j > length(t_fft_s)
            h_out[i] = h_t_fft_s[end]
        else
            t0 = t_fft_s[j-1]; t1 = t_fft_s[j]
            h0 = h_t_fft_s[j-1]; h1 = h_t_fft_s[j]
            h_out[i] = h0 + (h1 - h0) * (ti - t0) / (t1 - t0)
        end
    end
    return h_out
end
h_R_time_tabulated(t::Real; kwargs...) = h_R_time_tabulated([float(t)]; kwargs...)[1]

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Unified g_R dispatcher
# ═══════════════════════════════════════════════════════════════════════════════

"""
    g_R(Ω, γ; config=RamanConfig()) -> W⁻¹m⁻¹

Raman gain coefficient. `Ω` is the angular frequency shift from pump (rad/s);
`γ` is the fiber nonlinear coefficient (W⁻¹m⁻¹). Stokes (`Ω > 0`) gives
positive gain; anti-Stokes (`Ω < 0`) gives negative.

The model choice (`:bw`, `:tabulated`, `:hc`) and the fractional Raman
contribution `fR` are carried in [`RamanConfig`](@ref). Pass it explicitly at
every call site — there is no module-level default.
"""
function g_R(Ω, γ; config::RamanConfig=RamanConfig())
    im_hR = if config.model === :bw
        imag.(h_R_freq(Ω))
    elseif config.model === :tabulated
        im_h_R_tabulated(Ω)
    else  # :hc — validated by RamanConfig
        im_h_R_hc(Ω)
    end
    return 2 .* γ .* config.fR .* im_hR
end

"`g_R` between two wavelengths (m); shift is `2π c (1/λ_pump − 1/λ_signal)`."
function g_R_from_wavelengths(λ_pump, λ_signal, γ; config::RamanConfig=RamanConfig())
    Ω = 2π * SPEED_OF_LIGHT_M_PER_S * (1/λ_pump - 1/λ_signal)
    return float(g_R(Ω, γ; config=config))
end

"Raman gain using the tabulated spectrum (thin wrapper over `g_R`)."
g_R_tabulated(Ω, γ; fR::Real=RAMAN_FR_SI) =
    g_R(Ω, γ; config=RamanConfig(model=:tabulated, fR=fR))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Thermal phonon occupancy
# ═══════════════════════════════════════════════════════════════════════════════

"""
    thermal_phonon_number(Ω, T_K) -> dimensionless

Bose-Einstein occupancy `n_th(|Ω|, T) = 1 / (exp(ħ|Ω|/k_B T) − 1)`.
Clamps for very large arguments to avoid Inf.
"""
function thermal_phonon_number(Ω::Real, T_K::Real)
    x = _HBAR_RAMAN * abs(Ω) / (_KB_RAMAN * T_K)
    return x > 500.0 ? 0.0 : 1.0 / expm1(x)
end
thermal_phonon_number(Ω::AbstractArray, T_K::Real) = thermal_phonon_number.(Ω, T_K)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Spontaneous Raman scattering noise
# ═══════════════════════════════════════════════════════════════════════════════

function _sprs_rate_density_core(gR, P_pump, L, factor; pump_depletion::Bool=false)
    if !pump_depletion
        return @. P_pump * L * gR * factor / (2π)
    else
        x = @. gR * P_pump * L
        corrected_PL = @. ifelse(gR > 1e-25, expm1(x) / gR, P_pump * L)
        return @. corrected_PL * factor / (2π)
    end
end

"""
    sprs_photon_rate_density(Ω, ω_pump, P_pump, L, γ, T_K;
                             sideband=:stokes, pump_depletion=false,
                             config=RamanConfig())

Spectral density `dṄ/dΩ` of spontaneous Raman photons
[photons s⁻¹ (rad/s)⁻¹]. Convention: `Ω = ω_pump − ω_signal`; positive `Ω`
is Stokes.
"""
function sprs_photon_rate_density(Ω, ω_pump, P_pump, L, γ, T_K;
                                    sideband::Symbol=:stokes,
                                    pump_depletion::Bool=false,
                                    config::RamanConfig=RamanConfig())
    Ω_for_gR = sideband === :antistokes ? abs.(Ω) : Ω
    gR  = g_R(Ω_for_gR, γ; config=config)
    nth = thermal_phonon_number(Ω, T_K)
    factor = sideband === :stokes ? (nth .+ 1.0) : nth
    rate = _sprs_rate_density_core(gR, P_pump, L, factor; pump_depletion=pump_depletion)
    return max.(rate, 0.0)
end

"""
    _validate_raman_channel(λ_pump, λ_channel, Δλ, sideband)

Sanity-check the spectral window `[λ_channel − Δλ/2, λ_channel + Δλ/2]`
before integrating the spontaneous-Raman rate density across it. Catches:

  * non-positive wavelengths or width;
  * a width so large that the lower edge would fall to ≤ 0;
  * a channel that overlaps or straddles the pump line (the integrand is
    not defined on that interval);
  * a channel on the wrong side of the pump for the requested sideband
    (Stokes ⇒ longer λ than pump, anti-Stokes ⇒ shorter λ).

Emits an `@warn` (rather than throwing) if any part of the channel maps to
`|Ω| > _RAMAN_OMEGA_MAX`, where the silica gain spectrum is unconstrained
and the integrand evaluates to ~0 for the tabulated/HC models.
"""
function _validate_raman_channel(λ_pump::Real, λ_channel::Real, Δλ::Real,
                                   sideband::Symbol)
    λ_pump    > 0 || throw(ArgumentError("λ_pump must be positive, got $λ_pump"))
    λ_channel > 0 || throw(ArgumentError("λ_channel must be positive, got $λ_channel"))
    Δλ        > 0 || throw(ArgumentError("Δλ must be positive, got $Δλ"))
    Δλ < 2λ_channel || throw(ArgumentError(
        "Δλ = $Δλ m would push the lower channel edge to ≤ 0 (require Δλ < 2·λ_channel = $(2λ_channel) m)"))

    λ_lo, λ_hi = λ_channel - Δλ/2, λ_channel + Δλ/2
    λ_lo ≤ λ_pump ≤ λ_hi && throw(ArgumentError(
        "Channel [$λ_lo, $λ_hi] m overlaps the pump (λ_pump = $λ_pump m); the spontaneous-Raman integrand is not defined on this interval."))

    if sideband === :stokes && λ_hi < λ_pump
        throw(ArgumentError(
            "sideband=:stokes but channel ∈ [$λ_lo, $λ_hi] m sits at shorter wavelength than λ_pump = $λ_pump m. Did you mean :antistokes?"))
    elseif sideband === :antistokes && λ_lo > λ_pump
        throw(ArgumentError(
            "sideband=:antistokes but channel ∈ [$λ_lo, $λ_hi] m sits at longer wavelength than λ_pump = $λ_pump m. Did you mean :stokes?"))
    end

    ω_pump = 2π * SPEED_OF_LIGHT_M_PER_S / λ_pump
    Ω_max  = max(abs(ω_pump - 2π * SPEED_OF_LIGHT_M_PER_S / λ_lo),
                 abs(ω_pump - 2π * SPEED_OF_LIGHT_M_PER_S / λ_hi))
    if Ω_max > _RAMAN_OMEGA_MAX
        @warn @sprintf(
            "Channel reaches %.2f THz from the pump, beyond the silica Raman gain support (%.2f THz). The tabulated and HC integrands are zero past this point; the BW form is unvalidated.",
            Ω_max / (2π * 1e12), _RAMAN_OMEGA_MAX / (2π * 1e12))
    end
    return nothing
end

"""
    sprs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, γ, T_K;
                          sideband=:stokes, pump_depletion=false,
                          n_points=1001, config=RamanConfig())
    sprs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                          T_K=fiber.T_ref_K, length_m=nothing, kwargs...)

Total spontaneous Raman photon rate into a quantum channel [photons / s].

The first form takes raw `(L, γ, T_K)`. The second derives them from a
`Fiber`: `γ = (2π/λ)·n₂/A_eff` from `fiber.cross_section` and
`L = fiber.s_end − fiber.s_start` (override with `length_m`).

The `(λ_pump, λ_channel, Δλ, sideband)` tuple is validated via
[`_validate_raman_channel`](@ref) before integration; see there for the
exact geometric and spectral-support constraints enforced.
"""
function sprs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, γ, T_K;
                                 sideband::Symbol=:stokes,
                                 pump_depletion::Bool=false,
                                 n_points::Int=1001,
                                 config::RamanConfig=RamanConfig())
    _validate_raman_channel(λ_pump, λ_channel, Δλ, sideband)
    ω_pump = 2π * SPEED_OF_LIGHT_M_PER_S / λ_pump
    ω_ch   = 2π * SPEED_OF_LIGHT_M_PER_S / λ_channel
    Δω     = 2π * SPEED_OF_LIGHT_M_PER_S * Δλ / λ_channel^2

    ω_lo = ω_ch - Δω/2
    ω_hi = ω_ch + Δω/2
    ω_s_arr = collect(range(ω_lo, ω_hi; length=n_points))
    Ω_arr   = ω_pump .- ω_s_arr

    rate = sprs_photon_rate_density(Ω_arr, ω_pump, P_pump, L, γ, T_K;
                                      sideband=sideband,
                                      pump_depletion=pump_depletion,
                                      config=config)
    return float(_simpson(rate, ω_s_arr))
end

# `nonlinear_coefficient` was removed from the cross-section layer when it was split into
# `fiber-cross-section/step-index.jl` (main commit c57af29). Reconstruct it locally from the
# accessors that remain so the Nonlinear module stays self-contained:
# γ = (2π/λ)·n₂/A_eff  [W⁻¹m⁻¹].
function _nonlinear_coefficient(xs, λ, T_K)
    n2 = core_nonlinear_refractive_index(xs, λ, T_K)
    A_eff = effective_mode_area(xs, λ, T_K)
    return (2π / λ) * n2 / A_eff
end

function sprs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                                 T_K=fiber.T_ref_K,
                                 length_m=nothing,
                                 kwargs...)
    γ = _nonlinear_coefficient(fiber.cross_section, λ_pump, T_K)
    L = length_m === nothing ? fiber.s_end - fiber.s_start : length_m
    return sprs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, γ, T_K; kwargs...)
end

"""
    check_raman_depletion_validity(P_pump, L, γ; fR=RAMAN_FR_SI, threshold=0.1)

Warn if the pump-depletion correction `[exp(x)−1]/x` (with `x = g_R·P·L`)
exceeds `threshold` (10 % default). Uses the Blow-Wood peak (`:bw` at
`RAMAN_OMEGA_R`) as a cheap upper-bound estimate of `g_R`; the answer is
indicative, not model-exact. Returns the ratio.
"""
function check_raman_depletion_validity(P_pump, L, γ; fR::Real=RAMAN_FR_SI, threshold=0.1)
    gR_peak = g_R(RAMAN_OMEGA_R, γ; config=RamanConfig(model=:bw, fR=fR))
    x = gR_peak * P_pump * L
    ratio = x > 1e-10 ? expm1(x) / x : 1.0
    if ratio - 1.0 > threshold
        @warn @sprintf("Raman pump-depletion correction is %.1f%% (x = gR·P·L = %.3f). Pass pump_depletion=true.",
                       (ratio - 1) * 100, x)
    end
    return ratio
end
