"""
raman.jl — Spontaneous Raman scattering noise.

This file is included into the top-level `BIFROST` module by `BIFROST.jl`,
**not** wrapped in its own module — so it shares scope with the rest of BIFROST
(types like `Fiber`, `FiberCrossSection`, helpers like `nonlinear_coefficient`).

Three Raman response-function models are exposed via the module-level
`RAMAN_MODEL` ref or a per-call `model=` keyword:

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
sprs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump; T_K=fiber.T_ref_K, ...)
```

derives γ via `nonlinear_coefficient(fiber.cross_section, λ_pump, T_K)` and the
fiber length from `fiber.s_end - fiber.s_start`, then forwards to the legacy
raw-parameter form. The raw form is also kept available for reuse outside a
`Fiber` context.

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

# ── Model selector ────────────────────────────────────────────────────────────
"""Active Raman model (`:bw`, `:tabulated`, or `:hc`). Set with `set_raman_model!`."""
const RAMAN_MODEL = Ref(:bw)
const VALID_RAMAN_MODELS = (:bw, :tabulated, :hc)

"""
    set_raman_model!(model::Symbol) -> Symbol

Set the module-level default Raman model.
"""
function set_raman_model!(model::Symbol)
    model ∈ VALID_RAMAN_MODELS ||
        throw(ArgumentError("model must be one of $VALID_RAMAN_MODELS, got $(repr(model))"))
    RAMAN_MODEL[] = model
    return model
end

# ── Blow-Wood constants ───────────────────────────────────────────────────────
const RAMAN_TAU1_SI = 12.2e-15      # s
const RAMAN_TAU2_SI = 32.0e-15      # s
const RAMAN_FR_SI   = 0.18          # fractional Raman contribution
const RAMAN_OMEGA_R = 2π * 13.2e12  # rad/s, Raman peak

# ── Numerical helpers (private; reused by brillouin.jl) ───────────────────────
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

"Raman gain coefficient using the tabulated spectrum."
function g_R_tabulated(Ω, γ; fR=RAMAN_FR_SI)
    return 2 .* γ .* fR .* im_h_R_tabulated(Ω)
end

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
    g_R(Ω, γ; fR=RAMAN_FR_SI, model=nothing) -> W⁻¹m⁻¹

Raman gain coefficient. `Ω` is the angular frequency shift from pump (rad/s);
`γ` is the fiber nonlinear coefficient (W⁻¹m⁻¹). Stokes (`Ω > 0`) gives
positive gain; anti-Stokes (`Ω < 0`) gives negative.
"""
function g_R(Ω, γ; fR=RAMAN_FR_SI, model=nothing)
    m = something(model, RAMAN_MODEL[])
    if m === :bw
        im_hR = imag.(h_R_freq(Ω))
    elseif m === :tabulated
        im_hR = im_h_R_tabulated(Ω)
    elseif m === :hc
        im_hR = im_h_R_hc(Ω)
    else
        throw(ArgumentError("Unknown Raman model $(repr(m)). Choose from $VALID_RAMAN_MODELS."))
    end
    return 2 .* γ .* fR .* im_hR
end

"`g_R` between two wavelengths (m); shift is `2π c (1/λ_pump − 1/λ_signal)`."
function g_R_from_wavelengths(λ_pump, λ_signal, γ; fR=RAMAN_FR_SI, model=nothing)
    Ω = 2π * SPEED_OF_LIGHT_M_PER_S * (1/λ_pump - 1/λ_signal)
    return float(g_R(Ω, γ; fR=fR, model=model))
end

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
                             sideband=:stokes, pump_depletion=false, model=nothing)

Spectral density `dṄ/dΩ` of spontaneous Raman photons
[photons s⁻¹ (rad/s)⁻¹]. Convention: `Ω = ω_pump − ω_signal`; positive `Ω`
is Stokes.
"""
function sprs_photon_rate_density(Ω, ω_pump, P_pump, L, γ, T_K;
                                    sideband::Symbol=:stokes,
                                    pump_depletion::Bool=false,
                                    model=nothing)
    Ω_for_gR = sideband === :antistokes ? abs.(Ω) : Ω
    gR  = g_R(Ω_for_gR, γ; model=model)
    nth = thermal_phonon_number(Ω, T_K)
    factor = sideband === :stokes ? (nth .+ 1.0) : nth
    rate = _sprs_rate_density_core(gR, P_pump, L, factor; pump_depletion=pump_depletion)
    return max.(rate, 0.0)
end

"""
    sprs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, γ, T_K;
                          sideband=:stokes, pump_depletion=false,
                          n_points=1001, model=nothing)
    sprs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                          T_K=fiber.T_ref_K, length_m=nothing, kwargs...)

Total spontaneous Raman photon rate into a quantum channel [photons / s].

The first form takes raw `(L, γ, T_K)`. The second derives them from a
`Fiber`: `γ = nonlinear_coefficient(fiber.cross_section, λ_pump, T_K)` and
`L = fiber.s_end − fiber.s_start` (override with `length_m`).
"""
function sprs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, γ, T_K;
                                 sideband::Symbol=:stokes,
                                 pump_depletion::Bool=false,
                                 n_points::Int=1001,
                                 model=nothing)
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
                                      model=model)
    return float(_simpson(rate, ω_s_arr))
end

function sprs_noise_in_channel(fiber::Fiber, λ_pump, λ_channel, Δλ, P_pump;
                                 T_K=fiber.T_ref_K,
                                 length_m=nothing,
                                 kwargs...)
    γ = nonlinear_coefficient(fiber.cross_section, λ_pump, T_K)
    L = length_m === nothing ? fiber.s_end - fiber.s_start : length_m
    return sprs_noise_in_channel(λ_pump, λ_channel, Δλ, P_pump, L, γ, T_K; kwargs...)
end

"""
    check_raman_depletion_validity(P_pump, L, γ; fR=RAMAN_FR_SI, threshold=0.1)

Warn if the pump-depletion correction `[exp(x)−1]/x` (with `x = g_R·P·L`)
exceeds `threshold` (10 % default). Returns the ratio.
"""
function check_raman_depletion_validity(P_pump, L, γ; fR=RAMAN_FR_SI, threshold=0.1)
    gR_peak = g_R(RAMAN_OMEGA_R, γ; fR=fR, model=:bw)
    x = gR_peak * P_pump * L
    ratio = x > 1e-10 ? expm1(x) / x : 1.0
    if ratio - 1.0 > threshold
        @warn @sprintf("Raman pump-depletion correction is %.1f%% (x = gR·P·L = %.3f). Pass pump_depletion=true.",
                       (ratio - 1) * 100, x)
    end
    return ratio
end
