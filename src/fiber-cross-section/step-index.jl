#################################################
#
# Base Constants and Structure
#
#################################################

"""
    STEP_INDEX_GLASS

Union of the glass types accepted as step-index core or cladding material.
"""
const STEP_INDEX_GLASS = Union{SilicaGermaniaGlass, SilicaFluorinatedGlass}

"""
    LP11_CUTOFF_V

Normalized frequency `V = 2.405` of the LP11 cutoff; the fiber is single-mode
below this value.
"""
const LP11_CUTOFF_V = 2.405

"""
    MARCUSE_V_MIN

Lower bound of the `V` range for which the modified Marcuse mode-area
approximation in [`effective_mode_area`](@ref) is calibrated.
"""
const MARCUSE_V_MIN = 1.2

"""
    MARCUSE_V_MAX

Upper bound of the `V` range for which the modified Marcuse mode-area
approximation in [`effective_mode_area`](@ref) is calibrated.
"""
const MARCUSE_V_MAX = 2.4

"""
    StepIndexCrossSection(core_material, cladding_material, core_diameter_m,
                          cladding_diameter_m; manufacturer = nothing,
                          model_number = nothing, ellipticity_axis_ratio = 1.0,
                          ellipticity_axis_angle = 0.0)

Ideal step-index fiber cross section.

The baseline object is a circular step-index cross section described by
core/cladding materials ([`STEP_INDEX_GLASS`](@ref)) and diameters (m); the core
diameter must be smaller than the cladding diameter. Environmental perturbations
such as bending, axial tension, and twist are handled by response functions with
explicit arguments rather than stored on the type.

Intrinsic transverse core ellipticity is described by two fields:

- `ellipticity_axis_ratio` is the major/minor magnitude and must be `≥ 1`
  (`1` ⇒ circular ⇒ no ellipticity birefringence). Requiring `≥ 1` makes the
  major axis — and hence the angle below — unambiguous: an ellipse with ratio
  `r < 1` is the same object as ratio `1/r` rotated 90°, so callers express it
  canonically as ratio `≥ 1`. The cross-section returns only the birefringence
  *magnitude*; the fiber generator orients it.
- `ellipticity_axis_angle` (rad) is the **major-axis** orientation in the local
  transverse frame, canonical in `[0, π)` (an ellipse is 180°-symmetric). At
  angle `0` the major axis is aligned with the curvature normal (the same
  direction the bend birefringence picks out); increasing the angle rotates it
  toward the binormal, about the propagation tangent.

# Examples
```julia
fiber = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

λ = 1550e-9
T = 297.15

V = normalized_frequency(fiber, λ, T)
β = propagation_constant(fiber, λ, T)
Aeff = effective_mode_area(fiber, λ, T)
Δβ = bending_birefringence(fiber, λ, T; bend_radius_m = 0.03)
```
"""
struct StepIndexCrossSection{T<:Real} <: FiberCrossSection
    manufacturer::Union{Nothing, String}
    model_number::Union{Nothing, String}
    core_material::STEP_INDEX_GLASS
    cladding_material::STEP_INDEX_GLASS
    core_diameter_m::T
    cladding_diameter_m::T
    ellipticity_axis_ratio::T
    ellipticity_axis_angle::T

    function StepIndexCrossSection(
        core_material::STEP_INDEX_GLASS,
        cladding_material::STEP_INDEX_GLASS,
        core_diameter_m::Real,
        cladding_diameter_m::Real;
        manufacturer::Union{Nothing, AbstractString} = nothing,
        model_number::Union{Nothing, AbstractString} = nothing,
        ellipticity_axis_ratio::Real = 1.0,
        ellipticity_axis_angle::Real = 0.0
    )
        core_diameter = validate_positive_length(core_diameter_m, "core diameter")
        cladding_diameter = validate_positive_length(cladding_diameter_m, "cladding diameter")
        core_diameter < cladding_diameter || throw(ArgumentError(
            "core diameter must be smaller than cladding diameter"
        ))
        axis_ratio = validate_axis_ratio(ellipticity_axis_ratio)

        T = promote_type(typeof(core_diameter), typeof(cladding_diameter),
                         typeof(axis_ratio), typeof(ellipticity_axis_angle))
        new{T}(
            isnothing(manufacturer) ? nothing : String(manufacturer),
            isnothing(model_number) ? nothing : String(model_number),
            core_material,
            cladding_material,
            core_diameter,
            cladding_diameter,
            axis_ratio,
            ellipticity_axis_angle
        )
    end
end

#################################################
#
# Convenience Methods
#
#################################################

"""
    core_radius(fiber)

Return the core radius (m).
"""
core_radius(fiber::StepIndexCrossSection) = fiber.core_diameter_m / 2

"""
    cladding_radius(fiber)

Return the cladding radius (m).
"""
cladding_radius(fiber::StepIndexCrossSection) = fiber.cladding_diameter_m / 2

"""
    core_refractive_index([style], fiber, λ, T_K)

Return the refractive index of the core material.

Many spectral queries here take an optional leading `style::SpectralStyle`,
written `[style]` in their signatures. The omitted-argument form is
[`ValueOnly`](@ref) and returns the plain value; [`WithDerivative`](@ref)
returns a [`SpectralResponse`](@ref) carrying the value and its
angular-frequency derivative. Subsequent methods sharing this convention do not
repeat it.
"""
function core_refractive_index(style::SpectralStyle, fiber::StepIndexCrossSection, λ, T_K)
    return refractive_index(style, fiber.core_material, λ, T_K)
end

core_refractive_index(fiber::StepIndexCrossSection, λ, T_K) =
    core_refractive_index(ValueOnly(), fiber, λ, T_K)

"""
    cladding_refractive_index([style], fiber, λ, T_K)

Return the refractive index of the cladding material.
"""
function cladding_refractive_index(
    style::SpectralStyle,
    fiber::StepIndexCrossSection,
    λ,
    T_K
)
    return refractive_index(style, fiber.cladding_material, λ, T_K)
end

cladding_refractive_index(fiber::StepIndexCrossSection, λ, T_K) =
    cladding_refractive_index(ValueOnly(), fiber, λ, T_K)

#################################################
#
# Base Quantities (needed for subsequent calculations)
#
#################################################

"""
    waveguide_factor(V)

Return the Gloge approximation to the LP01 transverse core parameter `u(V)`,
`u ≈ (1 + √2)·V / (1 + (4 + V⁴)^(1/4))`.

Gloge, "Weakly guiding fibers", doi:10.1364/AO.10.002252.
"""
function waveguide_factor(V)
    α = one(V) + sqrt(2 * one(V))
    t = (4 + V^4)^(one(V) / 4)
    return α * V / (one(V) + t)
end

"""
    waveguide_factor_prime(V)

Return `d/dV` of [`waveguide_factor`](@ref).
"""
function waveguide_factor_prime(V)
    α = one(V) + sqrt(2 * one(V))
    t = (4 + V^4)^(one(V) / 4)
    dt_dV = V^3 / (4 + V^4)^(3 * one(V) / 4)
    den = one(V) + t
    return α * (den - V * dt_dV) / den^2
end

"""
    modal_prefactor(V)

Return `1 - u²/V²` using the Gloge approximation [`waveguide_factor`](@ref) for
`u(V)`; this prefactor appears in several birefringence formulas.
"""
function modal_prefactor(V)
    α = one(V) + sqrt(2*one(V))
    t = one(V) + (4 + V^4)^(one(V) / 4)
    return one(V) - α^2 / t^2
end

"""
    modal_prefactor_prime(V)

Return `d/dV` of [`modal_prefactor`](@ref).
"""
function modal_prefactor_prime(V)
    α = one(V) + sqrt(2*one(V))
    t = (4 + V^4)^(one(V) / 4)
    dt_dV = V^3 / (4 + V^4)^(3 * one(V) / 4)
    den = one(V) + t
    return 2 * α^2 * dt_dV / den^3
end

"""
    mode_terms(style, fiber, λ, T_K) -> NamedTuple

Return the guided-mode quantities shared by the downstream formulas: radii,
vacuum wavenumber `k0`, core/cladding indices, numerical aperture `na`,
normalized frequency `V`, Gloge `waveguide_factor` and `modal_prefactor`, and
propagation constant `β`, each paired with its angular-frequency derivative
(`dk0_dω`, `dn_core_dω`, …).

With `ValueOnly()` the derivative slots are zero; with `WithDerivative()` they
carry the chain-ruled ω-derivatives. Throws an `ArgumentError` unless
`n_core > n_cladding`.
"""
function mode_terms(::ValueOnly, fiber::StepIndexCrossSection, λ, T_K)
    n_core = core_refractive_index(fiber, λ, T_K)
    n_clad = cladding_refractive_index(fiber, λ, T_K)
    n_core > n_clad || throw(ArgumentError(
        "guided-mode calculations require n_core > n_cladding; got " *
        "n_core=$(n_core), n_cladding=$(n_clad)"
    ))

    r_core = core_radius(fiber)
    r_clad = cladding_radius(fiber)
    k0 = 2π / λ
    dk0_dω = one(λ) / SPEED_OF_LIGHT_M_PER_S
    na = sqrt(n_core^2 - n_clad^2)
    V = r_core * k0 * na
    g = waveguide_factor(V)
    q = modal_prefactor(V)
    β = sqrt((n_core^2 * k0^2) - g^2 / r_core^2)
    z = zero(β)

    return (
        core_radius = r_core,
        cladding_radius = r_clad,
        k0 = k0,
        dk0_dω = dk0_dω,
        n_core = n_core,
        dn_core_dω = z,
        n_clad = n_clad,
        dn_clad_dω = z,
        na = na,
        dna_dω = z,
        V = V,
        dV_dω = z,
        waveguide_factor = g,
        dwaveguide_factor_dω = z,
        modal_prefactor = q,
        dmodal_prefactor_dω = z,
        β = β,
        dβ_dω = z
    )
end

function mode_terms(::WithDerivative, fiber::StepIndexCrossSection, λ, T_K)
    n_core_resp = core_refractive_index(WithDerivative(), fiber, λ, T_K)
    n_clad_resp = cladding_refractive_index(WithDerivative(), fiber, λ, T_K)
    n_core = n_core_resp.value
    n_clad = n_clad_resp.value
    n_core > n_clad || throw(ArgumentError(
        "guided-mode calculations require n_core > n_cladding; got " *
        "n_core=$(n_core), n_cladding=$(n_clad)"
    ))

    r_core = core_radius(fiber)
    r_clad = cladding_radius(fiber)
    k0 = 2π / λ
    dk0_dω = one(λ) / SPEED_OF_LIGHT_M_PER_S
    na = sqrt(n_core^2 - n_clad^2)
    dna_dω = (n_core * n_core_resp.dω - n_clad * n_clad_resp.dω) / na
    V = r_core * k0 * na
    dV_dω = r_core * (dk0_dω * na + k0 * dna_dω)
    g = waveguide_factor(V)
    dg_dω = waveguide_factor_prime(V) * dV_dω
    q = modal_prefactor(V)
    dq_dω = modal_prefactor_prime(V) * dV_dω
    β = sqrt((n_core^2) * k0^2 - g^2 / r_core^2)
    dβ_inner_dω = 2 * n_core * n_core_resp.dω * k0^2 +
                  2 * n_core^2 * k0 * dk0_dω -
                  2 * g * dg_dω / r_core^2
    dβ_dω = dβ_inner_dω / (2 * β)

    return (
        core_radius = r_core,
        cladding_radius = r_clad,
        k0 = k0,
        dk0_dω = dk0_dω,
        n_core = n_core,
        dn_core_dω = n_core_resp.dω,
        n_clad = n_clad,
        dn_clad_dω = n_clad_resp.dω,
        na = na,
        dna_dω = dna_dω,
        V = V,
        dV_dω = dV_dω,
        waveguide_factor = g,
        dwaveguide_factor_dω = dg_dω,
        modal_prefactor = q,
        dmodal_prefactor_dω = dq_dω,
        β = β,
        dβ_dω = dβ_dω
    )
end

"""
    relative_index_difference([style], fiber, λ, T_K)

Return the relative index difference `Δ = (n_core - n_clad) / n_clad`
(dimensionless).
"""
function relative_index_difference(style::ValueOnly, fiber::StepIndexCrossSection, λ, T_K)
    n_core = core_refractive_index(fiber, λ, T_K)
    n_clad = cladding_refractive_index(fiber, λ, T_K)
    return (n_core - n_clad) / n_clad
end

relative_index_difference(fiber::StepIndexCrossSection, λ, T_K) =
    relative_index_difference(ValueOnly(), fiber, λ, T_K)

"""
    numerical_aperture([style], fiber, λ, T_K)

Return the numerical aperture `√(n_core² - n_clad²)`.
"""
numerical_aperture(style::ValueOnly, fiber::StepIndexCrossSection, λ, T_K) =
    mode_terms(style, fiber, λ, T_K).na

function numerical_aperture(style::WithDerivative, fiber::StepIndexCrossSection, λ, T_K)
    terms = mode_terms(style, fiber, λ, T_K)
    return SpectralResponse(terms.na, terms.dna_dω)
end

numerical_aperture(fiber::StepIndexCrossSection, λ, T_K) =
    numerical_aperture(ValueOnly(), fiber, λ, T_K)

"""
    normalized_frequency([style], fiber, λ, T_K)

Return the normalized frequency `V = r_core · k0 · NA`.
"""
normalized_frequency(style::ValueOnly, fiber::StepIndexCrossSection, λ, T_K) =
    mode_terms(style, fiber, λ, T_K).V

function normalized_frequency(style::WithDerivative, fiber::StepIndexCrossSection, λ, T_K)
    terms = mode_terms(style, fiber, λ, T_K)
    return SpectralResponse(terms.V, terms.dV_dω)
end

normalized_frequency(fiber::StepIndexCrossSection, λ, T_K) =
    normalized_frequency(ValueOnly(), fiber, λ, T_K)

"""
    propagation_constant([style], fiber, λ, T_K)

Return the LP01 propagation constant `β` (rad/m).
"""
propagation_constant(style::ValueOnly, fiber::StepIndexCrossSection, λ, T_K) =
    mode_terms(style, fiber, λ, T_K).β

function propagation_constant(style::WithDerivative, fiber::StepIndexCrossSection, λ, T_K)
    terms = mode_terms(style, fiber, λ, T_K)
    return SpectralResponse(terms.β, terms.dβ_dω)
end

propagation_constant(fiber::StepIndexCrossSection, λ, T_K) =
    propagation_constant(ValueOnly(), fiber, λ, T_K)

"""
    effective_mode_index([style], fiber, λ, T_K)

Return the effective mode index `n_eff = β / k0`.
"""
function effective_mode_index(style::ValueOnly, fiber::StepIndexCrossSection, λ, T_K)
    terms = mode_terms(style, fiber, λ, T_K)
    return terms.β / terms.k0
end

function effective_mode_index(style::WithDerivative, fiber::StepIndexCrossSection, λ, T_K)
    terms = mode_terms(style, fiber, λ, T_K)
    value = terms.β / terms.k0
    dω = (terms.dβ_dω * terms.k0 - terms.β * terms.dk0_dω) / terms.k0^2
    return SpectralResponse(value, dω)
end

effective_mode_index(fiber::StepIndexCrossSection, λ, T_K) =
    effective_mode_index(ValueOnly(), fiber, λ, T_K)

"""
    effective_group_index(fiber, λ, T_K)

Return the effective group index `n_g = n_eff + ω · dn_eff/dω`.
"""
function effective_group_index(
    fiber::StepIndexCrossSection,
    λ,
    T_K
)

    n_eff = effective_mode_index(WithDerivative(), fiber, λ, T_K)
    ω = 2*pi*SPEED_OF_LIGHT_M_PER_S/λ
    return n_eff.value + ω*n_eff.dω
end

"""
    effective_mode_area(fiber, λ, T_K)

Return the effective mode area (m²) by the 1/e² criterion, `π·w²` with the mode
field radius `w` from a modified Marcuse approximation that is better at lower
`V` (accurate to within 1% for `1.5 < V < 2.5`).

Warns when `V` is outside the calibrated range `[MARCUSE_V_MIN, MARCUSE_V_MAX]`
or above the single-mode cutoff `LP11_CUTOFF_V` (the form applies only to the
fundamental mode).
"""
function effective_mode_area(fiber::StepIndexCrossSection, λ, T_K)
    V = normalized_frequency(fiber, λ, T_K)
    if !(MARCUSE_V_MIN <= V <= MARCUSE_V_MAX)
        @warn "Marcuse effective-area approximation is calibrated for " *
              "1.2 <= V <= 2.4; got V=$(V)"
    end
    if (V > LP11_CUTOFF_V)
        @warn "Marcuse approximation only applies to the fundamental mode. " *
              "V=$(V) is above the single-mode cutoff."
    end
    w_over_r = 0.65 + 1.619 / V^1.5 + 2.879 / V^6 - 0.016 - 1.561 / V^7
    w = w_over_r * core_radius(fiber)
    return π * w^2
end

"""
    core_nonlinear_refractive_index(fiber, λ, T_K)

Return the nonlinear refractive index `n2` of the fiber core material in m²/W.
"""
core_nonlinear_refractive_index(fiber::StepIndexCrossSection, λ, T_K) =
    nonlinear_refractive_index(fiber.core_material, λ, T_K)

"""
    chromatic_dispersion_parameter(fiber, λ, T_K; dλ = 0.1e-9)

Return the chromatic dispersion parameter `D_CD = -2πc/λ² · d²k/dω²` in
ps/(nm·km), the wavelength-derivative form commonly quoted for fibers.

Uses a central finite difference of the effective mode index with wavelength
step `dλ` (m).
"""
function chromatic_dispersion_parameter(
    fiber::StepIndexCrossSection,
    λ,
    T_K;
    dλ = 0.1e-9
)
    n_center = effective_mode_index(fiber, λ, T_K)
    n_minus = effective_mode_index(fiber, λ - dλ, T_K)
    n_plus = effective_mode_index(fiber, λ + dλ, T_K)
    return -λ / SPEED_OF_LIGHT_M_PER_S * (n_plus - 2 * n_center + n_minus) / dλ^2 * 1e6
end

"""
    group_velocity_dispersion_parameter(fiber, λ, T_K; dλ = 0.1e-9)

Return the group-velocity dispersion `β₂ = d²k/dω²` in ps²/km, derived from
[`chromatic_dispersion_parameter`](@ref) with the same finite-difference step
`dλ` (m).
"""
function group_velocity_dispersion_parameter(
    fiber::StepIndexCrossSection,
    λ,
    T_K;
    dλ = 0.1e-9
)
    D_SI = chromatic_dispersion_parameter(fiber, λ, T_K; dλ = dλ) * 1e-6
    return -(λ^2 / (2π * SPEED_OF_LIGHT_M_PER_S)) * D_SI * 1e27
end

"""
    is_single_mode(fiber, λ, T_K) -> Bool

Return whether the fiber is single-mode at `λ`, i.e. `V < LP11_CUTOFF_V`.
"""
is_single_mode(fiber::StepIndexCrossSection, λ, T_K) =
    normalized_frequency(fiber, λ, T_K) < LP11_CUTOFF_V

"""
    cutoff_wavelength(fiber, T_K; λ_min, λ_max, atol, maxiter)

Bisect for the LP11 cutoff wavelength.

!!! note "MCM compatibility"
    `T_K` must be a scalar (not `Particles`). The bisection's `signbit` branching is
    not defined per-particle; call with `pmean(T_K_uncertain)` when the temperature
    is uncertain. The return is always `Float64`.
"""
function cutoff_wavelength(
    fiber::StepIndexCrossSection,
    T_K;
    λ_min::Real = MIN_VALID_WAVELENGTH_M,
    λ_max::Real = MAX_VALID_WAVELENGTH_M,
    atol::Real = 1e-12,
    maxiter::Integer = 200
)
    λ_lo = validate_positive_length(λ_min, "λ_min")
    λ_hi = validate_positive_length(λ_max, "λ_max")
    λ_lo < λ_hi || throw(ArgumentError("λ_min must be smaller than λ_max"))

    V_lo = normalized_frequency(fiber, λ_lo, T_K) - LP11_CUTOFF_V
    V_hi = normalized_frequency(fiber, λ_hi, T_K) - LP11_CUTOFF_V

    if V_lo == 0.0
        return λ_lo
    elseif V_hi == 0.0
        return λ_hi
    elseif signbit(V_lo) == signbit(V_hi)
        throw(ArgumentError(
            "cutoff wavelength is not bracketed in [$(λ_lo), $(λ_hi)] m"
        ))
    end

    a = λ_lo
    b = λ_hi
    Va = V_lo

    for _ in 1:maxiter
        mid = (a + b) / 2
        Vm = normalized_frequency(fiber, mid, T_K) - LP11_CUTOFF_V
        if abs(Vm) <= atol || (b - a) / 2 <= atol
            return mid
        elseif signbit(Vm) == signbit(Va)
            a = mid
            Va = Vm
        else
            b = mid
        end
    end

    return (a + b) / 2
end

"""
    eccentricity_squared(axis_ratio)

Return the squared eccentricity `e² = 1 - 1/ε²` of an ellipse with canonical
major/minor axis ratio `ε ≥ 1` (validated).
"""
function eccentricity_squared(axis_ratio)
    ε = validate_axis_ratio(axis_ratio)
    return one(ε) - inv(ε)^2
end

#################################################
#
# Birefringences
#
#################################################

# The cross-section returns birefringence **magnitudes** only; the fiber
# generator orients the eigen-axes (from `ellipticity_axis_angle` plus
# spin/twist phase for the intrinsic terms, or the curvature normal for the
# bend/tension terms).

"""
    core_noncircularity_dω(style, fiber, λ, T_K;
                           axis_ratio = fiber.ellipticity_axis_ratio)

Return the [`BirefringenceResponse`](@ref) of geometric core-noncircularity
linear birefringence for an elliptical core with major/minor ratio `axis_ratio`.

`axis_ratio` is canonical (`≥ 1`), so the magnitude is a single nonnegative
branch in `eccentricity_squared(ε) = 1 - 1/ε²`; orientation lives in the
major-axis angle applied by the fiber generator, not in the sign of the
magnitude. A circular core (`axis_ratio == 1`) returns zero.
"""
function core_noncircularity_dω(style::SpectralStyle, fiber::StepIndexCrossSection, λ, T_K;
                                axis_ratio = fiber.ellipticity_axis_ratio)
    terms = mode_terms(style, fiber, λ, T_K)
    ε = validate_axis_ratio(axis_ratio)
    ε == one(ε) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    χ = one(terms.n_core) - terms.n_clad^2 / terms.n_core^2
    dχ_dω = -2 * terms.n_clad * terms.dn_clad_dω / terms.n_core^2 +
            2 * terms.n_clad^2 * terms.dn_core_dω / terms.n_core^3
    V = terms.V
    h = 4 * log(V)^3 / (V^3 * (one(V) + log(V)))
    h_prime = h / V * (3 / log(V) - 3 - inv(one(V) + log(V)))
    prefactor = eccentricity_squared(ε) / terms.core_radius
    Δβ = prefactor * (2*χ)^(3 / 2) * h
    dω = prefactor * ((3 / 2) * sqrt(2*χ) * dχ_dω * h + (2*χ)^(3 / 2) * h_prime * terms.dV_dω)
    return BirefringenceResponse(Δβ, dω)
end

"""
    asymmetric_thermal_stress_dω(style, fiber, λ, T_K;
                                 axis_ratio = fiber.ellipticity_axis_ratio)

Return the [`BirefringenceResponse`](@ref) of the linear birefringence from
asymmetric thermal stress frozen into an elliptical core on cooling from the
softening temperature (response ∝ `|T_soft - T_K|`).

The ellipse-asymmetry factor `(ε - 1)/(ε + 1)` is nonnegative because
`axis_ratio` is canonical (`≥ 1`); the fiber generator supplies the orientation
from the major-axis angle, shared with [`core_noncircularity_dω`](@ref). A
circular core returns zero.
"""
function asymmetric_thermal_stress_dω(
    style::SpectralStyle,
    fiber::StepIndexCrossSection,
    λ,
    T_K;
    axis_ratio = fiber.ellipticity_axis_ratio
)
    terms = mode_terms(style, fiber, λ, T_K)
    ε = validate_axis_ratio(axis_ratio)
    ε == one(ε) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    α_core = cte(fiber.core_material, T_K)
    α_clad = cte(fiber.cladding_material, T_K)
    T_soft = softening_temperature(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    const_factor = 0.5 * (p11 - p12) * (α_clad - α_core) * abs(T_soft - T_K) /
                   (1 - ν^2) * ((ε - 1) / (ε + 1))
    Δβ = terms.k0 * terms.modal_prefactor * terms.n_core^3 * const_factor
    dω = const_factor * (
        terms.dk0_dω * terms.modal_prefactor * terms.n_core^3 +
        terms.k0 * terms.dmodal_prefactor_dω * terms.n_core^3 +
        terms.k0 * terms.modal_prefactor * 3 * terms.n_core^2 * terms.dn_core_dω
    )
    return BirefringenceResponse(Δβ, dω)
end

"""
    bending_dω(style, fiber, λ, T_K; bend_radius_m)

Return the [`BirefringenceResponse`](@ref) of photoelastic bend-induced linear
birefringence at bend radius `bend_radius_m` (m); the response scales as
`(r_clad / R)²` and is zero for `R = Inf` (straight).
"""
function bending_dω(
    style::SpectralStyle,
    fiber::StepIndexCrossSection,
    λ,
    T_K;
    bend_radius_m
)
    terms = mode_terms(style, fiber, λ, T_K)
    R = validate_bend_radius(bend_radius_m)
    isinf(R) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    geom = 0.5 * (terms.cladding_radius^2 / R^2)
    const_factor = (p11 - p12) * (1 + ν) * geom / 2
    Δβ = terms.k0 * terms.n_core^3 * const_factor
    dω = const_factor * (terms.dk0_dω * terms.n_core^3 + terms.k0 * 3 * terms.n_core^2 * terms.dn_core_dω)
    return BirefringenceResponse(Δβ, dω)
end

"""
    axial_tension_dω(style, fiber, λ, T_K; bend_radius_m, axial_tension_N)

Return the [`BirefringenceResponse`](@ref) of photoelastic linear birefringence
from axial tension `axial_tension_N` (N, nonnegative) on a fiber bent at radius
`bend_radius_m` (m).

The response scales as `(r_clad / R) · F / (π·r_clad²·E)` and vanishes for a
straight fiber (`R = Inf`) or zero tension. Throws an `ArgumentError` for
negative tension.
"""
function axial_tension_dω(
    style::SpectralStyle,
    fiber::StepIndexCrossSection,
    λ,
    T_K;
    bend_radius_m,
    axial_tension_N
)
    terms = mode_terms(style, fiber, λ, T_K)
    R = validate_bend_radius(bend_radius_m)
    axial_tension_N >= zero(axial_tension_N) ||
        throw(ArgumentError("axial_tension_N must be nonnegative"))
    (isinf(R) || axial_tension_N == zero(axial_tension_N)) &&
        return BirefringenceResponse(zero(terms.β), zero(terms.β))

    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    E = youngs_modulus(fiber.core_material, T_K)
    geom = ((2 - 3 * ν) / (1 - ν)) * (terms.cladding_radius / R) *
           (axial_tension_N / (π * terms.cladding_radius^2 * E))
    const_factor = (p11 - p12) * (1 + ν) * geom / 2
    Δβ = terms.k0 * terms.n_core^3 * const_factor
    dω = const_factor * (terms.dk0_dω * terms.n_core^3 + terms.k0 * 3 * terms.n_core^2 * terms.dn_core_dω)
    return BirefringenceResponse(Δβ, dω)
end

"""
    twisting_dω(style, fiber, λ, T_K; twist_rate_rad_per_m)

Return the [`BirefringenceResponse`](@ref) of twist-induced **circular**
birefringence (optical activity) from the photoelastic effect: the difference in
propagation constants of the two circular polarizations, `βLC − βRC`, per unit
mechanical-twist rate.

The polarization tracks the geometric rotation of the medium (the leading `1`)
reduced by the photoelastic slip `n²(p₁₁−p₁₂)/2` (negative for silica), so the
net rotation rate is `(1 + n²(p₁₁−p₁₂)/2)·τ_m`. The fiber generator places this
on the real antisymmetric (rotation) part of `K`; see
`circular_birefringence_generator`.

This convention agrees with the experimentally validated treatment in
doi:10.1016/j.yofte.2011.10.001 and differs from the seminal Rashleigh paper,
doi:10.1109/JLT.1983.1072121.

Throws an `ArgumentError` for a non-finite twist rate; a zero rate returns zero.
"""
function twisting_dω(
    style::SpectralStyle,
    fiber::StepIndexCrossSection,
    λ,
    T_K;
    twist_rate_rad_per_m
)
    terms = mode_terms(style, fiber, λ, T_K)
    tr = float(twist_rate_rad_per_m)
    isfinite(tr) || throw(ArgumentError("twist_rate_rad_per_m must be finite"))
    tr == zero(tr) && return BirefringenceResponse(zero(terms.β), zero(terms.β))
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    coeff = (p11 - p12) / 2
    Δβ = (one(terms.n_core) + coeff * terms.n_core^2) * tr
    # d/dω of the leading `1·tr` term is zero; only the n² term carries dispersion.
    dω = 2 * coeff * terms.n_core * terms.dn_core_dω * tr
    return BirefringenceResponse(Δβ, dω)
end

"""
    core_noncircularity_birefringence([style], fiber, λ, T_K;
                                      axis_ratio = fiber.ellipticity_axis_ratio)

Return the core-noncircularity linear birefringence magnitude (rad/m); with
`WithDerivative()`, return the full [`BirefringenceResponse`](@ref). See
[`core_noncircularity_dω`](@ref).
"""
core_noncircularity_birefringence(::ValueOnly, fiber::StepIndexCrossSection, λ, T_K;
                                  axis_ratio = fiber.ellipticity_axis_ratio) =
    core_noncircularity_dω(ValueOnly(), fiber, λ, T_K; axis_ratio = axis_ratio).Δβ

core_noncircularity_birefringence(::WithDerivative, fiber::StepIndexCrossSection, λ, T_K;
                                  axis_ratio = fiber.ellipticity_axis_ratio) =
    core_noncircularity_dω(WithDerivative(), fiber, λ, T_K; axis_ratio = axis_ratio)

core_noncircularity_birefringence(fiber::StepIndexCrossSection, λ, T_K;
                                  axis_ratio = fiber.ellipticity_axis_ratio) =
    core_noncircularity_birefringence(ValueOnly(), fiber, λ, T_K; axis_ratio = axis_ratio)

"""
    asymmetric_thermal_stress_birefringence([style], fiber, λ, T_K;
                                            axis_ratio = fiber.ellipticity_axis_ratio)

Return the asymmetric-thermal-stress linear birefringence magnitude (rad/m);
with `WithDerivative()`, return the full [`BirefringenceResponse`](@ref). See
[`asymmetric_thermal_stress_dω`](@ref).
"""
asymmetric_thermal_stress_birefringence(::ValueOnly, fiber::StepIndexCrossSection, λ, T_K;
                                        axis_ratio = fiber.ellipticity_axis_ratio) =
    asymmetric_thermal_stress_dω(ValueOnly(), fiber, λ, T_K; axis_ratio = axis_ratio).Δβ

asymmetric_thermal_stress_birefringence(::WithDerivative, fiber::StepIndexCrossSection,
                                        λ, T_K;
                                        axis_ratio = fiber.ellipticity_axis_ratio) =
    asymmetric_thermal_stress_dω(WithDerivative(), fiber, λ, T_K; axis_ratio = axis_ratio)

asymmetric_thermal_stress_birefringence(fiber::StepIndexCrossSection, λ, T_K;
                                        axis_ratio = fiber.ellipticity_axis_ratio) =
    asymmetric_thermal_stress_birefringence(ValueOnly(), fiber, λ, T_K;
                                            axis_ratio = axis_ratio)

"""
    bending_birefringence([style], fiber, λ, T_K; bend_radius_m)

Return the bend-induced linear birefringence magnitude (rad/m); with
`WithDerivative()`, return the full [`BirefringenceResponse`](@ref). See
[`bending_dω`](@ref).
"""
bending_birefringence(::ValueOnly, fiber::StepIndexCrossSection, λ, T_K; bend_radius_m) =
    bending_dω(ValueOnly(), fiber, λ, T_K; bend_radius_m = bend_radius_m).Δβ

bending_birefringence(::WithDerivative, fiber::StepIndexCrossSection, λ, T_K;
                      bend_radius_m) =
    bending_dω(WithDerivative(), fiber, λ, T_K; bend_radius_m = bend_radius_m)

bending_birefringence(fiber::StepIndexCrossSection, λ, T_K; bend_radius_m) =
    bending_birefringence(ValueOnly(), fiber, λ, T_K; bend_radius_m = bend_radius_m)

"""
    axial_tension_birefringence([style], fiber, λ, T_K; bend_radius_m,
                                axial_tension_N)

Return the axial-tension linear birefringence magnitude (rad/m); with
`WithDerivative()`, return the full [`BirefringenceResponse`](@ref). See
[`axial_tension_dω`](@ref).
"""
axial_tension_birefringence(::ValueOnly, fiber::StepIndexCrossSection, λ, T_K;
                            bend_radius_m, axial_tension_N) =
    axial_tension_dω(ValueOnly(), fiber, λ, T_K;
                     bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N).Δβ

axial_tension_birefringence(::WithDerivative, fiber::StepIndexCrossSection, λ, T_K;
                            bend_radius_m, axial_tension_N) =
    axial_tension_dω(WithDerivative(), fiber, λ, T_K;
                     bend_radius_m = bend_radius_m, axial_tension_N = axial_tension_N)

axial_tension_birefringence(fiber::StepIndexCrossSection, λ, T_K;
                            bend_radius_m, axial_tension_N) =
    axial_tension_birefringence(ValueOnly(), fiber, λ, T_K;
                                bend_radius_m = bend_radius_m,
                                axial_tension_N = axial_tension_N)

"""
    twisting_birefringence([style], fiber, λ, T_K; twist_rate_rad_per_m)

Return the twist-induced circular birefringence magnitude (rad/m); with
`WithDerivative()`, return the full [`BirefringenceResponse`](@ref). See
[`twisting_dω`](@ref).
"""
twisting_birefringence(::ValueOnly, fiber::StepIndexCrossSection, λ, T_K;
                       twist_rate_rad_per_m) =
    twisting_dω(ValueOnly(), fiber, λ, T_K;
                twist_rate_rad_per_m = twist_rate_rad_per_m).Δβ

twisting_birefringence(::WithDerivative, fiber::StepIndexCrossSection, λ, T_K;
                       twist_rate_rad_per_m) =
    twisting_dω(WithDerivative(), fiber, λ, T_K;
                twist_rate_rad_per_m = twist_rate_rad_per_m)

twisting_birefringence(fiber::StepIndexCrossSection, λ, T_K; twist_rate_rad_per_m) =
    twisting_birefringence(ValueOnly(), fiber, λ, T_K;
                           twist_rate_rad_per_m = twist_rate_rad_per_m)
