"""
Material properties for optical glasses.

This file defines base types and common algebraic methods that will
be available to all specific materials defined in the material folder.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- x_ge, x_f             dopant molar fraction (dimensionless, 0..1)
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa
- nonlinear_refractive_index (n_2)  m²/W


Specific materials should define all of the following.

refractive_index(::ValueOnly, material::AbstractMaterial, λ, T_K)
refractive_index(::WithDerivative, material::AbstractMaterial, λ, T_K)
cte(material::AbstractMaterial, T_K)
softening_temperature(material::AbstractMaterial, T_K)
poisson_ratio(material::AbstractMaterial, T_K)
photoelastic_constants(material::AbstractMaterial, T_K)
youngs_modulus(material::AbstractMaterial, T_K)
nonlinear_refractive_index(material::AbstractMaterial, λ, T_K)

"""

const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0

abstract type AbstractMaterial end
abstract type SpectralStyle end

# dω is the ω-derivative at the parameters that produce value.
struct SpectralResponse{T}
    value::T
    dω::T
end

struct ValueOnly <: SpectralStyle end
struct WithDerivative <: SpectralStyle end

# Shared by all materials to simplify the ValueOnly() case;
# Users do not need to copy or override.
refractive_index(material::AbstractMaterial, λ, T_K) =
    refractive_index(ValueOnly(), material, λ, T_K)

#################################################
#
# Validity ranges
#
#################################################

"""
    ValidityRange(lo, hi, name)

Closed interval `[lo, hi]` bounding the domain over which a material model is
valid, for a quantity labelled `name` (used only in error messages).
"""
struct ValidityRange
    lo::Float64
    hi::Float64
    name::String
end

"""
    check_range(value, range::ValidityRange) -> value
    check_range(values::NamedTuple, ranges::NamedTuple) -> nothing

Validate `value` against `range`, returning it (handy inside constructors) or
throwing an `ArgumentError` when it is non-finite or outside `[range.lo,
range.hi]`. The `NamedTuple` method validates several values at once, matching
each to the range that shares its key.

The check is MCM-safe: it uses only `isfinite` and `lo <= value <= hi`, with no
coercion or branching on uncertain values, so `Particles` pass through unchanged.
"""
function check_range(value, r::ValidityRange)
    isfinite(value) && r.lo <= value <= r.hi ||
        throw(ArgumentError("$(r.name) = $(value) outside [$(r.lo), $(r.hi)]"))
    return value
end

check_range(values::NamedTuple, ranges::NamedTuple) =
    foreach(k -> check_range(getfield(values, k), getfield(ranges, k)), keys(ranges))

"""
    runtime_ranges(material) -> NamedTuple

Per-material runtime validity envelope as a `NamedTuple` of [`ValidityRange`](@ref)s
keyed by quantity (for example `(; T_K, λ)`). Defaults to empty; each material
declares its own. Consuming layers read these bounds directly — the step-index
cutoff search, for instance, takes its bisection limits from
`runtime_ranges(material).λ`.
"""
runtime_ranges(::AbstractMaterial) = NamedTuple()

#################################################
#
# Common interpolation utilities
#
#################################################

# Used for simple scalar interpolation
interpolate_scalar(a, b, x) = (one(x) - x) * a + x * b

# Used for photoelastic_constants which are paired
function interpolate_pair(a::Tuple, b::Tuple, x)
    return (
        interpolate_scalar(a[1], b[1], x),
        interpolate_scalar(a[2], b[2], x)
    )
end

#################################################
#
# Sellmeier coefficient calculation methods
#
#################################################

"""
Most materials are modeled with the Sellmeier equation:
``n^2 = 1 + \\sum_{i=1}^n B_i\\lambda^2/(\\lambda^2 - C_i^2)`` where ``B_i`` and ``C_i`` are
strength and wavelength properties of each resonance used in the
calculation. The functions

    sellmeier_index_from_coefficients(coeffs, λ) -> Float64
    sellmeier_index_from_coefficients_dω(coeffs, λ) -> SpectralResponse

provide the refractive index given the Sellmeier coefficients B and C and the wavelength λ, with
coeffs specified as an n-tuple of 2-tuples (B, C). Thus, if appropriate, your material could
simply implement

    refractive_index(::ValueOnly, material::YourMaterial, λ, T_K) = 
        sellmeier_index_from_coefficients(YOUR_COEFFICIENTS, λ)
    refractive_index(::WithDerivative, material::YourMaterial, λ, T_K) = 
        sellmeier_index_from_coefficients_dω(YOUR_COEFFICIENTS, λ)

Sellmeier coefficients are often functions of other parameters such as temperature or molar 
fraction of a dopant. We provide the _evaluate_sellmeier_polynomials(B_coeffs, C_coeffs, x)
and _evaluate_sellmeier_constants(coeffs, x) utilities; see silica.jl and germania.jl for
examples of their use. 

These evaluators are pure numerics and do no domain validation: callers validate their
inputs at the `refractive_index` entry via `check_range` against the material's
`runtime_ranges`.

Note also that any implemented material must ensure compatibility with `Particles` to allow
Monte Carlo calculation. This happens naturally through the Sellmeier structure included here.
"""

function _evaluate_sellmeier_polynomials(B_coeffs, C_coeffs, x)
    return map((B, C) -> (evalpoly(x, B), evalpoly(x, C)), B_coeffs, C_coeffs)
end

function _evaluate_sellmeier_constants(coeffs, x)
    scale = one(x)
    return map(coeffs) do (B, C)
        (B * scale, C * scale)
    end
end

function sellmeier_index_from_coefficients(coeffs, λ)
    λ_um = λ * 1e6
    total = one(λ_um)
    for (B, C) in coeffs
        total += B * λ_um^2 / (λ_um^2 - C^2)
    end
    return sqrt(total)
end

function sellmeier_index_from_coefficients_dω(coeffs, λ)
    λ_um = λ * 1e6
    total = one(λ_um)
    dtotal_dλm = zero(λ_um)
    for (B, C) in coeffs
        denom = λ_um^2 - C^2
        total += B * λ_um^2 / denom
        dterm_dλum = -2 * B * λ_um * C^2 / denom^2
        dtotal_dλm += dterm_dλum * 1e6
    end
    n = sqrt(total)
    dn_dλ = dtotal_dλm / (2 * n)
    dλ_dω = -(λ^2) / (2π * SPEED_OF_LIGHT_M_PER_S)
    return SpectralResponse(n, dn_dλ * dλ_dω)
end
