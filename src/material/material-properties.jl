# Base types and shared spectral/interpolation utilities for the material layer.
# Concrete materials live in the sibling files of `src/material/`.

"""
    SPEED_OF_LIGHT_M_PER_S

Speed of light in vacuum (m/s).
"""
const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0

"""
    AbstractMaterial

Supertype for optical glass materials.

# Implementation

A concrete material must implement the following methods, with wavelength `λ` in
metres and temperature `T_K` in kelvin:

- `refractive_index(::ValueOnly, material, λ, T_K)` — refractive index
  (dimensionless).
- `refractive_index(::WithDerivative, material, λ, T_K)` — [`SpectralResponse`](@ref)
  carrying the index and its angular-frequency derivative.
- `cte(material, T_K)` — linear coefficient of thermal expansion (1/K).
- `softening_temperature(material, T_K)` — softening temperature (K).
- `poisson_ratio(material, T_K)` — Poisson ratio (dimensionless).
- `photoelastic_constants(material, T_K)` — `(p11, p12)` (dimensionless).
- `youngs_modulus(material, T_K)` — Young's modulus (Pa).
- `nonlinear_refractive_index(material, λ, T_K)` — nonlinear index `n₂` (m²/W).

Property methods must lift through `MonteCarloMeasurements.Particles` on their
uncertain-input slots (this happens naturally through the Sellmeier utilities in
this file).
"""
abstract type AbstractMaterial end

"""
    SpectralStyle

Dispatch flag selecting whether a spectral quantity is returned as a plain value
([`ValueOnly`](@ref)) or together with its angular-frequency derivative
([`WithDerivative`](@ref)).
"""
abstract type SpectralStyle end

"""
    SpectralResponse(value, dω)

Pair a spectral quantity with its angular-frequency derivative.

`dω` is the ω-derivative evaluated at the same parameters that produce `value`.
"""
struct SpectralResponse{T}
    value::T
    dω::T
end

"""
    ValueOnly()

`SpectralStyle` selecting the plain value of a spectral quantity.
"""
struct ValueOnly <: SpectralStyle end

"""
    WithDerivative()

`SpectralStyle` selecting a [`SpectralResponse`](@ref) that carries the value and
its angular-frequency derivative.
"""
struct WithDerivative <: SpectralStyle end

"""
    refractive_index(material, λ, T_K)

Return the refractive index of `material` at wavelength `λ` (m) and temperature
`T_K` (K).

Equivalent to `refractive_index(ValueOnly(), material, λ, T_K)`. Provided for all
materials, so implementations define only the style-explicit methods.
"""
refractive_index(material::AbstractMaterial, λ, T_K) =
    refractive_index(ValueOnly(), material, λ, T_K)

#################################################
#
# Validation utilities
#
#################################################

"""
    MIN_VALID_TEMPERATURE_K

Lower bound (K) of the temperature range the material models are valid for.
"""
const MIN_VALID_TEMPERATURE_K = 243.0

"""
    MAX_VALID_TEMPERATURE_K

Upper bound (K) of the temperature range the material models are valid for.
"""
const MAX_VALID_TEMPERATURE_K = 373.0

"""
    MIN_VALID_WAVELENGTH_M

Lower bound (m) of the wavelength range the material models are valid for.
"""
const MIN_VALID_WAVELENGTH_M = 1300e-9

"""
    MAX_VALID_WAVELENGTH_M

Upper bound (m) of the wavelength range the material models are valid for.
"""
const MAX_VALID_WAVELENGTH_M = 1700e-9

"""
    validate_molar_fraction(x) -> Float64

Return `x` as a `Float64`, throwing an `ArgumentError` unless it is a finite
dopant molar fraction in `[0, 1]`.
"""
function validate_molar_fraction(x::Real)
    xf = Float64(x)
    if !(isfinite(xf) && 0.0 <= xf <= 1.0)
        throw(ArgumentError("molar fraction must be between 0 and 1 inclusive"))
    end
    return xf
end

# TODO: Edit these functions to make them material-specific, with args for max/min
"""
    validate_model_temperature(T_K)

Return `T_K`, throwing an `ArgumentError` unless it is a finite positive
temperature (K) inside `[MIN_VALID_TEMPERATURE_K, MAX_VALID_TEMPERATURE_K]`.
"""
function validate_model_temperature(T_K)
    if !(isfinite(T_K) && T_K > 0.0)
        throw(ArgumentError("temperature must be a finite positive value in kelvin"))
    end
    if !(MIN_VALID_TEMPERATURE_K <= T_K <= MAX_VALID_TEMPERATURE_K)
        throw(ArgumentError(
            "temperature is outside the current model validity range " *
            "[$(MIN_VALID_TEMPERATURE_K), $(MAX_VALID_TEMPERATURE_K)] K: got $(T_K)"
        ))
    end
    return T_K
end

"""
    validate_model_wavelength(λ)

Return `float(λ)`, throwing an `ArgumentError` unless it is a finite positive
wavelength (m) inside `[MIN_VALID_WAVELENGTH_M, MAX_VALID_WAVELENGTH_M]`.
"""
function validate_model_wavelength(λ)
    λ = float(λ)
    if !(isfinite(λ) && λ > 0.0)
        throw(ArgumentError("wavelength must be a finite positive value in meters"))
    end
    if !(MIN_VALID_WAVELENGTH_M <= λ <= MAX_VALID_WAVELENGTH_M)
        throw(ArgumentError(
            "wavelength is outside the current model validity range " *
            "[$(MIN_VALID_WAVELENGTH_M), $(MAX_VALID_WAVELENGTH_M)] m: got $(λ)"
        ))
    end
    return λ
end

#################################################
#
# Common interpolation utilities
#
#################################################

"""
    interpolate_scalar(a, b, x)

Return the linear interpolation `(1 - x) * a + x * b`.

Used to interpolate scalar properties between two endmember glasses by dopant
molar fraction `x`.
"""
interpolate_scalar(a, b, x) = (one(x) - x) * a + x * b

"""
    interpolate_pair(a::Tuple, b::Tuple, x)

Return the elementwise linear interpolation of two 2-tuples.

Used for paired properties such as the photoelastic constants `(p11, p12)`.
"""
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
    _evaluate_sellmeier_polynomials(B_coeffs, C_coeffs, x) -> coeffs

Return an n-tuple of `(B, C)` Sellmeier pairs by evaluating per-resonance
polynomials in `x` (e.g. temperature or dopant fraction). See `silica.jl` for an
example of use.
"""
function _evaluate_sellmeier_polynomials(B_coeffs, C_coeffs, x)
    return map((B, C) -> (evalpoly(x, B), evalpoly(x, C)), B_coeffs, C_coeffs)
end

"""
    _evaluate_sellmeier_constants(coeffs, x) -> coeffs

Return the `(B, C)` Sellmeier pairs of `coeffs` unchanged in value, scaled by
`one(x)` so the element type lifts through `x` (e.g. `Particles`). See
`germania.jl` for an example of use.
"""
function _evaluate_sellmeier_constants(coeffs, x)
    scale = one(x)
    return map(coeffs) do (B, C)
        (B * scale, C * scale)
    end
end

"""
    sellmeier_index_from_coefficients(coeffs, λ)

Return the Sellmeier refractive index at wavelength `λ` (m).

The index follows ``n^2 = 1 + \\sum_i B_i \\lambda^2 / (\\lambda^2 - C_i^2)``,
where ``B_i`` and ``C_i`` are the strength and wavelength (μm) of each resonance.
`coeffs` is an n-tuple of `(B, C)` 2-tuples. `λ` is validated against the model
wavelength range.

A material whose coefficients are known can implement its index directly:

    refractive_index(::ValueOnly, material::YourMaterial, λ, T_K) =
        sellmeier_index_from_coefficients(YOUR_COEFFICIENTS, λ)
    refractive_index(::WithDerivative, material::YourMaterial, λ, T_K) =
        sellmeier_index_from_coefficients_dω(YOUR_COEFFICIENTS, λ)

Coefficients that depend on temperature or dopant fraction can be produced with
[`_evaluate_sellmeier_polynomials`](@ref) or [`_evaluate_sellmeier_constants`](@ref).
The arithmetic lifts through `MonteCarloMeasurements.Particles`.
"""
function sellmeier_index_from_coefficients(coeffs, λ)
    λ_m = validate_model_wavelength(λ)
    λ_um = λ_m * 1e6
    total = one(λ_um)
    for (B, C) in coeffs
        total += B * λ_um^2 / (λ_um^2 - C^2)
    end
    return sqrt(total)
end

"""
    sellmeier_index_from_coefficients_dω(coeffs, λ) -> SpectralResponse

Return the Sellmeier refractive index at wavelength `λ` (m) together with its
angular-frequency derivative.

See [`sellmeier_index_from_coefficients`](@ref) for the coefficient convention.
"""
function sellmeier_index_from_coefficients_dω(coeffs, λ)
    λ_m = validate_model_wavelength(λ)
    λ_um = λ_m * 1e6
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
    dλ_dω = -(λ_m^2) / (2π * SPEED_OF_LIGHT_M_PER_S)
    return SpectralResponse(n, dn_dλ * dλ_dω)
end
