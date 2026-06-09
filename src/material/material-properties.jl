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
# Validation utilities
#
#################################################

const MIN_VALID_TEMPERATURE_K = 243.0
const MAX_VALID_TEMPERATURE_K = 373.0
const MIN_VALID_WAVELENGTH_M = 1300e-9
const MAX_VALID_WAVELENGTH_M = 1700e-9

function validate_molar_fraction(x::Real)
    xf = Float64(x)
    if !(isfinite(xf) && 0.0 <= xf <= 1.0)
        throw(ArgumentError("molar fraction must be between 0 and 1 inclusive"))
    end
    return xf
end

# TODO: Edit these functions to make them material-specific, with args for max/min
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
    λ_m = validate_model_wavelength(λ)
    λ_um = λ_m * 1e6
    total = one(λ_um)
    for (B, C) in coeffs
        total += B * λ_um^2 / (λ_um^2 - C^2)
    end
    return sqrt(total)
end

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
