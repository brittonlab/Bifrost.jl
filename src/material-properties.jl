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
"""

const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0

abstract type AbstractMaterial end
abstract type SpectralStyle end

# dω is the ω-derivative of the function that produces value at the parameters that produce value
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

# TODO Banner... validate ranges
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
Every material must implement

    refractive_index(::ValueOnly, material, λ, T_K) -> Float64
    refractive_index(::WithDerivative, material, λ, T_K) -> SpectralResponse
    
That's it. Implement those however makes sense for your material.

Most materials are modeled with the Sellmeier equation, which states
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

To store Sellmeier coefficients and any dependence on any parameters, we provide the 
SellmeierTerm struct. It has two properties, B_law and C_law, which must be either numbers
or callables that depend on some parameter of interest. For example, you may define constants:

    germania_term_1 = SellmeierTerm(0.80686642, 0.068972606)

or functions of something, such as temperature:

    const SILICA_TERM_1 = SellmeierTerm(
        T -> 1.10127 - 4.94251e-5*T + 5.27414e-7*T^2 - 1.59700e-9*T^3 + 1.75949e-12*T^4,
        T -> -8.906e-2 + 9.0873e-6*T - 6.53638e-8*T^2 + 7.77072e-11*T^3 + 6.84605e-14*T^4
    )

One can then define a material by these terms, e.g.

    struct SiO2 <: AbstractMaterial
        sellmeier_terms::NTuple{3, SellmeierTerm}
    end
    const PURE_SILICA = SiO2((SILICA_TERM_1, SILICA_TERM_2, SILICA_TERM_3))

Then one can use the generic sellmeier_coefficients(material::AbstractMaterial, param) provided 
here to obtain the tuple of 2-tuples needed for sellmeier_index_from_coefficients(). NOTE: the
generic function requires `material` to have a field `sellmeier_terms``, and it also has no
validation for `param`` because it doesn't know what `param` is. It's strongly recommended to
override sellmeier_coefficients() to do your own parameter validation before calling the
map(term -> evaluate()) line on the terms.)

We also provide the SellmeierCorrectionTerm struct which is purely an API convenience with
the same behavior as SellmeierTerm but with different names. This reduces confusion for
materials whose refractive indices are defined by corrections to the indices of another material.

Note also that any implemented material must ensure compatibility with `Particles` to allow
Monte Carlo calculation. This happens naturally through the Sellmeier structure if all B and C
coefficients are numbers or callable polynomials. Anything more complicated than polynomials
needs to be carefully vetted with Particles.
"""

function _validate_law(law, name::String)
    # Accept anything callable
    if applicable(law, 273.15)
        return nothing
    end
    
    throw(ArgumentError(
        "$name must be a number or callable, got $(typeof(law))"
    ))
end

struct SellmeierTerm{TB, TC}
    B_law::TB
    C_law::TC

    function SellmeierTerm(B_law, C_law)
        B_law_norm = B_law isa Number ? (x -> B_law * one(x)) : B_law
        C_law_norm = C_law isa Number ? (x -> C_law * one(x)) : C_law
        _validate_law(B_law_norm, "B_law")
        _validate_law(C_law_norm, "C_law")
        return new{typeof(B_law_norm), typeof(C_law_norm)}(B_law_norm, C_law_norm)
    end
end

function _evaluate_law(law, param)
    return law isa Number ? law : law(param)
end

evaluate(term::SellmeierTerm, param) = (_evaluate_law(term.B_law, param), _evaluate_law(term.C_law, param))

struct SellmeierCorrectionTerm{TB, TC}
    ΔB_law::TB
    ΔC_law::TC

    function SellmeierCorrectionTerm(ΔB_law, ΔC_law)
        _validate_law(ΔB_law, "ΔB_law")
        _validate_law(ΔC_law, "ΔC_law")
        return new{typeof(ΔB_law), typeof(ΔC_law)}(ΔB_law, ΔC_law)
    end
end

evaluate(term::SellmeierCorrectionTerm, param) =
        (_evaluate_law(term.ΔB_law, param), _evaluate_law(term.ΔC_law, param))

function sellmeier_coefficients(material::AbstractMaterial, param)
    # NO validation of param! Validate param before calling this, or override this.
    # The material must have a sellmeier_terms field to use this function.
    if !hasfield(typeof(material), :sellmeier_terms)
        throw(ArgumentError(
            "$(typeof(material)) does not have a sellmeier_terms field."
        ))
    end
    return map(term -> evaluate(term, param), material.sellmeier_terms)
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

"""

Generic interface documentation:
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
