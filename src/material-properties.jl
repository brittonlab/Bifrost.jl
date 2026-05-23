"""
Material properties for silica-based optical glasses.

This file defines endpoint materials (`SiO2`, `GeO2`) and binary silica-glass
models with either germania or fluorine as the single dopant.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- x_ge, x_f             dopant molar fraction (dimensionless, 0..1)
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa
- nonlinear_refractive_index (n2)  m²/W
"""

# Example usage
"""
glass = GermaniaSilicaGlass(0.036)   # 3.6 mol% GeO2 in SiO2
T_K = 297.15
λ = 1550e-9
n = refractive_index(glass, λ, T_K)
cte_value = cte(glass, T_K)
"""

abstract type AbstractMaterial end
abstract type AbstractBinarySilicaGlass <: AbstractMaterial end
abstract type SpectralStyle end

struct ValueOnly <: SpectralStyle end
struct WithDerivative <: SpectralStyle end

struct SpectralResponse{T}
    value::T
    dω::T
end

const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0

struct TemperaturePolynomial
    coeffs::NTuple{5, Float64}
    TemperaturePolynomial(coeffs::NTuple{5, <:Real}) = new(map(Float64, coeffs))
end

function (poly::TemperaturePolynomial)(T_K)
    value = zero(T_K)
    power = one(T_K)
    for coeff in poly.coeffs
        value += coeff * power
        power *= T_K
    end
    return value
end

struct SellmeierConstantLaw
    value::Float64
    SellmeierConstantLaw(value::Real) = new(Float64(value))
end

(law::SellmeierConstantLaw)(x) = law.value * one(x)

struct SellmeierQuadraticMolarLaw
    quadratic::Float64
    linear::Float64
    SellmeierQuadraticMolarLaw(quadratic::Real, linear::Real) = new(Float64(quadratic), Float64(linear))
end

(law::SellmeierQuadraticMolarLaw)(x) = law.quadratic * x^2 + law.linear * x

struct SellmeierTerm{TB, TC}
    B_law::TB
    C_law::TC
end

evaluate(term::SellmeierTerm, temperature_like) = (term.B_law(temperature_like), term.C_law(temperature_like))

struct SellmeierCorrectionTerm
    ΔB_law::SellmeierQuadraticMolarLaw
    ΔC_law::SellmeierQuadraticMolarLaw
end

evaluate(term::SellmeierCorrectionTerm, molar_fraction) = (term.ΔB_law(molar_fraction), term.ΔC_law(molar_fraction))

struct SiO2 <: AbstractMaterial
    sellmeier_terms::NTuple{3, SellmeierTerm}
end

struct GeO2 <: AbstractMaterial
    sellmeier_terms::NTuple{3, SellmeierTerm}
end

struct GermaniaSilicaGlass <: AbstractBinarySilicaGlass
    x_ge::Float64
    function GermaniaSilicaGlass(x_ge::Real)
        xf = validate_molar_fraction(x_ge)
        return new(xf)
    end
end

struct FluorinatedSilicaGlass <: AbstractBinarySilicaGlass
    x_f::Float64
    function FluorinatedSilicaGlass(x_f::Real)
        xf = validate_molar_fraction(x_f)
        return new(xf)
    end
end

const SILICA_TERM_1 = SellmeierTerm(
    TemperaturePolynomial((1.10127, -4.94251e-5, 5.27414e-7, -1.59700e-9, 1.75949e-12)),
    TemperaturePolynomial((-8.906e-2, 9.0873e-6, -6.53638e-8, 7.77072e-11, 6.84605e-14))
)

const SILICA_TERM_2 = SellmeierTerm(
    TemperaturePolynomial((1.78752e-5, 4.76391e-5, -4.49019e-7, 1.44546e-9, -1.57223e-12)),
    TemperaturePolynomial((2.97562e-1, -8.59578e-4, 6.59069e-6, -1.09482e-8, 7.85145e-13))
)

const SILICA_TERM_3 = SellmeierTerm(
    TemperaturePolynomial((7.93552e-1, -1.27815e-3, 1.84595e-5, -9.20275e-8, 1.48829e-10)),
    TemperaturePolynomial((9.34454, -70.9788e-3, 1.01968e-4, -5.07660e-7, 8.21348e-10))
)

const GERMANIA_TERM_1 = SellmeierTerm(SellmeierConstantLaw(0.80686642), SellmeierConstantLaw(0.068972606))
const GERMANIA_TERM_2 = SellmeierTerm(SellmeierConstantLaw(0.71815848), SellmeierConstantLaw(0.15396605))
const GERMANIA_TERM_3 = SellmeierTerm(SellmeierConstantLaw(0.85416831), SellmeierConstantLaw(11.841931))

const FLUORINE_TERM_1 = SellmeierCorrectionTerm(
    SellmeierQuadraticMolarLaw(-61.25, 0.2565),
    SellmeierQuadraticMolarLaw(-23.0, 0.101)
)

const FLUORINE_TERM_2 = SellmeierCorrectionTerm(
    SellmeierQuadraticMolarLaw(73.9, -1.836),
    SellmeierQuadraticMolarLaw(10.7, -0.005)
)

const FLUORINE_TERM_3 = SellmeierCorrectionTerm(
    SellmeierQuadraticMolarLaw(233.5, -5.82),
    SellmeierQuadraticMolarLaw(1090.5, -24.695)
)

const PURE_SILICA = SiO2((SILICA_TERM_1, SILICA_TERM_2, SILICA_TERM_3))
const PURE_GERMANIA = GeO2((GERMANIA_TERM_1, GERMANIA_TERM_2, GERMANIA_TERM_3))
const FLUORINE_CORRECTION_TERMS = (FLUORINE_TERM_1, FLUORINE_TERM_2, FLUORINE_TERM_3)

SiO2() = PURE_SILICA
GeO2() = PURE_GERMANIA

const GERMANIA_REFERENCE_TEMPERATURE_K = 297.15

const SILICA_CTE = 5.4e-7
const GERMANIA_CTE = 10e-6

const SILICA_SOFTENING_TEMPERATURE_K = 1100.0 + 273.15
const GERMANIA_SOFTENING_TEMPERATURE_K = 300.0 + 273.15

const SILICA_POISSON_RATIO = 0.170
const GERMANIA_POISSON_RATIO = 0.212

const SILICA_PHOTOELASTIC_CONSTANTS = (0.121, 0.270)
const GERMANIA_PHOTOELASTIC_CONSTANTS = (0.130, 0.288)

const SILICA_YOUNGS_MODULUS = 74e9
const GERMANIA_YOUNGS_MODULUS = 45.5e9

const SILICA_N2 = 2.2e-20
const GERMANIA_N2 = 4.6e-20

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

# used for AbstractBinarySilicaGlass 
interpolate_scalar(a, b, x) = (one(x) - x) * a + x * b

# used for photoelastic_constants
function interpolate_pair(a::Tuple, b::Tuple, x)
    return (
        interpolate_scalar(a[1], b[1], x),
        interpolate_scalar(a[2], b[2], x)
    )
end

function sellmeier_coefficients(material::SiO2, T_K)
    T = validate_model_temperature(T_K)
    return map(term -> evaluate(term, T), material.sellmeier_terms)
end

function sellmeier_coefficients(glass::FluorinatedSilicaGlass, T_K)
    silica_coeffs = sellmeier_coefficients(PURE_SILICA, T_K)
    x_f = glass.x_f
    return ntuple(i -> begin
        B0, C0 = silica_coeffs[i]
        ΔB, ΔC = evaluate(FLUORINE_CORRECTION_TERMS[i], x_f)
        (B0 + ΔB, C0 + ΔC)
    end, 3)
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

function thermo_optic_index_shift(material::GeO2, T_K)
    T = validate_model_temperature(T_K)
    Tref = GERMANIA_REFERENCE_TEMPERATURE_K
    return 6.2153e-13 / 4 * (T^4 - Tref^4) -
           5.3387e-10 / 3 * (T^3 - Tref^3) +
           1.6654e-7 / 2 * (T^2 - Tref^2)
end

function reference_refractive_index(material::GeO2, λ, T_K)
    T = validate_model_temperature(T_K)
    base_coeffs = map(term -> evaluate(term, T), material.sellmeier_terms)
    n_ref = sellmeier_index_from_coefficients(base_coeffs, λ)
    return n_ref + thermo_optic_index_shift(material, T)
end

function reference_refractive_index(::WithDerivative, material::GeO2, λ, T_K)
    T = validate_model_temperature(T_K)
    base_coeffs = map(term -> evaluate(term, T), material.sellmeier_terms)
    base = sellmeier_index_from_coefficients_dω(base_coeffs, λ)
    return SpectralResponse(base.value + thermo_optic_index_shift(material, T), base.dω)
end

function refractive_index(::ValueOnly, material::SiO2, λ, T_K)
    coeffs = sellmeier_coefficients(material, T_K)
    return sellmeier_index_from_coefficients(coeffs, λ)
end

function refractive_index(::WithDerivative, material::SiO2, λ, T_K)
    coeffs = sellmeier_coefficients(material, T_K)
    return sellmeier_index_from_coefficients_dω(coeffs, λ)
end

refractive_index(style::ValueOnly, material::GeO2, λ, T_K) =
    reference_refractive_index(material, λ, T_K)

refractive_index(style::WithDerivative, material::GeO2, λ, T_K) =
    reference_refractive_index(style, material, λ, T_K)

function refractive_index(::ValueOnly, glass::GermaniaSilicaGlass, λ, T_K)
    n_silica = refractive_index(ValueOnly(), PURE_SILICA, λ, T_K)
    n_germania = refractive_index(ValueOnly(), PURE_GERMANIA, λ, T_K)
    return interpolate_scalar(n_silica, n_germania, glass.x_ge)
end

function refractive_index(::WithDerivative, glass::GermaniaSilicaGlass, λ, T_K)
    n_silica = refractive_index(WithDerivative(), PURE_SILICA, λ, T_K)
    n_germania = refractive_index(WithDerivative(), PURE_GERMANIA, λ, T_K)
    return SpectralResponse(
        interpolate_scalar(n_silica.value, n_germania.value, glass.x_ge),
        interpolate_scalar(n_silica.dω, n_germania.dω, glass.x_ge)
    )
end

function refractive_index(::ValueOnly, glass::FluorinatedSilicaGlass, λ, T_K)
    coeffs = sellmeier_coefficients(glass, T_K)
    return sellmeier_index_from_coefficients(coeffs, λ)
end

function refractive_index(::WithDerivative, glass::FluorinatedSilicaGlass, λ, T_K)
    coeffs = sellmeier_coefficients(glass, T_K)
    return sellmeier_index_from_coefficients_dω(coeffs, λ)
end

refractive_index(material::AbstractMaterial, λ, T_K) =
    refractive_index(ValueOnly(), material, λ, T_K)

cte(::SiO2, _) = SILICA_CTE
cte(::GeO2, _) = GERMANIA_CTE
cte(glass::GermaniaSilicaGlass, _) = interpolate_scalar(SILICA_CTE, GERMANIA_CTE, glass.x_ge)
cte(::FluorinatedSilicaGlass, _) = unsupported_fluorine_property("cte")

softening_temperature(::SiO2, _) = SILICA_SOFTENING_TEMPERATURE_K
softening_temperature(::GeO2, _) = GERMANIA_SOFTENING_TEMPERATURE_K
softening_temperature(glass::GermaniaSilicaGlass, _) = interpolate_scalar(SILICA_SOFTENING_TEMPERATURE_K, GERMANIA_SOFTENING_TEMPERATURE_K, glass.x_ge)
softening_temperature(::FluorinatedSilicaGlass, _) = unsupported_fluorine_property("softening_temperature")

poisson_ratio(::SiO2, _) = SILICA_POISSON_RATIO
poisson_ratio(::GeO2, _) = GERMANIA_POISSON_RATIO
poisson_ratio(glass::GermaniaSilicaGlass, _) = interpolate_scalar(SILICA_POISSON_RATIO, GERMANIA_POISSON_RATIO, glass.x_ge)
poisson_ratio(::FluorinatedSilicaGlass, _) = unsupported_fluorine_property("poisson_ratio")

photoelastic_constants(::SiO2, _) = SILICA_PHOTOELASTIC_CONSTANTS
photoelastic_constants(::GeO2, _) = GERMANIA_PHOTOELASTIC_CONSTANTS
photoelastic_constants(glass::GermaniaSilicaGlass, _) = interpolate_pair(SILICA_PHOTOELASTIC_CONSTANTS, GERMANIA_PHOTOELASTIC_CONSTANTS, glass.x_ge)
photoelastic_constants(::FluorinatedSilicaGlass, _) = unsupported_fluorine_property("photoelastic_constants")

youngs_modulus(::SiO2, _) = SILICA_YOUNGS_MODULUS
youngs_modulus(::GeO2, _) = GERMANIA_YOUNGS_MODULUS
youngs_modulus(glass::GermaniaSilicaGlass, _) = interpolate_scalar(SILICA_YOUNGS_MODULUS, GERMANIA_YOUNGS_MODULUS, glass.x_ge)
youngs_modulus(::FluorinatedSilicaGlass, _) = unsupported_fluorine_property("youngs_modulus")

function nonlinear_refractive_index(::SiO2, λ, T_K)
    validate_model_wavelength(λ)
    validate_model_temperature(T_K)
    return SILICA_N2
end

function nonlinear_refractive_index(::GeO2, λ, T_K)
    validate_model_wavelength(λ)
    validate_model_temperature(T_K)
    return GERMANIA_N2
end

function nonlinear_refractive_index(glass::GermaniaSilicaGlass, λ, T_K)
    validate_model_wavelength(λ)
    validate_model_temperature(T_K)
    return interpolate_scalar(SILICA_N2, GERMANIA_N2, glass.x_ge)
end

nonlinear_refractive_index(::FluorinatedSilicaGlass, _, _) = unsupported_fluorine_property("nonlinear_refractive_index")

function unsupported_fluorine_property(name::AbstractString)
    throw(ArgumentError("$(name) is not defined for fluorine-doped silica in the current model"))
end
