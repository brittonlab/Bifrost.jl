"""
Material properties for pure germania GeO2 glass.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa

[Example usage]

glass = GeO2()
T_K = 297.15
λ = 1550e-9
n = refractive_index(glass, λ, T_K)
cte_value = cte(glass, T_K)
"""

#################################################
#
# Material constants
#
#################################################

const GERMANIA_TERM_1 = SellmeierTerm(SellmeierConstantLaw(0.80686642), SellmeierConstantLaw(0.068972606))
const GERMANIA_TERM_2 = SellmeierTerm(SellmeierConstantLaw(0.71815848), SellmeierConstantLaw(0.15396605))
const GERMANIA_TERM_3 = SellmeierTerm(SellmeierConstantLaw(0.85416831), SellmeierConstantLaw(11.841931))

const GERMANIA_REFERENCE_TEMPERATURE_K = 297.15

const GERMANIA_CTE = 10e-6

const GERMANIA_SOFTENING_TEMPERATURE_K = 300.0 + 273.15

const GERMANIA_POISSON_RATIO = 0.212

const GERMANIA_PHOTOELASTIC_CONSTANTS = (0.130, 0.288)

const GERMANIA_YOUNGS_MODULUS = 45.5e9

#################################################
#
# Structures and Utility Methods
#
#################################################

struct GeO2 <: AbstractMaterial
    sellmeier_terms::NTuple{3, SellmeierTerm}
end

const PURE_GERMANIA = GeO2((GERMANIA_TERM_1, GERMANIA_TERM_2, GERMANIA_TERM_3))
GeO2() = PURE_GERMANIA

#################################################
#
# Pure Germania GeO2 Refractive Index
#
#################################################

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

refractive_index(style::ValueOnly, material::GeO2, λ, T_K) =
    reference_refractive_index(material, λ, T_K)

refractive_index(style::WithDerivative, material::GeO2, λ, T_K) =
    reference_refractive_index(style, material, λ, T_K)

#################################################
#
# Other Material Properties
#
#################################################

cte(::GeO2, _) = GERMANIA_CTE
softening_temperature(::GeO2, _) = GERMANIA_SOFTENING_TEMPERATURE_K
poisson_ratio(::GeO2, _) = GERMANIA_POISSON_RATIO
photoelastic_constants(::GeO2, _) = GERMANIA_PHOTOELASTIC_CONSTANTS
youngs_modulus(::GeO2, _) = GERMANIA_YOUNGS_MODULUS