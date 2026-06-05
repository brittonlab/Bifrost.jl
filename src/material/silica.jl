"""
Material properties for pure silica SiO2 glass.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa
- nonlinear_refractive_index (n_2)  m²/W

[Example usage]

glass = SiO2()
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

# Sellmeier polynomial coefficients from Leviton and Frey,
# Optomechanical Technologies for Astronomy 2006.

const SILICA_TERM_1 = SellmeierTerm(
    T -> 1.10127 - 4.94251e-5*T + 5.27414e-7*T^2 - 1.59700e-9*T^3 + 1.75949e-12*T^4,
    T -> -8.906e-2 + 9.0873e-6*T - 6.53638e-8*T^2 + 7.77072e-11*T^3 + 6.84605e-14*T^4
)

const SILICA_TERM_2 = SellmeierTerm(
    T -> 1.78752e-5 + 4.76391e-5*T - 4.49019e-7*T^2 + 1.44546e-9*T^3 - 1.57223e-12*T^4,
    T -> 2.97562e-1 - 8.59578e-4*T + 6.59069e-6*T^2 - 1.09482e-8*T^3 + 7.85145e-13*T^4
)

const SILICA_TERM_3 = SellmeierTerm(
    T -> 7.93552e-1 - 1.27815e-3*T + 1.84595e-5*T^2 - 9.20275e-8*T^3 + 1.48829e-10*T^4,
    T -> 9.34454 - 7.09788e-3*T + 1.01968e-4*T^2 - 5.07660e-7*T^3 + 8.21348e-10*T^4
)

const SILICA_CTE = 5.4e-7

const SILICA_SOFTENING_TEMPERATURE_K = 1100.0 + 273.15

const SILICA_POISSON_RATIO = 0.170

const SILICA_PHOTOELASTIC_CONSTANTS = (0.121, 0.270)

const SILICA_YOUNGS_MODULUS = 74e9

const SILICA_N2 = 2.2e-20

#################################################
#
# Structures and Utility Methods
#
#################################################

struct SiO2 <: AbstractMaterial
    sellmeier_terms::NTuple{3, SellmeierTerm}
end

const PURE_SILICA = SiO2((SILICA_TERM_1, SILICA_TERM_2, SILICA_TERM_3))
SiO2() = PURE_SILICA

#################################################
#
# Refractive Index
#
#################################################

# Override sellmeier_coefficients_generic to do validation of T_K
function sellmeier_coefficients(material::SiO2, T_K)
    T = validate_model_temperature(T_K)
    return map(term -> evaluate(term, T), material.sellmeier_terms)
end

refractive_index(::ValueOnly, material::SiO2, λ, T_K) = 
        sellmeier_index_from_coefficients(sellmeier_coefficients(material, T_K), λ)
refractive_index(::WithDerivative, material::SiO2, λ, T_K) = 
        sellmeier_index_from_coefficients_dω(sellmeier_coefficients(material, T_K), λ)

#################################################
#
# Other Material Properties
#
#################################################

cte(::SiO2, _) = SILICA_CTE
softening_temperature(::SiO2, _) = SILICA_SOFTENING_TEMPERATURE_K
poisson_ratio(::SiO2, _) = SILICA_POISSON_RATIO
photoelastic_constants(::SiO2, _) = SILICA_PHOTOELASTIC_CONSTANTS
youngs_modulus(::SiO2, _) = SILICA_YOUNGS_MODULUS

function nonlinear_refractive_index(::SiO2, λ, T_K)
    validate_model_wavelength(λ)
    validate_model_temperature(T_K)
    return SILICA_N2
end