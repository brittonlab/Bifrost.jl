#################################################
#
# Material constants
#
#################################################

# Sellmeier polynomial coefficients from Leviton and Frey, doi:10.1117/12.672853.
const _SILICA_SELLMEIER_B_COEFFS = (
    (1.10127, -4.94251e-5, 5.27414e-7, -1.59700e-9, 1.75949e-12),
    (1.78752e-5, 4.76391e-5, -4.49019e-7, 1.44546e-9, -1.57223e-12),
    (7.93552e-1, -1.27815e-3, 1.84595e-5, -9.20275e-8, 1.48829e-10)
)

const _SILICA_SELLMEIER_C_COEFFS = (
    (-8.906e-2, 9.0873e-6, -6.53638e-8, 7.77072e-11, 6.84605e-14),
    (2.97562e-1, -8.59578e-4, 6.59069e-6, -1.09482e-8, 7.85145e-13),
    (9.34454, -7.09788e-3, 1.01968e-4, -5.07660e-7, 8.21348e-10)
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

"""
    SiO2()

Pure fused-silica (SiO₂) glass.

Implements the [`AbstractMaterial`](@ref) property interface (SI units). The
refractive index uses temperature-dependent Sellmeier polynomial coefficients
from Leviton and Frey, doi:10.1117/12.672853.

# Examples
```julia
glass = SiO2()
n = refractive_index(glass, 1550e-9, 297.15)
α = cte(glass, 297.15)
```
"""
struct SiO2 <: AbstractMaterial end

"""
    PURE_SILICA

Shared singleton `SiO2()` instance, used as the base glass of doped mixtures.
"""
const PURE_SILICA = SiO2()

#################################################
#
# Refractive Index
#
#################################################

function _sellmeier_coefficients(::SiO2, T_K)
    T = validate_model_temperature(T_K)
    return _evaluate_sellmeier_polynomials(
        _SILICA_SELLMEIER_B_COEFFS,
        _SILICA_SELLMEIER_C_COEFFS,
        T
    )
end

function refractive_index(::ValueOnly, material::SiO2, λ, T_K)
    return sellmeier_index_from_coefficients(_sellmeier_coefficients(material, T_K), λ)
end

function refractive_index(::WithDerivative, material::SiO2, λ, T_K)
    return sellmeier_index_from_coefficients_dω(_sellmeier_coefficients(material, T_K), λ)
end

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
