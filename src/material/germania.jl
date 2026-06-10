"""
Material properties for pure germania GeO2 glass.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa
- nonlinear_refractive_index (n_2)  m²/W

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

# Sellmeier coefficients from Fleming, Applied Optics (1984).
# doi:10.1364/AO.23.004486
const _GERMANIA_SELLMEIER_COEFFICIENTS = (
    (0.80686642, 0.068972606),
    (0.71815848, 0.15396605),
    (0.85416831, 11.841931)
)

const GERMANIA_REFERENCE_TEMPERATURE_K = 297.15

const GERMANIA_CTE = 10e-6

const GERMANIA_SOFTENING_TEMPERATURE_K = 300.0 + 273.15

const GERMANIA_POISSON_RATIO = 0.212

const GERMANIA_PHOTOELASTIC_CONSTANTS = (0.130, 0.288)

const GERMANIA_YOUNGS_MODULUS = 45.5e9

const GERMANIA_N2 = 4.6e-20

#################################################
#
# Structures and Utility Methods
#
#################################################

"""
    GeO2()

Pure germania (GeO₂) glass.

Implements the [`AbstractMaterial`](@ref) property interface (SI units). The
refractive index combines room-temperature Sellmeier coefficients from Fleming,
doi:10.1364/AO.23.004486, with the thermo-optic shift of
[`thermo_optic_index_shift`](@ref).

# Examples
```julia
glass = GeO2()
n = refractive_index(glass, 1550e-9, 297.15)
α = cte(glass, 297.15)
```
"""
struct GeO2 <: AbstractMaterial end

"""
    PURE_GERMANIA

Shared singleton `GeO2()` instance, used as the dopant endmember of
silica–germania mixtures.
"""
const PURE_GERMANIA = GeO2()

#################################################
#
# Pure Germania GeO2 Refractive Index
#
#################################################

"""
    thermo_optic_index_shift(material::GeO2, T_K)

Return the thermo-optic refractive-index shift of germania relative to the
reference temperature `GERMANIA_REFERENCE_TEMPERATURE_K`.

Polynomial model from G. M. Rego, Sensors (2024), doi:10.3390/s24154857.
"""
function thermo_optic_index_shift(material::GeO2, T_K)
    T = validate_model_temperature(T_K)
    Tref = GERMANIA_REFERENCE_TEMPERATURE_K
    return 6.2153e-13 / 4 * (T^4 - Tref^4) -
           5.3387e-10 / 3 * (T^3 - Tref^3) +
           1.6654e-7 / 2 * (T^2 - Tref^2)
end

function _sellmeier_coefficients(::GeO2, T_K)
    T = validate_model_temperature(T_K)
    return _evaluate_sellmeier_constants(_GERMANIA_SELLMEIER_COEFFICIENTS, T)
end

function refractive_index(::ValueOnly, material::GeO2, λ, T_K)
    T = validate_model_temperature(T_K)
    base_coeffs = _sellmeier_coefficients(material, T)
    n_ref = sellmeier_index_from_coefficients(base_coeffs, λ)
    return n_ref + thermo_optic_index_shift(material, T)
end

function refractive_index(::WithDerivative, material::GeO2, λ, T_K)
    T = validate_model_temperature(T_K)
    base_coeffs = _sellmeier_coefficients(material, T)
    base = sellmeier_index_from_coefficients_dω(base_coeffs, λ)
    return SpectralResponse(base.value + thermo_optic_index_shift(material, T), base.dω)
end

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
