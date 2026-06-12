"""
Material properties for silica-germania binary glasses.

This file defines mixtures of pure silica SiO2 and pure germania GeO2
defined by a doping fraction of germania into a silica base.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- x_ge                  germania molar fraction (dimensionless, 0..0.05)
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa
- nonlinear_refractive_index (n_2)  m²/W

[Example usage]

glass = SilicaGermaniaGlass(0.036)   # 3.6 mol% GeO2 in SiO2
T_K = 297.15
λ = 1550e-9
n = refractive_index(glass, λ, T_K)
cte_value = cte(glass, T_K)
"""

# Caution: the validity range is an estimate; proper ranges tracked in issue #4.
const GERMANIA_FRACTION_RANGE = ValidRange(0.0, 0.05, "germania molar fraction")

struct SilicaGermaniaGlass <: AbstractMaterial
    x_ge::Float64

    SilicaGermaniaGlass(x_ge::Real) =
        new(_check_range(Float64(x_ge), GERMANIA_FRACTION_RANGE))
end

runtime_range(::SilicaGermaniaGlass) = runtime_range((SiO2(), GeO2()))

#################################################
#
# Refractive Index
#
#################################################

function refractive_index(::ValueOnly, glass::SilicaGermaniaGlass, λ, T_K)
    _check_range((; T_K, λ), runtime_range(glass))
    n_silica = refractive_index(ValueOnly(), SiO2(), λ, T_K)
    n_germania = refractive_index(ValueOnly(), GeO2(), λ, T_K)
    return _interpolate_scalar(n_silica, n_germania, glass.x_ge)
end

function refractive_index(::WithDerivative, glass::SilicaGermaniaGlass, λ, T_K)
    _check_range((; T_K, λ), runtime_range(glass))
    n_silica = refractive_index(WithDerivative(), SiO2(), λ, T_K)
    n_germania = refractive_index(WithDerivative(), GeO2(), λ, T_K)
    return SpectralResponse(
        _interpolate_scalar(n_silica.value, n_germania.value, glass.x_ge),
        _interpolate_scalar(n_silica.dω, n_germania.dω, glass.x_ge)
    )
end

#################################################
#
# Other Material Properties
#
#################################################

cte(glass::SilicaGermaniaGlass, _) = _interpolate_scalar(SILICA_CTE, GERMANIA_CTE, glass.x_ge)

softening_temperature(glass::SilicaGermaniaGlass, _) = _interpolate_scalar(SILICA_SOFTENING_TEMPERATURE_K, GERMANIA_SOFTENING_TEMPERATURE_K, glass.x_ge)

poisson_ratio(glass::SilicaGermaniaGlass, _) = _interpolate_scalar(SILICA_POISSON_RATIO, GERMANIA_POISSON_RATIO, glass.x_ge)

photoelastic_constants(glass::SilicaGermaniaGlass, _) = _interpolate_pair(SILICA_PHOTOELASTIC_CONSTANTS, GERMANIA_PHOTOELASTIC_CONSTANTS, glass.x_ge)

youngs_modulus(glass::SilicaGermaniaGlass, _) = _interpolate_scalar(SILICA_YOUNGS_MODULUS, GERMANIA_YOUNGS_MODULUS, glass.x_ge)