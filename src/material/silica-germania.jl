"""
Material properties for silica-germania binary glasses.

This file defines mixtures of pure silica SiO2 and pure germania GeO2
defined by a doping fraction of germania into a silica base.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- x_ge                  dopant molar fraction (dimensionless, 0..1)
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

#################################################
#
# Material constants:
# All from !silica.jl and !germania.jl
#
#################################################

#################################################
#
# Structures and Utility Methods
#
#################################################

struct SilicaGermaniaGlass <: AbstractMaterial
    x_ge::Float64
    function SilicaGermaniaGlass(x_ge::Real)
        xf = validate_molar_fraction(x_ge)
        return new(xf)
    end
end

#################################################
#
# Refractive Index
#
#################################################

function refractive_index(::ValueOnly, glass::SilicaGermaniaGlass, λ, T_K)
    n_silica = refractive_index(ValueOnly(), PURE_SILICA, λ, T_K)
    n_germania = refractive_index(ValueOnly(), PURE_GERMANIA, λ, T_K)
    return interpolate_scalar(n_silica, n_germania, glass.x_ge)
end

function refractive_index(::WithDerivative, glass::SilicaGermaniaGlass, λ, T_K)
    n_silica = refractive_index(WithDerivative(), PURE_SILICA, λ, T_K)
    n_germania = refractive_index(WithDerivative(), PURE_GERMANIA, λ, T_K)
    return SpectralResponse(
        interpolate_scalar(n_silica.value, n_germania.value, glass.x_ge),
        interpolate_scalar(n_silica.dω, n_germania.dω, glass.x_ge)
    )
end

#################################################
#
# Other Material Properties
#
#################################################

cte(glass::SilicaGermaniaGlass, _) = interpolate_scalar(SILICA_CTE, GERMANIA_CTE, glass.x_ge)

softening_temperature(glass::SilicaGermaniaGlass, _) = interpolate_scalar(SILICA_SOFTENING_TEMPERATURE_K, GERMANIA_SOFTENING_TEMPERATURE_K, glass.x_ge)

poisson_ratio(glass::SilicaGermaniaGlass, _) = interpolate_scalar(SILICA_POISSON_RATIO, GERMANIA_POISSON_RATIO, glass.x_ge)

photoelastic_constants(glass::SilicaGermaniaGlass, _) = interpolate_pair(SILICA_PHOTOELASTIC_CONSTANTS, GERMANIA_PHOTOELASTIC_CONSTANTS, glass.x_ge)

youngs_modulus(glass::SilicaGermaniaGlass, _) = interpolate_scalar(SILICA_YOUNGS_MODULUS, GERMANIA_YOUNGS_MODULUS, glass.x_ge)

function nonlinear_refractive_index(glass::SilicaGermaniaGlass, λ, T_K)
    validate_model_wavelength(λ)
    validate_model_temperature(T_K)
    return interpolate_scalar(SILICA_N2, GERMANIA_N2, glass.x_ge)
end