# Binary silica-germania glasses. All material constants come from silica.jl
# and germania.jl.

"""
    SilicaGermaniaGlass(x_ge)

Binary silica–germania glass with germania molar fraction `x_ge`.

Every property is the linear interpolation between [`PURE_SILICA`](@ref
PURE_SILICA) and [`PURE_GERMANIA`](@ref PURE_GERMANIA) at fraction `x_ge`
(validated to lie in `[0, 1]`). Implements the [`AbstractMaterial`](@ref)
property interface (SI units).

# Examples
```julia
glass = SilicaGermaniaGlass(0.036)   # 3.6 mol% GeO2 in SiO2
n = refractive_index(glass, 1550e-9, 297.15)
α = cte(glass, 297.15)
```
"""
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

cte(glass::SilicaGermaniaGlass, _) =
    interpolate_scalar(SILICA_CTE, GERMANIA_CTE, glass.x_ge)

softening_temperature(glass::SilicaGermaniaGlass, _) = interpolate_scalar(
    SILICA_SOFTENING_TEMPERATURE_K, GERMANIA_SOFTENING_TEMPERATURE_K, glass.x_ge)

poisson_ratio(glass::SilicaGermaniaGlass, _) =
    interpolate_scalar(SILICA_POISSON_RATIO, GERMANIA_POISSON_RATIO, glass.x_ge)

photoelastic_constants(glass::SilicaGermaniaGlass, _) = interpolate_pair(
    SILICA_PHOTOELASTIC_CONSTANTS, GERMANIA_PHOTOELASTIC_CONSTANTS, glass.x_ge)

youngs_modulus(glass::SilicaGermaniaGlass, _) =
    interpolate_scalar(SILICA_YOUNGS_MODULUS, GERMANIA_YOUNGS_MODULUS, glass.x_ge)