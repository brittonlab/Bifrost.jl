# Fluorine-doped silica. Loaded after silica.jl, whose PURE_SILICA Sellmeier
# coefficients are the base that the fluorine corrections perturb.

const _FLUORINE_SELLMEIER_B_CORRECTION_COEFFS = (
    (0.0, 0.2565, -61.25),
    (0.0, -1.836, 73.9),
    (0.0, -5.82, 233.5)
)

const _FLUORINE_SELLMEIER_C_CORRECTION_COEFFS = (
    (0.0, 0.101, -23.0),
    (0.0, -0.005, 10.7),
    (0.0, -24.695, 1090.5)
)

"""
    SilicaFluorinatedGlass(x_f)

Fluorine-doped silica glass with fluorine molar fraction `x_f`.

The refractive index applies fluorine-dependent corrections to the pure-silica
Sellmeier coefficients; `x_f` is validated to lie in `[0, 1]`. Only the optical
properties are modeled: the mechanical and thermal properties (`cte`,
`softening_temperature`, `poisson_ratio`, `photoelastic_constants`,
`youngs_modulus`) throw an `ArgumentError` in the current model.

# Examples
```julia
glass = SilicaFluorinatedGlass(0.01)   # 1.0 mol% F in SiO2
n = refractive_index(glass, 1550e-9, 297.15)
```
"""
struct SilicaFluorinatedGlass <: AbstractMaterial
    x_f::Float64

    function SilicaFluorinatedGlass(x_f::Real)
        xf = validate_molar_fraction(x_f)
        return new(xf)
    end
end

function _sellmeier_coefficients(glass::SilicaFluorinatedGlass, T_K)
    silica_coeffs = _sellmeier_coefficients(PURE_SILICA, T_K)
    x_f = glass.x_f
    corrections = _evaluate_sellmeier_polynomials(
        _FLUORINE_SELLMEIER_B_CORRECTION_COEFFS,
        _FLUORINE_SELLMEIER_C_CORRECTION_COEFFS,
        x_f
    )
    return map(silica_coeffs, corrections) do (B0, C0), (ΔB, ΔC)
        (B0 + ΔB, C0 + ΔC)
    end
end

function refractive_index(::ValueOnly, material::SilicaFluorinatedGlass, λ, T_K)
    return sellmeier_index_from_coefficients(_sellmeier_coefficients(material, T_K), λ)
end

function refractive_index(::WithDerivative, material::SilicaFluorinatedGlass, λ, T_K)
    return sellmeier_index_from_coefficients_dω(_sellmeier_coefficients(material, T_K), λ)
end

cte(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("cte")

softening_temperature(::SilicaFluorinatedGlass, _) =
    unsupported_fluorine_property("softening_temperature")

poisson_ratio(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("poisson_ratio")

photoelastic_constants(::SilicaFluorinatedGlass, _) =
    unsupported_fluorine_property("photoelastic_constants")

youngs_modulus(::SilicaFluorinatedGlass, _) =
    unsupported_fluorine_property("youngs_modulus")

"""
    unsupported_fluorine_property(name)

Throw the `ArgumentError` reporting that property `name` is not defined for
fluorine-doped silica in the current model.
"""
function unsupported_fluorine_property(name::AbstractString)
    msg = "$(name) is not defined for fluorine-doped silica in the current model"
    throw(ArgumentError(msg))
end
