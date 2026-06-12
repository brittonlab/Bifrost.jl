"""
Material properties for fluorinated silica glass.
File name prefixed with 20 so it is loaded after 10-silica-germania-binary glasses, as it
depends on PURE_SILICA and other siilica material constants defined there.

This file defines fluorinated silica glass, defined by a doping fraction
of fluorine into a silica base.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- x_f                   dopant molar fraction (dimensionless, 0..1)
- refractive indices, Poisson ratio, photoelastic constants: dimensionless
- cte                   1/K
- softening_temperature  K
- youngs_modulus        Pa
- nonlinear_refractive_index (n_2)  m²/W

[Example usage]

glass = SilicaFluorinatedGlass(0.01)   # 1.0 mol% F in SiO2
T_K = 297.15
λ = 1550e-9
n = refractive_index(glass, λ, T_K)
cte_value = cte(glass, T_K)
"""

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

# Fluorine enters as a molar fraction over the physical interval [0, 1], matching the
# binary-glass coefficient source's fraction validation.
const FLUORINE_MOLAR_FRACTION_RANGE = ValidityRange(0.0, 1.0, "fluorine molar fraction")

struct SilicaFluorinatedGlass <: AbstractMaterial
    x_f::Float64

    SilicaFluorinatedGlass(x_f::Real) =
        new(check_range(Float64(x_f), FLUORINE_MOLAR_FRACTION_RANGE))
end

# Fluorine doping does not shift the validity window inherited from pure silica.
runtime_ranges(::SilicaFluorinatedGlass) = runtime_ranges(PURE_SILICA)

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
    check_range((; T_K, λ), runtime_ranges(material))
    return sellmeier_index_from_coefficients(_sellmeier_coefficients(material, T_K), λ)
end

function refractive_index(::WithDerivative, material::SilicaFluorinatedGlass, λ, T_K)
    check_range((; T_K, λ), runtime_ranges(material))
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

function unsupported_fluorine_property(name::AbstractString)
    msg = "$(name) is not defined for fluorine-doped silica in the current model"
    throw(ArgumentError(msg))
end
