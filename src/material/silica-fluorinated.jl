"""
Material properties for fluorinated silica glass.

This file defines fluorinated silica glass, defined by a doping fraction
of fluorine into a silica base.

Units (SI unless noted):
- λ                     wavelength in m
- T_K                   temperature in K
- x_f                   fluorine molar fraction (dimensionless, 0..0.05)
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

# Caution: the validity range is an estimate.
const FLUORINE_MOLAR_FRACTION_RANGE = ValidRange(0.0, 0.05, "fluorine molar fraction")

struct SilicaFluorinatedGlass <: AbstractMaterial
    x_f::Float64

    SilicaFluorinatedGlass(x_f::Real) =
        new(_check_range(Float64(x_f), FLUORINE_MOLAR_FRACTION_RANGE))
end

runtime_range(::SilicaFluorinatedGlass) = runtime_range(SiO2())

function _sellmeier_coefficients(glass::SilicaFluorinatedGlass, T_K)
    silica_coeffs = _sellmeier_coefficients(SiO2(), T_K)
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
    _check_range((; T_K, λ), SILICA_VALIDITY)
    return _sellmeier_index_from_coefficients(_sellmeier_coefficients(material, T_K), λ)
end

function refractive_index(::WithDerivative, material::SilicaFluorinatedGlass, λ, T_K)
    _check_range((; T_K, λ), SILICA_VALIDITY)
    return _sellmeier_index_from_coefficients_dω(_sellmeier_coefficients(material, T_K), λ)
end

cte(::SilicaFluorinatedGlass, _) = _unsupported_fluorine_property("cte")

softening_temperature(::SilicaFluorinatedGlass, _) =
    _unsupported_fluorine_property("softening_temperature")

poisson_ratio(::SilicaFluorinatedGlass, _) = _unsupported_fluorine_property("poisson_ratio")

photoelastic_constants(::SilicaFluorinatedGlass, _) =
    _unsupported_fluorine_property("photoelastic_constants")

youngs_modulus(::SilicaFluorinatedGlass, _) =
    _unsupported_fluorine_property("youngs_modulus")

function _unsupported_fluorine_property(name::AbstractString)
    msg = "$(name) is not defined for fluorine-doped silica in the current model"
    throw(ArgumentError(msg))
end
