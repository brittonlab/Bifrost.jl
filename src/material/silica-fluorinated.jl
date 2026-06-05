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

[Example usage]

glass = SilicaFluorinatedGlass(0.01)   # 1.0 mol% F in SiO2
T_K = 297.15
λ = 1550e-9
n = refractive_index(glass, λ, T_K)
cte_value = cte(glass, T_K)
"""

const FLUORINE_TERM_1 = SellmeierCorrectionTerm(
    SellmeierQuadraticMolarLaw(-61.25, 0.2565),
    SellmeierQuadraticMolarLaw(-23.0, 0.101)
)

const FLUORINE_TERM_2 = SellmeierCorrectionTerm(
    SellmeierQuadraticMolarLaw(73.9, -1.836),
    SellmeierQuadraticMolarLaw(10.7, -0.005)
)

const FLUORINE_TERM_3 = SellmeierCorrectionTerm(
    SellmeierQuadraticMolarLaw(233.5, -5.82),
    SellmeierQuadraticMolarLaw(1090.5, -24.695)
)

const FLUORINE_CORRECTION_TERMS = (FLUORINE_TERM_1, FLUORINE_TERM_2, FLUORINE_TERM_3)

struct SilicaFluorinatedGlass <: AbstractMaterial
    x_f::Float64
    function SilicaFluorinatedGlass(x_f::Real)
        xf = validate_molar_fraction(x_f)
        return new(xf)
    end
end

function sellmeier_coefficients(glass::SilicaFluorinatedGlass, T_K)
    silica_coeffs = sellmeier_coefficients(PURE_SILICA, T_K)
    x_f = glass.x_f
    return ntuple(i -> begin
        B0, C0 = silica_coeffs[i]
        ΔB, ΔC = evaluate(FLUORINE_CORRECTION_TERMS[i], x_f)
        (B0 + ΔB, C0 + ΔC)
    end, 3)
end

function refractive_index(::ValueOnly, glass::SilicaFluorinatedGlass, λ, T_K)
    coeffs = sellmeier_coefficients(glass, T_K)
    return sellmeier_index_from_coefficients(coeffs, λ)
end

function refractive_index(::WithDerivative, glass::SilicaFluorinatedGlass, λ, T_K)
    coeffs = sellmeier_coefficients(glass, T_K)
    return sellmeier_index_from_coefficients_dω(coeffs, λ)
end

cte(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("cte")

softening_temperature(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("softening_temperature")

poisson_ratio(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("poisson_ratio")

photoelastic_constants(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("photoelastic_constants")

youngs_modulus(::SilicaFluorinatedGlass, _) = unsupported_fluorine_property("youngs_modulus")

function unsupported_fluorine_property(name::AbstractString)
    throw(ArgumentError("$(name) is not defined for fluorine-doped silica in the current model"))
end