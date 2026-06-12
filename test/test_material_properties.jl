using Test
using Bifrost

const MP = Bifrost.MaterialProperties

const SILICA_B_COEFFS = (
    (1.10127, -4.94251e-5, 5.27414e-7, -1.59700e-9, 1.75949e-12),
    (1.78752e-5, 4.76391e-5, -4.49019e-7, 1.44546e-9, -1.57223e-12),
    (7.93552e-1, -1.27815e-3, 1.84595e-5, -9.20275e-8, 1.48829e-10)
)

const SILICA_C_COEFFS = (
    (-8.906e-2, 9.0873e-6, -6.53638e-8, 7.77072e-11, 6.84605e-14),
    (2.97562e-1, -8.59578e-4, 6.59069e-6, -1.09482e-8, 7.85145e-13),
    (9.34454, -7.09788e-3, 1.01968e-4, -5.07660e-7, 8.21348e-10)
)

const GERMANIA_BC = (
    (0.80686642, 0.068972606),
    (0.71815848, 0.15396605),
    (0.85416831, 11.841931)
)

const FLUORINE_BC = (
    ((-61.25, 0.2565), (-23.0, 0.101)),
    ((73.9, -1.836), (10.7, -0.005)),
    ((233.5, -5.82), (1090.5, -24.695))
)

reference_polynomial(coeffs, T_kelvin) =
    sum(coeffs[i] * T_kelvin^(i - 1) for i in eachindex(coeffs))
reference_constant(value, _) = value
reference_quadratic(quadratic, linear, x) = quadratic * x^2 + linear * x
reference_scalar_mix(a, b, x) = (1 - x) * a + x * b
reference_pair_mix(a, b, x) =
    (reference_scalar_mix(a[1], b[1], x), reference_scalar_mix(a[2], b[2], x))

function reference_sellmeier_index(coeffs, λ_meters)
    λ_um = λ_meters * 1e6
    total = 1.0
    for (B, C) in coeffs
        total += B * λ_um^2 / (λ_um^2 - C^2)
    end
    return sqrt(total)
end

function reference_silica_coefficients(T_kelvin)
    return ntuple(i -> (
        reference_polynomial(SILICA_B_COEFFS[i], T_kelvin),
        reference_polynomial(SILICA_C_COEFFS[i], T_kelvin)
    ), 3)
end

reference_silica_index(λ_meters, T_kelvin) =
    reference_sellmeier_index(reference_silica_coefficients(T_kelvin), λ_meters)

function reference_germania_base_coefficients(T_kelvin)
    return ntuple(i -> (
        reference_constant(GERMANIA_BC[i][1], T_kelvin),
        reference_constant(GERMANIA_BC[i][2], T_kelvin)
    ), 3)
end

function reference_germania_thermo_optic_shift(T_kelvin)
    Tref = 297.15
    return 6.2153e-13 / 4 * (T_kelvin^4 - Tref^4) -
           5.3387e-10 / 3 * (T_kelvin^3 - Tref^3) +
           1.6654e-7 / 2 * (T_kelvin^2 - Tref^2)
end

function reference_germania_index(λ_meters, T_kelvin)
    n_ref = reference_sellmeier_index(
        reference_germania_base_coefficients(T_kelvin),
        λ_meters
    )
    return n_ref + reference_germania_thermo_optic_shift(T_kelvin)
end

function reference_fluorinated_coefficients(x_f, T_kelvin)
    silica = reference_silica_coefficients(T_kelvin)
    return ntuple(i -> begin
        B0, C0 = silica[i]
        ΔB = reference_quadratic(FLUORINE_BC[i][1][1], FLUORINE_BC[i][1][2], x_f)
        ΔC = reference_quadratic(FLUORINE_BC[i][2][1], FLUORINE_BC[i][2][2], x_f)
        (B0 + ΔB, C0 + ΔC)
    end, 3)
end

reference_fluorinated_index(λ_meters, T_kelvin, x_f) =
    reference_sellmeier_index(reference_fluorinated_coefficients(x_f, T_kelvin), λ_meters)

function unsupported_message(f)
    try
        f()
        return nothing
    catch err
        return sprint(showerror, err)
    end
end

reference_dλ_dω(λ_meters) = -(λ_meters^2) / (2π * SPEED_OF_LIGHT_M_PER_S)

function finite_difference_dω(f, λ_meters; dλ = 1e-12)
    df_dλ = (f(λ_meters + dλ) - f(λ_meters - dλ)) / (2dλ)
    return df_dλ * reference_dλ_dω(λ_meters)
end

@testset "Constructors and Sellmeier internals" begin
    @test SilicaGermaniaGlass(0.0).x_ge == 0.0
    @test SilicaGermaniaGlass(1.0).x_ge == 1.0
    @test SilicaFluorinatedGlass(0.0).x_f == 0.0
    @test SilicaFluorinatedGlass(1.0).x_f == 1.0

    @test_throws ArgumentError SilicaGermaniaGlass(-1e-6)
    @test_throws ArgumentError SilicaGermaniaGlass(1.000001)
    @test_throws ArgumentError SilicaFluorinatedGlass(-1e-6)
    @test_throws ArgumentError SilicaFluorinatedGlass(1.000001)

    exported_names = Set(names(MP))
    old_sellmeier_helpers = (
        :TemperaturePolynomial,
        :SellmeierConstantLaw,
        :SellmeierQuadraticMolarLaw,
        :SellmeierTerm,
        :SellmeierCorrectionTerm,
        :evaluate,
        :sellmeier_coefficients,
        :SILICA_TERM_1,
        :GERMANIA_TERM_1,
        :FLUORINE_TERM_1,
    )
    # T-GUARDRAIL: Sellmeier coefficient representation is private data, not API.
    for name in old_sellmeier_helpers
        @test !(name in exported_names)
    end
end

@testset "SiO2" begin
    silica = SiO2()

    # T-VALIDATION: Leviton and Frey fused-silica Sellmeier data,
    # doi:10.1117/12.672853.
    @test SILICA_C_COEFFS[3][2] == -7.09788e-3

    for T in (243.0, 297.15, 373.0)
        actual = MP._sellmeier_coefficients(silica, T)
        expected = reference_silica_coefficients(T)
        @test length(actual) == 3
        for i in 1:3
            @test actual[i][1] ≈ expected[i][1] atol = 1e-12 rtol = 1e-12
            @test actual[i][2] ≈ expected[i][2] atol = 1e-12 rtol = 1e-12
        end
    end

    for (λ, T) in ((1300e-9, 243.0), (1550e-9, 297.15), (1700e-9, 373.0))
        @test refractive_index(silica, λ, T) ≈
              reference_silica_index(λ, T) atol = 1e-12 rtol = 1e-12
    end

    λs = range(1300e-9, 1700e-9; length = 4)
    Ts = range(243.0, 373.0; length = 4)
    for T in Ts, λ in λs
        @test isfinite(refractive_index(silica, λ, T))
    end
end

@testset "GeO2" begin
    germania = GeO2()

    for (λ, T) in ((1300e-9, 243.0), (1550e-9, 297.15), (1700e-9, 373.0))
        @test refractive_index(germania, λ, T) ≈
              reference_germania_index(λ, T) atol = 1e-12 rtol = 1e-12
    end

    λ = 1550e-9
    n_ref = refractive_index(germania, λ, 297.15)
    n_warm = refractive_index(germania, λ, 350.0)
    n_ref_expected = reference_sellmeier_index(
        reference_germania_base_coefficients(297.15),
        λ
    )
    @test n_ref ≈ n_ref_expected atol = 1e-12 rtol = 1e-12
    @test reference_germania_thermo_optic_shift(297.15) == 0.0
    @test n_warm > n_ref

    @test cte(germania, 297.15) == GERMANIA_CTE
    @test softening_temperature(germania, 297.15) == GERMANIA_SOFTENING_TEMPERATURE_K
    @test poisson_ratio(germania, 297.15) == GERMANIA_POISSON_RATIO
    @test photoelastic_constants(germania, 297.15) == GERMANIA_PHOTOELASTIC_CONSTANTS
    @test youngs_modulus(germania, 297.15) == GERMANIA_YOUNGS_MODULUS
end

@testset "SilicaGermaniaGlass" begin
    λ = 1550e-9
    T = 297.15

    pure_silica_glass = SilicaGermaniaGlass(0.0)
    pure_germania_glass = SilicaGermaniaGlass(1.0)
    @test refractive_index(pure_silica_glass, λ, T) == refractive_index(SiO2(), λ, T)
    @test refractive_index(pure_germania_glass, λ, T) == refractive_index(GeO2(), λ, T)
    @test cte(pure_silica_glass, T) == cte(SiO2(), T)
    @test cte(pure_germania_glass, T) == cte(GeO2(), T)
    @test softening_temperature(pure_silica_glass, T) == softening_temperature(SiO2(), T)
    @test softening_temperature(pure_germania_glass, T) == softening_temperature(GeO2(), T)
    @test poisson_ratio(pure_silica_glass, T) == poisson_ratio(SiO2(), T)
    @test poisson_ratio(pure_germania_glass, T) == poisson_ratio(GeO2(), T)
    @test photoelastic_constants(pure_silica_glass, T) == photoelastic_constants(SiO2(), T)
    @test photoelastic_constants(pure_germania_glass, T) ==
          photoelastic_constants(GeO2(), T)
    @test youngs_modulus(pure_silica_glass, T) == youngs_modulus(SiO2(), T)
    @test youngs_modulus(pure_germania_glass, T) == youngs_modulus(GeO2(), T)

    for x in (0.036, 0.25)
        glass = SilicaGermaniaGlass(x)
        n_expected = reference_scalar_mix(
            reference_silica_index(λ, T),
            reference_germania_index(λ, T),
            x
        )
        @test refractive_index(glass, λ, T) ≈ n_expected atol = 1e-12 rtol = 1e-12
        @test cte(glass, T) == reference_scalar_mix(SILICA_CTE, GERMANIA_CTE, x)
        @test softening_temperature(glass, T) ==
              reference_scalar_mix(
                  SILICA_SOFTENING_TEMPERATURE_K,
                  GERMANIA_SOFTENING_TEMPERATURE_K,
                  x
              )
        @test poisson_ratio(glass, T) ==
              reference_scalar_mix(SILICA_POISSON_RATIO, GERMANIA_POISSON_RATIO, x)
        @test photoelastic_constants(glass, T) ==
              reference_pair_mix(
                  SILICA_PHOTOELASTIC_CONSTANTS,
                  GERMANIA_PHOTOELASTIC_CONSTANTS,
                  x
              )
        @test youngs_modulus(glass, T) ==
              reference_scalar_mix(SILICA_YOUNGS_MODULUS, GERMANIA_YOUNGS_MODULUS, x)
    end

    glass0 = SilicaGermaniaGlass(0.0)
    glass_small = SilicaGermaniaGlass(0.036)
    glass_large = SilicaGermaniaGlass(1.0)
    @test refractive_index(glass_small, λ, T) > refractive_index(glass0, λ, T)
    @test cte(glass0, T) < cte(glass_small, T) < cte(glass_large, T)
    @test poisson_ratio(glass0, T) <
          poisson_ratio(glass_small, T) <
          poisson_ratio(glass_large, T)

    p11, p12 = photoelastic_constants(SilicaGermaniaGlass(0.036), T)
    @test p11 isa Float64
    @test p12 isa Float64
    @test p11 < p12
end

@testset "SilicaFluorinatedGlass" begin
    λ = 1550e-9
    T = 297.15

    undoped = SilicaFluorinatedGlass(0.0)
    @test refractive_index(undoped, λ, T) ≈
          refractive_index(SiO2(), λ, T) atol = 1e-12 rtol = 1e-12

    for x in (0.01, 0.02)
        glass = SilicaFluorinatedGlass(x)
        actual = MP._sellmeier_coefficients(glass, T)
        expected = reference_fluorinated_coefficients(x, T)
        for i in 1:3
            @test actual[i][1] ≈ expected[i][1] atol = 1e-12 rtol = 1e-12
            @test actual[i][2] ≈ expected[i][2] atol = 1e-12 rtol = 1e-12
        end
        @test refractive_index(glass, λ, T) ≈
              reference_fluorinated_index(λ, T, x) atol = 1e-12 rtol = 1e-12
    end

    @test refractive_index(SilicaFluorinatedGlass(0.01), λ, T) <
          refractive_index(SiO2(), λ, T)

    unsupported_calls = [
        () -> cte(SilicaFluorinatedGlass(0.01), T),
        () -> softening_temperature(SilicaFluorinatedGlass(0.01), T),
        () -> poisson_ratio(SilicaFluorinatedGlass(0.01), T),
        () -> photoelastic_constants(SilicaFluorinatedGlass(0.01), T),
        () -> youngs_modulus(SilicaFluorinatedGlass(0.01), T),
    ]

    for f in unsupported_calls
        @test_throws ArgumentError f()
        msg = unsupported_message(f)
        @test !isnothing(msg)
        @test occursin("fluorine-doped silica", lowercase(msg))
    end
end

@testset "Validity guards" begin
    silica = SiO2()
    germania = GeO2()
    ge_glass = SilicaGermaniaGlass(0.036)
    f_glass = SilicaFluorinatedGlass(0.01)

    for T in (243.0, 373.0)
        @test isfinite(refractive_index(silica, 1550e-9, T))
        @test isfinite(refractive_index(germania, 1550e-9, T))
        @test isfinite(refractive_index(ge_glass, 1550e-9, T))
        @test isfinite(refractive_index(f_glass, 1550e-9, T))
    end

    for λ in (1300e-9, 1700e-9)
        @test isfinite(refractive_index(silica, λ, 297.15))
        @test isfinite(refractive_index(germania, λ, 297.15))
        @test isfinite(refractive_index(ge_glass, λ, 297.15))
        @test isfinite(refractive_index(f_glass, λ, 297.15))
    end

    @test_throws ArgumentError refractive_index(silica, 1550e-9, 242.999)
    @test_throws ArgumentError refractive_index(germania, 1550e-9, 373.001)
    @test_throws ArgumentError refractive_index(ge_glass, 1550e-9, 242.999)
    @test_throws ArgumentError refractive_index(f_glass, 1550e-9, 373.001)

    @test_throws ArgumentError refractive_index(silica, 1299.999e-9, 297.15)
    @test_throws ArgumentError refractive_index(germania, 1700.001e-9, 297.15)
    @test_throws ArgumentError refractive_index(ge_glass, 1299.999e-9, 297.15)
    @test_throws ArgumentError refractive_index(f_glass, 1700.001e-9, 297.15)
end

@testset "Validity range primitives" begin
    # T-GUARDRAIL: a range checks one value, returning it in range and throwing outside.
    r = ValidityRange(1.0, 2.0, "demo")
    @test check_range(1.5, r) == 1.5
    @test check_range(1.0, r) == 1.0      # closed lower bound
    @test check_range(2.0, r) == 2.0      # closed upper bound
    @test_throws ArgumentError check_range(0.999, r)
    @test_throws ArgumentError check_range(2.001, r)
    @test_throws ArgumentError check_range(NaN, r)
    @test_throws ArgumentError check_range(Inf, r)

    # T-GUARDRAIL: the NamedTuple method validates each value against its keyed range.
    ranges = (a = ValidityRange(0.0, 1.0, "a"), b = ValidityRange(10.0, 20.0, "b"))
    @test check_range((a = 0.5, b = 15.0), ranges) === nothing
    @test_throws ArgumentError check_range((a = 0.5, b = 25.0), ranges)

    # T-GUARDRAIL: declared runtime windows carry the documented bounds.
    for m in (SiO2(), GeO2(), SilicaGermaniaGlass(0.036), SilicaFluorinatedGlass(0.01))
        rr = runtime_range(m)
        @test rr.T_K.lo == 243.0 && rr.T_K.hi == 373.0
        @test rr.λ.lo == 1300e-9 && rr.λ.hi == 1700e-9
    end
end

@testset "Spectral responses" begin
    λ = 1550e-9
    T = 297.15

    materials = (
        SiO2(),
        GeO2(),
        SilicaGermaniaGlass(0.036),
        SilicaFluorinatedGlass(0.01)
    )

    for material in materials
        scalar = refractive_index(material, λ, T)
        resp = refractive_index(WithDerivative(), material, λ, T)
        @test resp isa SpectralResponse
        @test resp.value ≈ scalar atol = 1e-14 rtol = 1e-14
        dω_expected = finite_difference_dω(λp -> refractive_index(material, λp, T), λ)
        @test resp.dω ≈ dω_expected atol = 1e-18 rtol = 1e-6
    end
end
