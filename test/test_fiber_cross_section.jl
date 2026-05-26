using Test
using Bifrost

"""
Why these tests are currently weak.
- They do not provide independent evidence from literature or fiber.py.
- They do not protect well against copy-pasted formula errors.
- They can give a false sense of coverage because they exercise the same algebra twice.
What could do better.
- Compare against literature values for a specific fiber, wavelength, and temperature.
- Compare against fiber.py for several specific fibers, wavelengths, and temperatures.
- Scaling tests like bend birefringence ~ 1/R^2, axial-tension term ~ tension/R, 
  twist birefringence linear in twist rate.
- Independently derived dispersion checks using a different numerical route than the 
  implementation.
"""


const TEST_C = 299_792_458.0
const TEST_V_CUTOFF = 2.405

reference_indices(fiber, λ_m, T_K) = (
    refractive_index(fiber.core_material, λ_m, T_K),
    refractive_index(fiber.cladding_material, λ_m, T_K)
)

reference_k0(λ_m) = 2π / λ_m

function reference_v(fiber, λ_m, T_K)
    n_core, n_clad = reference_indices(fiber, λ_m, T_K)
    return core_radius(fiber) * reference_k0(λ_m) * sqrt(n_core^2 - n_clad^2)
end

function reference_beta(fiber, λ_m, T_K)
    n_core, _ = reference_indices(fiber, λ_m, T_K)
    v = reference_v(fiber, λ_m, T_K)
    waveguide_term = ((1 + sqrt(2)) * v) / (1 + (4 + v^4)^(1 / 4))
    return sqrt((n_core^2) * reference_k0(λ_m)^2 - waveguide_term^2 / core_radius(fiber)^2)
end

reference_neff(fiber, λ_m, T_K) = reference_beta(fiber, λ_m, T_K) / reference_k0(λ_m)

function reference_ng(fiber, λ_m, T_K; dλ = 0.1e-9)
    n_minus = reference_neff(fiber, λ_m - dλ, T_K)
    n_plus = reference_neff(fiber, λ_m + dλ, T_K)
    return reference_neff(fiber, λ_m, T_K) - λ_m * (n_plus - n_minus) / (2dλ)
end

function reference_aeff(fiber, λ_m, T_K)
    v = reference_v(fiber, λ_m, T_K)
    w_over_r = 0.65 + 1.619 / v^1.5 + 2.879 / v^6
    return π * (w_over_r * core_radius(fiber))^2
end

function reference_gamma(fiber, λ_m, T_K)
    n2 = nonlinear_refractive_index(fiber.core_material, λ_m, T_K)
    return reference_k0(λ_m) * n2 / reference_aeff(fiber, λ_m, T_K)
end

function reference_d(fiber, λ_m, T_K; dλ = 0.1e-9)
    n_minus = reference_neff(fiber, λ_m - dλ, T_K)
    n_center = reference_neff(fiber, λ_m, T_K)
    n_plus = reference_neff(fiber, λ_m + dλ, T_K)
    return -λ_m / TEST_C * (n_plus - 2n_center + n_minus) / dλ^2 * 1e6
end

function reference_beta2(fiber, λ_m, T_K; dλ = 0.1e-9)
    return -(λ_m^2 / (2π * TEST_C)) * (reference_d(fiber, λ_m, T_K; dλ = dλ) * 1e-6)
end

function reference_eccentricity_squared(axis_ratio; signed = false)
    ε = Float64(axis_ratio)
    if ε >= 1.0
        return 1 - inv(ε)^2
    end
    value = 1 - ε^2
    return signed ? -value : value
end

function reference_core_noncircularity_birefringence(fiber, λ_m, T_K; axis_ratio)
    n_core, n_clad = reference_indices(fiber, λ_m, T_K)
    v = reference_v(fiber, λ_m, T_K)
    e2 = reference_eccentricity_squared(axis_ratio; signed = true)
    return (e2 * (1 - n_clad^2 / n_core^2)^(3 / 2)) / core_radius(fiber) *
           (4 / v^3) * log(v)^3 / (1 + log(v))
end

function reference_asymmetric_thermal_stress_birefringence(fiber, λ_m, T_K; axis_ratio)
    ε = Float64(axis_ratio)
    n_core, _ = reference_indices(fiber, λ_m, T_K)
    β = reference_beta(fiber, λ_m, T_K)
    v = reference_v(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    α_core = cte(fiber.core_material, T_K)
    α_clad = cte(fiber.cladding_material, T_K)
    T_soft = softening_temperature(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    k0 = reference_k0(λ_m)
    prefactor = k0 * (1 - core_radius(fiber)^2 * (n_core^2 * k0^2 - β^2) / v^2)
    stress_factor = 0.5 * n_core^3 * (p11 - p12) * (α_clad - α_core) *
                    abs(T_soft - T_K) / (1 - ν^2)
    return prefactor * stress_factor * ((ε - 1) / (ε + 1))
end

function reference_bending_birefringence(fiber, λ_m, T_K; bend_radius_m)
    n_core, _ = reference_indices(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    r = cladding_radius(fiber)
    return reference_k0(λ_m) * (n_core^3 / 2) * (p11 - p12) * (1 + ν) * 0.5 * (r^2 / bend_radius_m^2)
end

function reference_axial_tension_birefringence(fiber, λ_m, T_K; bend_radius_m, axial_tension_N)
    n_core, _ = reference_indices(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    ν = poisson_ratio(fiber.core_material, T_K)
    E = youngs_modulus(fiber.core_material, T_K)
    r = cladding_radius(fiber)
    tension_term = ((2 - 3ν) / (1 - ν)) * (r / bend_radius_m) * (axial_tension_N / (π * r^2 * E))
    return reference_k0(λ_m) * (n_core^3 / 2) * (p11 - p12) * (1 + ν) * tension_term
end

function reference_twisting_birefringence(fiber, λ_m, T_K; twist_rate_rad_per_m)
    n_core, _ = reference_indices(fiber, λ_m, T_K)
    p11, p12 = photoelastic_constants(fiber.core_material, T_K)
    return (1 + (n_core^2 / 2) * (p11 - p12)) * twist_rate_rad_per_m
end

reference_dλ_dω(λ_meters) = -(λ_meters^2) / (2π * TEST_C)

function finite_difference_dω(f, λ_meters; dλ = 1e-12)
    df_dλ = (f(λ_meters + dλ) - f(λ_meters - dλ)) / (2dλ)
    return df_dλ * reference_dλ_dω(λ_meters)
end

const SMF_LIKE_FIBER = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

@testset "StepIndexCrossSection reference formulas" begin
    fiber = SMF_LIKE_FIBER
    λ = 1550e-9
    T = 297.15

    @test fiber.manufacturer == "Corning"
    @test fiber.model_number == "SMF-like"

    n_core, n_clad = reference_indices(fiber, λ, T)
    @test n_core > n_clad

    @test numerical_aperture(fiber, λ, T) ≈ sqrt(n_core^2 - n_clad^2) atol = 1e-14
    @test relative_index_difference(fiber, λ, T) ≈ (n_core - n_clad) / n_clad atol = 1e-14
    @test normalized_frequency(fiber, λ, T) ≈ reference_v(fiber, λ, T) rtol = 1e-12
    @test propagation_constant(fiber, λ, T) ≈ reference_beta(fiber, λ, T) rtol = 1e-12
    @test effective_mode_index(fiber, λ, T) ≈ reference_neff(fiber, λ, T) rtol = 1e-12
    @test effective_group_index(fiber, λ, T) ≈ reference_ng(fiber, λ, T) rtol = 1e-10
    @test effective_mode_area(fiber, λ, T) ≈ reference_aeff(fiber, λ, T) rtol = 1e-12
    @test nonlinear_coefficient(fiber, λ, T) ≈ reference_gamma(fiber, λ, T) rtol = 1e-12
    @test chromatic_dispersion_parameter(fiber, λ, T) ≈ reference_d(fiber, λ, T) rtol = 1e-10
    @test group_velocity_dispersion_parameter(fiber, λ, T) ≈ reference_beta2(fiber, λ, T) rtol = 1e-10

    # SMF-like sanity checks inspired by the operating regime discussed in the paper.
    @test 1.5 < normalized_frequency(fiber, λ, T) < TEST_V_CUTOFF
    @test 70e-12 < effective_mode_area(fiber, λ, T) < 95e-12
    @test 8.0 < chromatic_dispersion_parameter(fiber, λ, T) < 20.0
    @test 9e-4 < nonlinear_coefficient(fiber, λ, T) < 1.5e-3
end

@testset "StepIndexCrossSection cutoff behavior" begin
    fiber = SMF_LIKE_FIBER
    T = 297.15
    λ_cutoff = cutoff_wavelength(fiber, T)

    @test normalized_frequency(fiber, λ_cutoff, T) ≈ TEST_V_CUTOFF atol = 1e-6
    @test !is_single_mode(fiber, λ_cutoff - 1e-9, T)
    @test is_single_mode(fiber, λ_cutoff + 1e-9, T)
end

@testset "StepIndexCrossSection spectral responses" begin
    fiber = SMF_LIKE_FIBER
    λ = 1550e-9
    T = 297.15

    n_resp = guided_refractive_indices(WithDerivative(), fiber, λ, T)
    @test n_resp[1] isa SpectralResponse
    @test n_resp[2] isa SpectralResponse
    @test n_resp[1].value ≈ core_refractive_index(fiber, λ, T) atol = 1e-14
    @test n_resp[2].value ≈ cladding_refractive_index(fiber, λ, T) atol = 1e-14

    V_resp = normalized_frequency(WithDerivative(), fiber, λ, T)
    β_resp = propagation_constant(WithDerivative(), fiber, λ, T)
    neff_resp = effective_mode_index(WithDerivative(), fiber, λ, T)
    @test V_resp.value ≈ normalized_frequency(fiber, λ, T) atol = 1e-14
    @test β_resp.value ≈ propagation_constant(fiber, λ, T) atol = 1e-14
    @test neff_resp.value ≈ effective_mode_index(fiber, λ, T) atol = 1e-14
    @test V_resp.dω ≈ finite_difference_dω(λp -> normalized_frequency(fiber, λp, T), λ) atol = 1e-16 rtol = 1e-6
    @test β_resp.dω ≈ finite_difference_dω(λp -> propagation_constant(fiber, λp, T), λ) atol = 1e-12 rtol = 1e-6
    @test neff_resp.dω ≈ finite_difference_dω(λp -> effective_mode_index(fiber, λp, T), λ) atol = 1e-16 rtol = 1e-6

    bend_resp = bending_birefringence(WithDerivative(), fiber, λ, T; bend_radius_m = 0.03)
    twist_resp = twisting_birefringence(WithDerivative(), fiber, λ, T; twist_rate_rad_per_m = 0.7)
    cnc_resp = core_noncircularity_birefringence(WithDerivative(), fiber, λ, T; axis_ratio = 1.01)
    @test bend_resp.Δβ ≈ bending_birefringence(fiber, λ, T; bend_radius_m = 0.03) atol = 1e-14
    @test twist_resp.Δβ ≈ twisting_birefringence(fiber, λ, T; twist_rate_rad_per_m = 0.7) atol = 1e-14
    @test cnc_resp.Δβ ≈ core_noncircularity_birefringence(fiber, λ, T; axis_ratio = 1.01) atol = 1e-14
    @test bend_resp.dω ≈ finite_difference_dω(λp -> bending_birefringence(fiber, λp, T; bend_radius_m = 0.03), λ) atol = 1e-18 rtol = 1e-6
    @test twist_resp.dω ≈ finite_difference_dω(λp -> twisting_birefringence(fiber, λp, T; twist_rate_rad_per_m = 0.7), λ) atol = 1e-18 rtol = 1e-6
    @test cnc_resp.dω ≈ finite_difference_dω(λp -> core_noncircularity_birefringence(fiber, λp, T; axis_ratio = 1.01), λ) atol = 1e-12 rtol = 1e-6
end

@testset "StepIndexCrossSection perturbation formulas and invariants" begin
    fiber = SMF_LIKE_FIBER
    λ = 1550e-9
    T = 297.15
    ε = 1.01
    R = 0.03
    tf = 0.5
    tr = 0.7

    @test core_noncircularity_birefringence(fiber, λ, T; axis_ratio = 1.0) == 0.0
    @test asymmetric_thermal_stress_birefringence(fiber, λ, T; axis_ratio = 1.0) == 0.0
    @test bending_birefringence(fiber, λ, T; bend_radius_m = Inf) == 0.0
    @test axial_tension_birefringence(fiber, λ, T; bend_radius_m = Inf, axial_tension_N = tf) == 0.0
    @test axial_tension_birefringence(fiber, λ, T; bend_radius_m = R, axial_tension_N = 0.0) == 0.0
    @test twisting_birefringence(fiber, λ, T; twist_rate_rad_per_m = 0.0) == 0.0

    @test core_noncircularity_birefringence(fiber, λ, T; axis_ratio = ε) ≈
          reference_core_noncircularity_birefringence(fiber, λ, T; axis_ratio = ε) rtol = 1e-12
    @test asymmetric_thermal_stress_birefringence(fiber, λ, T; axis_ratio = ε) ≈
          reference_asymmetric_thermal_stress_birefringence(fiber, λ, T; axis_ratio = ε) rtol = 1e-12
    @test bending_birefringence(fiber, λ, T; bend_radius_m = R) ≈
          reference_bending_birefringence(fiber, λ, T; bend_radius_m = R) rtol = 1e-12
    @test axial_tension_birefringence(fiber, λ, T; bend_radius_m = R, axial_tension_N = tf) ≈
          reference_axial_tension_birefringence(fiber, λ, T; bend_radius_m = R, axial_tension_N = tf) rtol = 1e-12
    @test twisting_birefringence(fiber, λ, T; twist_rate_rad_per_m = tr) ≈
          reference_twisting_birefringence(fiber, λ, T; twist_rate_rad_per_m = tr) rtol = 1e-12

    @test core_noncircularity_birefringence(fiber, λ, T; axis_ratio = inv(ε)) ≈
          -core_noncircularity_birefringence(fiber, λ, T; axis_ratio = ε) rtol = 1e-12
    @test asymmetric_thermal_stress_birefringence(fiber, λ, T; axis_ratio = inv(ε)) ≈
          -asymmetric_thermal_stress_birefringence(fiber, λ, T; axis_ratio = ε) rtol = 1e-12
end

@testset "StepIndexCrossSection guided and fluorinated edge cases" begin
    λ = 1550e-9
    T = 297.15

    unguided = StepIndexCrossSection(
        SilicaFluorinatedGlass(0.005),
        SilicaGermaniaGlass(0.0),
        8.2e-6,
        125e-6
    )
    @test_throws ArgumentError normalized_frequency(unguided, λ, T)

    fluorinated_cladding = StepIndexCrossSection(
        SilicaGermaniaGlass(0.036),
        SilicaFluorinatedGlass(0.005),
        8.2e-6,
        125e-6
    )

    @test isfinite(normalized_frequency(fluorinated_cladding, λ, T))
    @test core_noncircularity_birefringence(fluorinated_cladding, λ, T; axis_ratio = 1.0) +
          asymmetric_thermal_stress_birefringence(fluorinated_cladding, λ, T; axis_ratio = 1.0) +
          bending_birefringence(fluorinated_cladding, λ, T; bend_radius_m = Inf) +
          axial_tension_birefringence(
              fluorinated_cladding,
              λ,
              T;
              bend_radius_m = Inf,
              axial_tension_N = 0.0
          ) +
          twisting_birefringence(fluorinated_cladding, λ, T; twist_rate_rad_per_m = 0.0) == 0.0
    @test_throws ArgumentError asymmetric_thermal_stress_birefringence(
        fluorinated_cladding,
        λ,
        T;
        axis_ratio = 1.01
    )
end
