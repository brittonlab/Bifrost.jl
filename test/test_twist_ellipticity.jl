using Bifrost
using Bifrost.PathGeometry: twist_rate, torsion_phase, spin_phase, twist_phase
using Test
using LinearAlgebra

# Physics validation for the spin / mechanical-twist / geometric-torsion model
# (#8) and the ellipticity + asymmetric-thermal-stress generators (#3), wired in
# twist-phase-d. Expectations are derived from the physics, not from code output.

const _XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0), 8.2e-6, 125e-6)
const _T_REF = 297.15
const _λ = 1550e-9

# Elliptical core sharing the same materials/diameters (axes along local N at φ).
_elliptical(; ratio = 1.05, angle = 0.0) = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0), 8.2e-6, 125e-6;
    ellipticity_axis_ratio = ratio, ellipticity_axis_angle = angle)

_straight_fiber(xs; length = 1.0, twist = nothing, spin_rate = nothing,
                meta = AbstractMeta[], T_ref_K = _T_REF) = begin
    sb = SubpathBuilder(); start!(sb; spin_rate = spin_rate)
    straight!(sb; length = length, twist = twist, meta = meta)
    seal!(sb)
    Fiber(sb; cross_section = xs, T_ref_K = T_ref_K)
end

# Extract the linear-birefringence eigen-axis angle φ (mod π) from a generator
# value K = 0.5i·Δβ·[[c2φ, s2φ],[s2φ, −c2φ]] via 2φ = atan2(Im K₁₂, Im K₁₁).
_linear_axis_angle(K) = 0.5 * atan(imag(K[1, 2]), imag(K[1, 1]))

@testset "T-PHYSICS: straight unspun circular fiber → identity Jones" begin
    f = _straight_fiber(_XS; length = 1.5)
    J, _ = propagate_fiber(f; λ_m = _λ, verbose = false)
    @test isapprox(J, Matrix{ComplexF64}(I, 2, 2); atol = 1e-9)
end

@testset "T-PHYSICS: pure manufacturing spin adds no birefringence" begin
    # A circular (unstressed) core spun at any rate stays isotropic: spin rotates
    # nonexistent linear axes, so the Jones matrix is the identity. (Under the
    # pre-#8 model spin was wrongly routed into the circular term.)
    f = _straight_fiber(_XS; length = 1.5, spin_rate = 12.0)
    J, _ = propagate_fiber(f; λ_m = _λ, verbose = false)
    @test isapprox(J, Matrix{ComplexF64}(I, 2, 2); atol = 1e-9)
end

@testset "T-PHYSICS: mechanical twist → optical rotation g·τ_m" begin
    # K = circular generator with rate Δβc = twisting_birefringence(τ_m); over a
    # straight length L the Jones matrix is a real rotation by θ = ½·Δβc·L.
    τm = 5.0
    L  = 2.0
    Δβc = twisting_birefringence(_XS, _λ, _T_REF; twist_rate_rad_per_m = τm)
    f = _straight_fiber(_XS; length = L, twist = τm)
    J, _ = propagate_fiber(f; λ_m = _λ, verbose = false)
    θ = 0.5 * Δβc * L
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    @test isapprox(J, R; atol = 1e-7)
end

@testset "T-PHYSICS: helix rotates the bend axis as ∫τ_g" begin
    # On a torsioned (helical) centerline the curvature-direction birefringence
    # axis rotates relative to the Bishop frame at the geometric torsion rate, so
    # its angle advances by τ_g·Δs between two points. No twist/spin ⇒ the
    # generator is purely the (linear) bend term.
    R = 0.05; pitch = 0.02; turns = 1.0
    h = pitch / (2π)
    τg = h / (R^2 + h^2)
    sb = SubpathBuilder(); start!(sb)
    helix!(sb; radius = R, pitch = pitch, turns = turns)
    seal!(sb)
    f = Fiber(sb; cross_section = _XS, T_ref_K = _T_REF)
    K = generator_K(f, _λ)
    s1 = 0.1 * pitch
    s2 = 0.4 * pitch
    Δφ = _linear_axis_angle(K(s2)) - _linear_axis_angle(K(s1))
    @test isapprox(Δφ, τg * (s2 - s1); atol = 1e-6)
end

@testset "T-PHYSICS: constant ellipticity → analytic linear retarder" begin
    # Straight elliptical fiber, axes along N (angle 0): K = diag(½iΔβ, −½iΔβ)
    # with Δβ = |core ellipticity| + |asymmetric thermal stress|. Over length L
    # the Jones matrix is diag(exp(½iΔβL), exp(−½iΔβL)).
    xs = _elliptical(; ratio = 1.05, angle = 0.0)
    L  = 1.0
    Δβ = core_noncircularity_birefringence(xs, _λ, _T_REF) +
         asymmetric_thermal_stress_birefringence(xs, _λ, _T_REF)
    @test Δβ > 0
    f = _straight_fiber(xs; length = L)
    J, _ = propagate_fiber(f; λ_m = _λ, verbose = false)
    expected = ComplexF64[exp(0.5im * Δβ * L) 0; 0 exp(-0.5im * Δβ * L)]
    @test isapprox(J, expected; atol = 1e-7)
end

@testset "T-PHYSICS: :T_K shifts thermal-stress birefringence (local temperature)" begin
    # The asymmetric-thermal-stress magnitude ∝ |T_soft − T|. A segment carrying
    # :T_K = ΔT must be evaluated at T = T_ref + ΔT, so its intrinsic-linear Δβ
    # equals the cross-section evaluated there — different from the T_ref value.
    xs = _elliptical(; ratio = 1.05, angle = 0.0)
    ΔT = 50.0
    L  = 1.0
    f = _straight_fiber(xs; length = L, meta = [MCMadd(:T_K, ΔT)])
    s = 0.1 * L                                   # inside the (expanded) segment
    @test isapprox(local_temperature(f, s), _T_REF + ΔT; atol = 1e-9)

    K = generator_K(f, _λ)
    Δβ_at_s = 2 * imag(K(s)[1, 1])                # axis angle 0 ⇒ c2φ = 1
    Δβ_hot = core_noncircularity_birefringence(xs, _λ, _T_REF + ΔT) +
             asymmetric_thermal_stress_birefringence(xs, _λ, _T_REF + ΔT)
    Δβ_ref = core_noncircularity_birefringence(xs, _λ, _T_REF) +
             asymmetric_thermal_stress_birefringence(xs, _λ, _T_REF)
    @test isapprox(Δβ_at_s, Δβ_hot; rtol = 1e-10)
    @test !isapprox(Δβ_at_s, Δβ_ref; rtol = 1e-6)
end

@testset "T-GUARDRAIL: circular core contributes no ellipticity birefringence" begin
    # ellipticity_axis_ratio == 1 ⇒ the intrinsic-linear term is exactly zero, so
    # a straight circular fiber has an all-zero generator everywhere.
    f = _straight_fiber(_XS; length = 1.0)
    K = generator_K(f, _λ)
    Kval = K(0.5)
    @test size(Kval) == (2, 2)
    @test eltype(Kval) <: Complex
    @test isapprox(Kval, zeros(ComplexF64, 2, 2); atol = 1e-14)
end

@testset "T-GUARDRAIL: untwisted segment has zero twist rate" begin
    f = _straight_fiber(_XS; length = 1.0)
    @test twist_rate(f.path, 0.5) == 0.0
    @test twist_phase(f.path, 0.5) == 0.0
end

@testset "T-GUARDRAIL: no :T_K ⇒ local_temperature ≡ T_ref_K" begin
    f = _straight_fiber(_XS; length = 1.0, T_ref_K = 300.0)
    @test local_temperature(f, 0.0) == 300.0
    @test local_temperature(f, 0.5) == 300.0
    @test local_temperature(f, 1.0) == 300.0
end

@testset "T-GUARDRAIL: generator stays 2×2 ComplexF64 with twist + ellipticity" begin
    xs = _elliptical(; ratio = 1.05, angle = π / 6)
    f = _straight_fiber(xs; length = 1.0, twist = 3.0, spin_rate = 1.0)
    for (mk, label) in ((generator_K(f, _λ), "K"), (generator_Kω(f, _λ), "Kω"))
        v = mk(0.4)
        @test size(v) == (2, 2)
        @test eltype(v) <: Complex
        # Linear + circular generators are traceless (lossless SU(2)).
        @test isapprox(v[1, 1] + v[2, 2], 0; atol = 1e-12)
    end
end
