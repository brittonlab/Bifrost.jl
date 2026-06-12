using Bifrost
using Test
using LinearAlgebra

# Generator-level verification of the fiber birefringence assembly in
# `fiber-path.jl` (issue #11): the assembled `generator_Kω` is checked against a
# finite difference of `generator_K` over ω, the generator axis orientations are
# checked against the path phase accumulators, and the realized propagation is
# checked against closed-form physics.

const FP = Bifrost.FiberPath

const _GEN_T = 297.15
const _GEN_C = 299_792_458.0
const _GEN_λ = 1550e-9

const _GEN_XS_CIRC = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0), 8.2e-6, 125e-6)
const _GEN_XS_ELL = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0), 8.2e-6, 125e-6;
    ellipticity_axis_ratio = 1.01, ellipticity_axis_angle = π / 7)

@testset "generator_Kω matches finite difference of generator_K over ω" begin
    # T-GUARDRAIL: the assembled Kω closure must be the ω-derivative of the
    # assembled K closure. The path exercises every additive generator term:
    # a twisted straight (twist + intrinsic ellipticity), a twisted tensioned
    # planar bend (bend + tension + twist), and a helix (torsion-rotated bend
    # axis), all under continuous spin with an elliptical core. A central
    # difference with relative step 1e-7 in ω agrees to ~4e-7, so rtol = 1e-4
    # passes cleanly while still failing on a single-term coefficient error of
    # a few percent (the issue #11 chain-rule bug produced 4.6%).
    sb = SubpathBuilder(); start!(sb; spin_rate = 2π)
    straight!(sb; length = 0.3, twist = 1.5)
    bend!(sb; radius = 0.04, angle = π / 2, twist = 0.8, meta = [MCMadd(:tension, 0.5)])
    helix!(sb; radius = 0.03, pitch = 0.01, turns = 1.5)
    seal!(sb)
    f = Fiber(Subpath(sb); cross_section = _GEN_XS_ELL, T_ref_K = _GEN_T)

    ω0 = 2π * _GEN_C / _GEN_λ
    dω = ω0 * 1e-7
    K_plus  = generator_K(f, 2π * _GEN_C / (ω0 + dω))
    K_minus = generator_K(f, 2π * _GEN_C / (ω0 - dω))
    Kω = generator_Kω(f, _GEN_λ)

    bps = sort(unique(Float64.(fiber_breakpoints(f))))
    @test length(bps) >= 4   # straight | bend | helix boundaries present
    for (a, b) in zip(bps[begin:end-1], bps[begin+1:end])
        b - a < 1e-9 && continue
        s = (a + b) / 2
        fd = (K_plus(s) .- K_minus(s)) ./ (2 * dω)
        @test norm(fd) > 0
        @test Kω(s) ≈ fd rtol = 1e-4
    end
end

@testset "bend axis on a helix is oriented by the torsion phase" begin
    # T-PHYSICS: in the parallel-transport (Bishop) frame the curvature
    # direction rotates at the geometric-torsion rate τ_g, so the bend generator
    # on a helix equals an independently constructed linear birefringence
    # generator oriented at φ = ∫₀ˢ τ_g = torsion_phase(s). For a helix of
    # radius R and pitch p (h = p/2π), τ_g = h/(R² + h²) is constant.
    sb = SubpathBuilder(); start!(sb)
    helix!(sb; radius = 0.03, pitch = 0.01, turns = 1.5)
    seal!(sb)
    f = Fiber(Subpath(sb); cross_section = _GEN_XS_CIRC, T_ref_K = _GEN_T)

    s = 0.6 * Float64(arc_length(f.path))
    φ = torsion_phase(f.path, s)
    h = 0.01 / (2π)
    @test φ ≈ (h / (0.03^2 + h^2)) * s rtol = 1e-9
    @test abs(φ) > 0.1   # orientation test is non-trivial

    κ = curvature(f.path, s)
    Δβ = bending_birefringence(_GEN_XS_CIRC, _GEN_λ, _GEN_T; bend_radius_m = inv(κ))
    expected = FP.linear_birefringence_generator(Δβ, cos(2φ), sin(2φ))
    @test FP.bend_generator_K(f, s, _GEN_λ) ≈ expected rtol = 1e-12
end

@testset "intrinsic axes follow ellipse angle + spin phase + twist phase" begin
    # T-PHYSICS: the frozen ellipse axes co-rotate with the glass, so the
    # ellipticity generator on a spun, twisted straight fiber is oriented at
    # φ = ellipticity_axis_angle + spin_phase(s) + twist_phase(s).
    sb = SubpathBuilder(); start!(sb; spin_rate = 2π)
    straight!(sb; length = 0.5, twist = 1.5)
    jumpto!(sb; point = (0.0, 0.0, 0.5), incoming_tangent = (0.0, 0.0, 1.0))
    f = Fiber(Subpath(sb); cross_section = _GEN_XS_ELL, T_ref_K = _GEN_T)

    s = 0.2
    φ = _GEN_XS_ELL.ellipticity_axis_angle + spin_phase(f.path, s) + twist_phase(f.path, s)
    @test φ ≈ π / 7 + 2π * s + 1.5 * s rtol = 1e-9
    Δβ = core_noncircularity_birefringence(_GEN_XS_ELL, _GEN_λ, _GEN_T) +
         asymmetric_thermal_stress_birefringence(_GEN_XS_ELL, _GEN_λ, _GEN_T)
    expected = FP.linear_birefringence_generator(Δβ, cos(2φ), sin(2φ))
    @test FP.ellipticity_generator_K(f, s, _GEN_λ) ≈ expected rtol = 1e-12
end

@testset "planar bend propagates as a linear retarder with Δβ_b·L retardance" begin
    # T-PHYSICS: a planar bend has zero torsion phase, so K is diagonal,
    # K = diag(iΔβ_b/2, −iΔβ_b/2), and J(L) = diag(e^{iΔβ_b L/2}, e^{−iΔβ_b L/2})
    # — a pure linear retarder of retardance Δβ_b·L on the curvature axis.
    R = 0.05
    θ = π / 2
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = R, angle = θ)
    seal!(sb)
    f = Fiber(Subpath(sb); cross_section = _GEN_XS_CIRC, T_ref_K = _GEN_T)

    L = R * θ
    Δβ = bending_birefringence(_GEN_XS_CIRC, _GEN_λ, _GEN_T; bend_radius_m = R)
    J, _ = propagate_fiber(f; λ_m = _GEN_λ, verbose = false)
    @test J[1, 1] ≈ exp(0.5im * Δβ * L) atol = 1e-7
    @test J[2, 2] ≈ exp(-0.5im * Δβ * L) atol = 1e-7
    @test abs(J[1, 2]) < 1e-10
    @test abs(J[2, 1]) < 1e-10
end

@testset "twist rotation pinned at θ/(τL) = (1 + n²(p₁₁−p₁₂)/2)/2" begin
    # T-SIM-REGRESSION: a straight fiber twisted at rate τ rotates linear
    # polarization by θ = Δβc·L/2 with Δβc = (1 + n²(p₁₁−p₁₂)/2)·τ, i.e.
    # θ/(τL) ≈ 0.4216 for this cross section at 1550 nm (matches legacy
    # `fibers.py`). Whether the physical rotation is this value or twice it is
    # unresolved (issue #11; see the FLAG above `twisting_dω` and
    # doi:10.1016/j.yofte.2011.10.001) — update this pin intentionally if the
    # twist model is adjudicated to the doubled value.
    τ = 2.0
    Lf = 1.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = Lf, twist = τ)
    jumpto!(sb; point = (0.0, 0.0, Lf), incoming_tangent = (0.0, 0.0, 1.0))
    f = Fiber(Subpath(sb); cross_section = _GEN_XS_CIRC, T_ref_K = _GEN_T)

    J, _ = propagate_fiber(f; λ_m = _GEN_λ, verbose = false)
    # J is a real SO(2) rotation; extract the angle from the first column.
    @test abs(imag(J[1, 1])) < 1e-10
    @test abs(imag(J[2, 1])) < 1e-10
    θ = atan(real(J[2, 1]), real(J[1, 1]))
    @test θ / (τ * Lf) ≈ 0.4215565707 atol = 1e-8
end
