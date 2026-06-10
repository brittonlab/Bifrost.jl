using Test
using LinearAlgebra
using Bifrost
using Bifrost.PathGeometry: _qc_nominalize
using MonteCarloMeasurements

# Fiber-layer per-segment tension (`:tension`). Like `:T_K`, an
# absolute axial tension in Newtons plays a dual role: it elongates the segment
# by (1 + ε) with ε = F/(π·r_clad²·E) using the cladding stiffness, and it sets
# the segment's axial-tension photoelastic birefringence (a linear birefringence
# on the bend eigen-axis, ∝ 1/R, recovered on demand via `tension(f, s)`).

const _FN_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),     # pure silica cladding
    8.2e-6,
    125e-6,
)
const _FN_T_REF = 297.15
const _FN_LAMBDA = 1550e-9
const _FN_R_CLAD = cladding_radius(_FN_XS)
const _FN_E_CLAD = youngs_modulus(_FN_XS.cladding_material, _FN_T_REF)

# Axial strain ε = F / (π·r_clad²·E), computed independently of the source.
_fn_strain(F) = F / (π * _FN_R_CLAD^2 * _FN_E_CLAD)

function _fn_subpath(f::Function; spin_rate = nothing)
    sb = SubpathBuilder(); start!(sb; spin_rate = spin_rate)
    f(sb)
    isnothing(sb.jumpto_point) && seal!(sb)
    return Subpath(sb)
end

_fn_fiber(f; spin_rate = nothing) =
    Fiber(_fn_subpath(f; spin_rate = spin_rate);
          cross_section = _FN_XS, T_ref_K = _FN_T_REF)

# Arc length (nominal) of placed segment `i` of a fiber's path.
_fn_seg_len(f, i) =
    Float64(_qc_nominalize(arc_length(f.path.placed_segments[i].segment)))

# Midpoint arc length (nominal) of placed segment `i`.
function _fn_seg_mid(f, i)
    ps = f.path.placed_segments[i]
    return Float64(_qc_nominalize(ps.s_offset_eff)) + 0.5 * _fn_seg_len(f, i)
end

# -----------------------------------------------------------------------
# Birefringence generator (T-PHYSICS)
# -----------------------------------------------------------------------

@testset "Fiber :tension — straight segment ⇒ zero generator (∝ 1/R)" begin
    # T-PHYSICS: axial-tension birefringence scales as 1/R, so a tensioned but
    # unbent segment contributes nothing to K.
    f = _fn_fiber() do sb
        straight!(sb; length = 1.0, meta = [MCMadd(:tension, 0.5)])
    end
    s = _fn_seg_mid(f, 1)
    @test Bifrost.FiberPath.tension_generator_K(f, s, _FN_LAMBDA) == zeros(ComplexF64, 2, 2)
end

@testset "Fiber :tension — bent segment matches independent birefringence call" begin
    # T-PHYSICS: on a (planar) bend the tension generator equals a linear
    # birefringence generator built from an independent `axial_tension_birefringence`
    # call, on the bend eigen-axis (planar ⇒ c2φ = 1, s2φ = 0).
    F = 2.0
    f = _fn_fiber() do sb
        straight!(sb; length = 0.5)
        bend!(sb; radius = 0.05, angle = π / 2, axis_angle = 0.0,
              meta = [MCMadd(:tension, F)])
    end
    s = _fn_seg_mid(f, 2)
    κ = curvature(f.path, s)
    Δβ = axial_tension_birefringence(_FN_XS, _FN_LAMBDA, temperature(f, s);
                                     bend_radius_m = inv(κ), axial_tension_N = tension(f, s))
    expected = linear_birefringence_generator(Δβ, 1.0, 0.0)
    G = Bifrost.FiberPath.tension_generator_K(f, s, _FN_LAMBDA)
    @test G ≈ expected atol = 1e-18
    @test G[1, 1] != 0                       # non-trivial response
end

@testset "Fiber :tension — generator is the 4th additive term of K" begin
    # T-GUARDRAIL: the full step-index generator is the sum of its four terms, and
    # the tension term is a non-trivial contributor on a tensioned bend.
    f = _fn_fiber() do sb
        straight!(sb; length = 0.5)
        bend!(sb; radius = 0.05, angle = π / 2, axis_angle = 0.0,
              meta = [MCMadd(:tension, 2.0)])
    end
    s = _fn_seg_mid(f, 2)
    FP = Bifrost.FiberPath
    K = generator_K(f, _FN_LAMBDA)(s)
    K_sum = FP.bend_generator_K(f, s, _FN_LAMBDA) +
            FP.twist_generator_K(f, s, _FN_LAMBDA) +
            FP.ellipticity_generator_K(f, s, _FN_LAMBDA) +
            FP.tension_generator_K(f, s, _FN_LAMBDA)
    @test K ≈ K_sum atol = 1e-18
    @test FP.tension_generator_K(f, s, _FN_LAMBDA) != zeros(ComplexF64, 2, 2)
end

# -----------------------------------------------------------------------
# Tension-induced length change (T-PHYSICS)
# -----------------------------------------------------------------------

@testset "Fiber :tension — elongates a segment by (1 + ε)" begin
    # T-PHYSICS: a bend under tension F has placed arc length scaled by
    # (1 + F/(π·r_clad²·E)) vs the untensioned build; ε from the analytic formula.
    F = 3.0
    spec = sb -> begin
        straight!(sb; length = 0.5)
        bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:tension, F)])
    end
    spec0 = sb -> begin
        straight!(sb; length = 0.5)
        bend!(sb; radius = 0.05, angle = π / 2)
    end
    f  = _fn_fiber(spec)
    f0 = _fn_fiber(spec0)
    @test _fn_seg_len(f, 2) ≈ (1 + _fn_strain(F)) * _fn_seg_len(f0, 2) rtol = 1e-12
    # Lead-in (no :tension) is unchanged.
    @test _fn_seg_len(f, 1) ≈ _fn_seg_len(f0, 1) rtol = 1e-12
end

@testset "Fiber :tension — divides twist rate by (1 + ε); conserves turns" begin
    # T-PHYSICS: tension elongation scales arc length by (1 + ε), so mechanical
    # twist (an inverse-length rate) divides by (1 + ε) and total turns ∫τ_m ds are
    # conserved — the same inverse-length scaling as :T_K.
    F   = 3.0
    τm0 = 5.0
    f  = _fn_fiber(sb -> straight!(sb; length = 1.0, twist = τm0,
                                   meta = [MCMadd(:tension, F)]))
    f0 = _fn_fiber(sb -> straight!(sb; length = 1.0, twist = τm0))
    L  = Float64(_qc_nominalize(arc_length(f.path)))
    L0 = Float64(_qc_nominalize(arc_length(f0.path)))

    @test twist_rate(f.path, 0.5 * L) ≈ τm0 / (1 + _fn_strain(F)) rtol = 1e-12
    # Total twist conserved: Φ unchanged from the untensioned build.
    @test twist_phase(f.path, L) ≈ twist_phase(f0.path, L0) rtol = 1e-9
end

@testset "Fiber :tension — composes multiplicatively with :T_K" begin
    # T-PHYSICS: a segment carrying both scales by (1 + α_lin·ΔT)·(1 + ε).
    F  = 2.0
    ΔT = 10.0
    α  = cte(_FN_XS.cladding_material, _FN_T_REF)
    spec  = sb -> straight!(sb; length = 1.0, meta = [MCMadd(:T_K, ΔT), MCMadd(:tension, F)])
    spec0 = sb -> straight!(sb; length = 1.0)
    f  = _fn_fiber(spec)
    f0 = _fn_fiber(spec0)
    expected = (1 + α * ΔT) * (1 + _fn_strain(F))
    @test _fn_seg_len(f, 1) ≈ expected * _fn_seg_len(f0, 1) rtol = 1e-12
end

# -----------------------------------------------------------------------
# tension(f, s) query (T-GUARDRAIL)
# -----------------------------------------------------------------------

@testset "Fiber :tension — tension(f, s) is 0 everywhere without :tension" begin
    # T-GUARDRAIL
    f = _fn_fiber() do sb
        straight!(sb; length = 1.0)
        bend!(sb; radius = 0.05, angle = π / 2)
    end
    @test tension(f, _fn_seg_mid(f, 1)) == 0.0
    @test tension(f, _fn_seg_mid(f, 2)) == 0.0
end

@testset "Fiber :tension — tension(f, s) recovers the authored value per segment" begin
    # T-GUARDRAIL: derived on demand from segment meta; 0 on the untensioned lead-in.
    F = 0.75
    f = _fn_fiber() do sb
        straight!(sb; length = 1.0)
        bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:tension, F)])
    end
    @test tension(f, _fn_seg_mid(f, 1)) == 0.0
    @test tension(f, _fn_seg_mid(f, 2)) ≈ F
end

@testset "Fiber :tension — negative tension rejected through the generator" begin
    # T-GUARDRAIL: the cross-section nonnegative `axial_tension_N` guard fires when
    # the tension generator is evaluated on a bent, negatively-tensioned segment.
    f = _fn_fiber() do sb
        straight!(sb; length = 0.5)
        bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:tension, -1.0)])
    end
    s = _fn_seg_mid(f, 2)
    @test_throws ArgumentError Bifrost.FiberPath.tension_generator_K(f, s, _FN_LAMBDA)
end

# -----------------------------------------------------------------------
# Terminal-connector seal meta (T-GUARDRAIL)
# -----------------------------------------------------------------------

@testset "Fiber :tension — jumpto! seal accepts :tension; elongates connector" begin
    # T-GUARDRAIL: a :tension on the jumpto! seal scales the terminal
    # connector's arc length by (1 + ε), still landing at the fixed jumpto_point.
    F = 5.0
    P = (0.1, 0.0, 0.5)
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.5)
    jumpto!(sb; point = P, incoming_tangent = (1.0, 0.0, 0.0),
            meta = [MCMadd(:tension, F)])
    L0 = Float64(_qc_nominalize(arc_length(build(sb).jumpto_quintic_connector)))

    f = Fiber(sb; cross_section = _FN_XS, T_ref_K = _FN_T_REF)
    @test isapprox(
        Float64(_qc_nominalize(arc_length(f.path.jumpto_quintic_connector))),
        (1 + _fn_strain(F)) * L0; rtol = 1e-6)
    s_end = Float64(_qc_nominalize(arc_length(f.path)))
    @test isapprox(collect(position(f.path, s_end)), collect(P); atol = 1e-6)
end

@testset "Fiber :tension — unsupported seal meta still errors" begin
    # T-GUARDRAIL: the extended allow-list permits :T_K and :tension only; a
    # field-level MCMadd or any MCMmul on the seal is still rejected.
    mk(meta) = begin
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 0.5)
        jumpto!(sb; point = (0.1, 0.0, 0.5),
                incoming_tangent = (1.0, 0.0, 0.0), meta = meta)
        sb
    end
    @test Fiber(mk([MCMadd(:tension, 1.0)]);
                cross_section = _FN_XS, T_ref_K = _FN_T_REF) isa Fiber
    @test_throws ArgumentError Fiber(mk([MCMadd(:length, 0.01)]);
                                     cross_section = _FN_XS, T_ref_K = _FN_T_REF)
    @test_throws ArgumentError Fiber(mk([MCMmul(:tension, 1.1)]);
                                     cross_section = _FN_XS, T_ref_K = _FN_T_REF)
end

# -----------------------------------------------------------------------
# MCM: Particles flow through query, generator, and length scaling
# -----------------------------------------------------------------------

@testset "Fiber :tension — Particles propagate through query/generator/length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        F = 2.0 ± 0.1
        f = _fn_fiber() do sb
            straight!(sb; length = 0.5)
            bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:tension, F)])
        end
        s = _fn_seg_mid(f, 2)
        # Length scaling lifts the bend segment to Particles.
        @test arc_length(f.path.placed_segments[2].segment) isa Particles
        # The query carries Particles through.
        @test tension(f, s) isa Particles
        @test pmean(tension(f, s)) ≈ 2.0 rtol = 1e-6
        # The generator entries carry Particles.
        G = Bifrost.FiberPath.tension_generator_K(f, s, _FN_LAMBDA)
        @test eltype(G) <: Complex{<:Particles}
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end
