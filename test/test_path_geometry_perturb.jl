using Test
using LinearAlgebra
using Bifrost
using Bifrost.PathGeometry: _qc_nominalize
using MonteCarloMeasurements

# Geometry-layer perturbation: build(sub; perturb=true) applies field-level
# MCMadd/MCMmul to a segment's own fields, leaves foreign meta (e.g. a thermal
# annotation it cannot interpret) inert, and never errors on it.

# Build a sealed Subpath from a do-block that authors interior segments. Seals
# at the natural exit (via seal!) if no jumpto! was called.
function _subpath(f::Function)
    sb = SubpathBuilder(); start!(sb)
    f(sb)
    isnothing(sb.jumpto_point) && seal!(sb)
    return Subpath(sb)
end

_perturb(f::Function) = build(_subpath(f); perturb = true)

# -----------------------------------------------------------------------
# Field-level MCM on a segment's own fields
# -----------------------------------------------------------------------

@testset "perturb — direct field :length on StraightSegment" begin
    seg = _perturb() do sb
        straight!(sb; length = 1.0, meta = [MCMadd(:length, 0.05)])
    end.placed_segments[1].segment
    @test seg.length ≈ 1.05 atol = 1e-12
end

@testset "perturb — direct field :radius on BendSegment" begin
    seg = _perturb() do sb
        bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:radius, 0.01)])
    end.placed_segments[1].segment
    @test seg.radius ≈ 0.06 atol = 1e-12
    @test seg.angle ≈ π / 2
end

@testset "perturb — direct field :angle on BendSegment" begin
    seg = _perturb() do sb
        bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:angle, π / 12)])
    end.placed_segments[1].segment
    @test seg.angle ≈ π / 2 + π / 12 atol = 1e-12
    @test seg.radius ≈ 0.05
end

@testset "perturb — direct field :pitch on HelixSegment (MCMmul)" begin
    seg = _perturb() do sb
        helix!(sb; radius = 0.03, pitch = 0.01, turns = 2.0, meta = [MCMmul(:pitch, 1.1)])
    end.placed_segments[1].segment
    @test seg.pitch ≈ 0.01 * 1.1 atol = 1e-14
    @test seg.radius ≈ 0.03
end

@testset "perturb — direct field :a on CatenarySegment" begin
    seg = _perturb() do sb
        catenary!(sb; a = 0.2, length = 1.0, meta = [MCMadd(:a, 0.002)])
    end.placed_segments[1].segment
    @test seg.a ≈ 0.202 atol = 1e-12
    @test seg.length ≈ 1.0
end

@testset "perturb — non-matching MCM symbols are ignored" begin
    seg = _perturb() do sb
        straight!(sb; length = 1.0,
                  meta = [MCMadd(:not_a_field, 1e6), Nickname("labelled")])
    end.placed_segments[1].segment
    @test arc_length(seg) ≈ 1.0
end

# -----------------------------------------------------------------------
# Geometry is agnostic of foreign meta (e.g. thermal :T_K)
# -----------------------------------------------------------------------

@testset "perturb — foreign :T_K meta is carried inertly (no error, no change)" begin
    # T-GUARDRAIL: :T_K is not a geometry concern. build(perturb=true) must not
    # error on it and must not change geometry (no α_lin exists here). Only the
    # fiber layer interprets :T_K.
    seg = _perturb() do sb
        straight!(sb; length = 2.0, meta = [MCMadd(:T_K, 100.0)])
    end.placed_segments[1].segment
    @test seg.length == 2.0
    @test arc_length(seg) ≈ 2.0
end

# -----------------------------------------------------------------------
# Connector / spinning interactions under field-MCM
# -----------------------------------------------------------------------

@testset "perturb — upstream bend change recomputes connector K0" begin
    # T-GUARDRAIL: doubling the bend radius via MCMmul halves its curvature; the
    # terminal connector's incoming curvature must track it.
    path = _perturb() do sb
        bend!(sb; radius = 1.0, angle = π / 3, meta = [MCMmul(:radius, 2.0)])
        jumpby!(sb; delta = (1.0, 0.0, 0.2))
    end
    bend_seg  = path.placed_segments[1].segment
    connector = path.placed_segments[2].segment
    @test curvature(bend_seg, arc_length(bend_seg)) ≈ 0.5 atol = 1e-12
    @test curvature(connector, 0.0) ≈ 0.5 atol = 1e-12
end

@testset "perturb — Spinning anchors tolerate MCM-valued modified length" begin
    # T-GUARDRAIL
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        path = _perturb() do sb
            straight!(sb; length = 1.0,
                      meta = [Spinning(; rate = 2.0), MCMadd(:length, 0.0 ± 0.01)])
        end
        seg = path.placed_segments[1].segment
        L_seg = arc_length(seg)
        @test L_seg isa Particles
        # Integrated material spinning over the segment = 2.0 * L_seg.
        @test total_spinning(path; s_start = 0.0,
                             s_end = Float64(_qc_nominalize(L_seg))) ≈
              2.0 * pmean(L_seg) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end
