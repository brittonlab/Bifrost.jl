using Test
using LinearAlgebra
using Bifrost

# -----------------------------------------------------------------------
# Test setup: a small helper to build a Fiber with a pure-silica cladding
# (so α_lin = SILICA_CTE) and a reference temperature.
# -----------------------------------------------------------------------

const _MODIFY_TEST_XS = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),     # pure silica cladding → α_lin = SILICA_CTE
    8.2e-6,
    125e-6,
)

const _MODIFY_TEST_T_REF = 297.15
const _MODIFY_ALPHA_LIN  = cte(_MODIFY_TEST_XS.cladding_material, _MODIFY_TEST_T_REF)

# Helper: ΔT that produces a given total-length α = 1 + α_lin·ΔT.
_ΔT_for(α) = (α - 1) / _MODIFY_ALPHA_LIN
_mcm(α)    = [MCMadd(:T_K, _ΔT_for(α))]

_modify_fiber_path(path) = modify(Fiber(path;
    cross_section = _MODIFY_TEST_XS, T_ref_K = _MODIFY_TEST_T_REF))

# -----------------------------------------------------------------------
# Uniform thermal shrinkage via :T_K
# -----------------------------------------------------------------------

@testset "modify — :T_K uniform scales path length" begin
    # T-PHYSICS: uniform α on every segment multiplies total arc length by α
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0, meta = _mcm(1.05))
    bend!(spec; radius = 0.1, angle = π / 2, meta = _mcm(1.05))
    path = build(spec)

    path_scaled = _modify_fiber_path(path)
    @test path_length(path_scaled) ≈ path_length(path) * 1.05 atol = 1e-10
end

@testset "modify — :T_K preserves joint tangent continuity" begin
    # T-GUARDRAIL
    spec = PathSpecBuilder()
    straight!(spec; length = 0.5, meta = _mcm(0.98))
    bend!(spec; radius = 0.1, angle = π / 3, meta = _mcm(0.98))
    path = _modify_fiber_path(build(spec))

    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "modify — :T_K on BendSegment preserves angle, scales radius" begin
    # T-PHYSICS: α scales R_eff = α·R; swept angle preserved.
    spec = PathSpecBuilder()
    bend!(spec; radius = 0.1, angle = π / 3, axis_angle = 0.0, meta = _mcm(1.1))
    path1 = build(spec)
    path2 = _modify_fiber_path(path1)

    seg1 = path1.placed_segments[1].segment
    seg2 = path2.placed_segments[1].segment
    @test seg2.radius ≈ 1.1 * seg1.radius
    @test seg2.angle == seg1.angle
    @test arc_length(seg2) ≈ 1.1 * arc_length(seg1)
end

@testset "modify — :T_K on CatenarySegment scales arc length and parameter a" begin
    # T-PHYSICS
    spec = PathSpecBuilder()
    catenary!(spec; a = 0.2, length = 1.0, meta = _mcm(0.95))
    path = build(spec)
    path_s = _modify_fiber_path(path)

    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.a ≈ 0.95 * seg.a
    @test seg_s.length ≈ 0.95 * seg.length
    @test arc_length(seg_s) ≈ 0.95 * arc_length(seg)
end

@testset "modify — :T_K on HelixSegment scales arc length but preserves turns" begin
    spec = PathSpecBuilder()
    helix!(spec; radius = 0.03, pitch = 0.01, turns = 2.0, meta = _mcm(0.9))
    path = build(spec)
    path_s = _modify_fiber_path(path)

    seg = path.placed_segments[1].segment
    seg_s = path_s.placed_segments[1].segment
    @test seg_s.turns == seg.turns
    @test seg_s.radius ≈ 0.9 * seg.radius
    @test seg_s.pitch ≈ 0.9 * seg.pitch
    @test arc_length(seg_s) ≈ 0.9 * arc_length(seg)
end

@testset "modify — :T_K on JumpBy endpoint scales" begin
    # T-PHYSICS: uniform α scales all positions
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.0, 0.0, 1.0), meta = _mcm(0.8))
    path = build(spec)
    path_s = _modify_fiber_path(path)
    @test end_point(path_s) ≈ [0.0, 0.0, 0.8] atol = 1e-10
end

# -----------------------------------------------------------------------
# Per-segment independence (meta on some segments, not others)
# -----------------------------------------------------------------------

@testset "modify — per-segment scaling applies independently" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0, meta = _mcm(1.0))   # α = 1 via ΔT = 0
    straight!(spec; length = 1.0, meta = _mcm(0.5))
    path = build(spec)

    path_s = _modify_fiber_path(path)
    @test path_length(path_s) ≈ 1.5 atol = 1e-10
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 0.5
end

@testset "modify — segments without MCM annotations default to α = 1.0" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)                     # no meta → α = 1
    straight!(spec; length = 2.0, meta = _mcm(0.5))
    path = build(spec)

    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path_s.placed_segments[2].segment) ≈ 1.0
end

@testset "modify — non-matching MCM symbols are ignored" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0,
              meta = [MCMadd(:not_a_field, 1e6),
                                  Nickname("labelled")])
    path = build(spec)

    path_s = _modify_fiber_path(path)
    @test arc_length(path_s.placed_segments[1].segment) ≈ 1.0
end

# -----------------------------------------------------------------------
# Direct field-level MCMs (new in fiber-path-modify.jl)
# -----------------------------------------------------------------------

@testset "modify — direct field :length on StraightSegment" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0, meta = [MCMadd(:length, 0.05)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.length ≈ 1.05 atol = 1e-12
end

@testset "modify — direct field :radius on BendSegment" begin
    spec = PathSpecBuilder()
    bend!(spec; radius = 0.05, angle = π / 2,
          meta = [MCMadd(:radius, 0.01)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.radius ≈ 0.06 atol = 1e-12
    @test seg.angle ≈ π / 2
end

@testset "modify — direct field :angle on BendSegment" begin
    spec = PathSpecBuilder()
    bend!(spec; radius = 0.05, angle = π / 2,
          meta = [MCMadd(:angle, π / 12)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.angle ≈ π / 2 + π / 12 atol = 1e-12
    @test seg.radius ≈ 0.05
end

@testset "modify — direct field :pitch on HelixSegment (MCMmul)" begin
    spec = PathSpecBuilder()
    helix!(spec; radius = 0.03, pitch = 0.01, turns = 2.0,
           meta = [MCMmul(:pitch, 1.1)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.pitch ≈ 0.01 * 1.1 atol = 1e-14
    @test seg.radius ≈ 0.03
end

@testset "modify — direct field :a on CatenarySegment" begin
    spec = PathSpecBuilder()
    catenary!(spec; a = 0.2, length = 1.0,
              meta = [MCMadd(:a, 0.002)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    @test seg.a ≈ 0.202 atol = 1e-12
    @test seg.length ≈ 1.0
end

@testset "modify — combined :T_K and direct :length on StraightSegment" begin
    # Order: T_K first (multiplicative length scaling via CTE), then the direct
    # :length additive offset on top.
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0,
              meta = [MCMadd(:T_K, 10.0),
                                  MCMadd(:length, 0.001)])
    path = build(spec)
    seg = _modify_fiber_path(path).placed_segments[1].segment
    expected = 1.0 * (1 + _MODIFY_ALPHA_LIN * 10.0) + 0.001
    @test seg.length ≈ expected atol = 1e-14
end

# -----------------------------------------------------------------------
# Twist overlay remapping
# -----------------------------------------------------------------------

# TODO: twist refactor — pending per-segment-meta twist subsystem.
@testset "modify — :T_K preserves constant twist rate; total scales with length" begin
    @test_skip true
end

# -----------------------------------------------------------------------
# MCMadd: Particles flow through modify via a Particles-valued ΔT
# -----------------------------------------------------------------------

if !isdefined(Main, :Particles)
    using MonteCarloMeasurements
end

@testset "modify — MCMadd Particles ΔT scales arc length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # Particles-valued ΔT with zero mean → nominal arc length preserved,
        # output lifts to Particles.
        ΔT = 0.0 ± (0.01 / _MODIFY_ALPHA_LIN)    # σ_α = 0.01 around α = 1
        spec = PathSpecBuilder()
        bend!(spec; radius = 0.05, angle = π / 2,
              meta = [MCMadd(:T_K, ΔT)])
        path = build(spec)
        path_s = _modify_fiber_path(path)

        @test arc_length(path_s) isa Particles
        @test pmean(arc_length(path_s)) ≈ 0.05 * (π / 2) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end


# -----------------------------------------------------------------------
# QuinticConnector :T_K scaling
# -----------------------------------------------------------------------

@testset "modify — :T_K scales QuinticConnector arc length" begin
    """
    A JumpBy with :T_K MCM meta must scale through to the resolved
    QuinticConnector after build(). The geometric scale τ should multiply the
    coefficient matrix and the arc-length table uniformly so the post-modify
    path length scales by τ.
    """
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.4, 0.0, 0.4),
            tangent = (1.0, 0.0, 0.0),
            meta = _mcm(0.5))
    path = build(spec)
    L_before = path_length(path)
    path_s = _modify_fiber_path(path)
    L_after = path_length(path_s)
    @test isapprox(L_after, 0.5 * L_before; rtol = 1e-10, atol = 1e-12)
end

@testset "modify — :T_K on QuinticConnector preserves end-point shrinkage" begin
    """
    Geometric scaling of the connector by τ must move its endpoint to τ·delta
    in the local frame. Verifies that QuinticConnector :T_K modify uses the
    same coefficient-and-table scaling pattern as the prior HermiteConnector.
    """
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.0, 0.0, 1.0), meta = _mcm(0.5))
    path = build(spec)
    path_s = _modify_fiber_path(path)
    @test end_point(path_s) ≈ [0.0, 0.0, 0.5] atol = 1e-10
end

@testset "modify — T-GUARDRAIL: upstream bend changes recompute connector K0" begin
    spec = PathSpecBuilder()
    bend!(spec; radius = 1.0, angle = π / 3,
          meta = [MCMmul(:radius, 2.0)])
    jumpby!(spec; delta = (1.0, 0.0, 0.2))

    path_s = _modify_fiber_path(build(spec))
    bend_seg = path_s.placed_segments[1].segment
    connector = path_s.placed_segments[2].segment

    @test curvature(bend_seg, arc_length(bend_seg)) ≈ 0.5 atol = 1e-12
    @test curvature(connector, 0.0) ≈ 0.5 atol = 1e-12
end

@testset "modify — T-GUARDRAIL: Twist anchors tolerate MCM-valued modified length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        spec = PathSpecBuilder()
        straight!(spec; length = 1.0,
                  meta = [
                      Twist(; rate = 2.0),
                      MCMadd(:length, 0.0 ± 0.01),
                  ])

        path_s = _modify_fiber_path(build(spec))

        @test arc_length(path_s) isa Particles
        @test total_material_twist(path_s) ≈ 2.0 * pmean(arc_length(path_s)) rtol = 1e-6
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# -----------------------------------------------------------------------
# JumpTo + :T_K — length-constrained connector resolve
# -----------------------------------------------------------------------
# JumpTo destination is a lab-frame invariant. Thermal expansion on the
# connector itself must show up as arc length, not as a chord rescaling
# (which would move the endpoint). modify() folds this into the resolver
# via target_arc_length = τ_conn · baseline_L.

@testset "modify — :T_K on JumpTo preserves destination" begin
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (0.4, 0.0, 0.4),
            tangent = (1.0, 0.0, 0.0),
            meta = _mcm(1.05))
    path = build(spec)
    path_s = _modify_fiber_path(path)
    @test end_point(path_s) ≈ [0.4, 0.0, 0.4] atol = 1e-6
end

@testset "modify — :T_K on JumpTo scales arc length only" begin
    α = 1.05
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (0.4, 0.0, 0.4),
            tangent = (1.0, 0.0, 0.0),
            meta = _mcm(α))
    path = build(spec)
    L_baseline = arc_length(path.placed_segments[1].segment)
    path_s = _modify_fiber_path(path)
    @test end_point(path_s) ≈ [0.4, 0.0, 0.4] atol = 1e-6
    L_modified = arc_length(path_s.placed_segments[1].segment)
    @test L_modified ≈ α * L_baseline rtol = 1e-4
end

@testset "modify — JumpTo without :T_K is unchanged" begin
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (0.5, 0.0, 0.5),
            tangent = (1.0, 0.0, 0.0))
    path = build(spec)
    path_s = _modify_fiber_path(path)
    @test end_point(path_s) ≈ [0.5, 0.0, 0.5] atol = 1e-10
    @test path_length(path_s) ≈ path_length(path) atol = 1e-10
end

@testset "modify — JumpTo absorbs upstream thermal drift" begin
    α = 1.05
    spec = PathSpecBuilder()
    straight!(spec; length = 0.5, meta = _mcm(α))
    jumpto!(spec; destination = (0.0, 0.0, 1.0))
    path = build(spec)
    path_s = _modify_fiber_path(path)

    # Endpoint pinned at lab-frame destination despite upstream expansion.
    @test end_point(path_s) ≈ [0.0, 0.0, 1.0] atol = 1e-8

    # Straight expanded by α; connector has no own :T_K so no length target,
    # arc length is whatever the unconstrained smooth ramp produces. Sanity:
    # connector arc length ≈ |1 - α·0.5| (chord magnitude, since tangents
    # default to chord-aligned for a straight-on jump).
    seg_straight = path_s.placed_segments[1].segment
    @test seg_straight.length ≈ α * 0.5 atol = 1e-12
    L_conn = arc_length(path_s.placed_segments[2].segment)
    @test L_conn ≈ (1.0 - α * 0.5) rtol = 1e-3
end

@testset "modify — JumpTo with both upstream and own :T_K" begin
    α_up = 1.05    # upstream straight expansion
    α_jt = 1.10    # connector's own thermal expansion
    spec = PathSpecBuilder()
    straight!(spec; length = 0.5, meta = _mcm(α_up))
    jumpto!(spec; destination = (0.0, 0.0, 1.0), meta = _mcm(α_jt))
    path = build(spec)
    L_conn_baseline = arc_length(path.placed_segments[2].segment)
    path_s = _modify_fiber_path(path)

    @test end_point(path_s) ≈ [0.0, 0.0, 1.0] atol = 1e-6
    seg_straight = path_s.placed_segments[1].segment
    @test seg_straight.length ≈ α_up * 0.5 atol = 1e-12
    L_conn_modified = arc_length(path_s.placed_segments[2].segment)
    @test L_conn_modified ≈ α_jt * L_conn_baseline rtol = 1e-4
end

@testset "modify — JumpTo with infeasible :T_K shrinkage throws" begin
    # τ < 1 with a chord that already equals baseline arc length: target
    # falls below chord magnitude, which is geometrically infeasible.
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (0.0, 0.0, 1.0),
            tangent = (0.0, 0.0, 1.0),
            meta = _mcm(0.5))
    path = build(spec)
    @test_throws ArgumentError _modify_fiber_path(path)
end
