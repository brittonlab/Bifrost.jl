using Test
using LinearAlgebra
using Bifrost

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

is_unit(v)       = abs(norm(v) - 1.0) < 1e-12
is_orthonormal(T, N, B) = is_unit(T) && is_unit(N) && is_unit(B) &&
                           abs(dot(T, N)) < 1e-12 &&
                           abs(dot(T, B)) < 1e-12 &&
                           abs(dot(N, B)) < 1e-12 &&
                           norm(cross(T, N) - B) < 1e-12   # right-handed

# -----------------------------------------------------------------------
# StraightSegment
# -----------------------------------------------------------------------

@testset "StraightSegment — arc length" begin
    # T-PHYSICS: arc length of a straight segment equals its length
    seg = StraightSegment(3.0)
    @test arc_length(seg) ≈ 3.0
end

@testset "StraightSegment — geometry" begin
    # T-PHYSICS: straight segment has zero curvature and zero torsion
    seg = StraightSegment(5.0)
    @test curvature(seg, 0.0)          == 0.0
    @test curvature(seg, 2.5)          == 0.0
    @test geometric_torsion(seg, 0.0)  == 0.0

    # tangent is always along local z
    @test tangent_local(seg, 0.0) ≈ [0.0, 0.0, 1.0]
    @test tangent_local(seg, 2.5) ≈ [0.0, 0.0, 1.0]

    # position advances along z
    @test position_local(seg, 0.0) ≈ [0.0, 0.0, 0.0]
    @test position_local(seg, 3.0) ≈ [0.0, 0.0, 3.0]

    # frame is orthonormal at any s
    T = tangent_local(seg, 1.0)
    N = normal_local(seg, 1.0)
    B = binormal_local(seg, 1.0)
    @test is_orthonormal(T, N, B)
end

# -----------------------------------------------------------------------
# BendSegment
# -----------------------------------------------------------------------

@testset "BendSegment — arc length and curvature" begin
    # T-PHYSICS: arc_length = radius * |angle|, κ = 1 / radius
    seg = BendSegment(0.1, π / 2)
    @test arc_length(seg) ≈ 0.1 * π / 2
    @test curvature(seg, 0.0) ≈ 1.0 / 0.1
end

@testset "BendSegment — zero torsion (planar curve)" begin
    # T-PHYSICS: circular arc is planar, geometric torsion = 0
    seg = BendSegment(0.05, π)
    @test geometric_torsion(seg, 0.0) == 0.0
    @test geometric_torsion(seg, arc_length(seg) / 2) == 0.0
end

@testset "BendSegment — initial tangent along local z" begin
    # T-GUARDRAIL: all segments must start with tangent (0,0,1) in local coords
    seg = BendSegment(0.2, π / 3, π / 4)
    @test tangent_local(seg, 0.0) ≈ [0.0, 0.0, 1.0]
end

@testset "BendSegment — end position for quarter circle (axis_angle = 0)" begin
    # T-PHYSICS: quarter circle of radius R in x-z plane (axis_angle=0).
    # Start: (0,0,0), tangent z. End: (R, 0, R) tangent x.
    R = 0.1
    seg = BendSegment(R, π / 2, 0.0)
    pos_end = end_position_local(seg)
    @test pos_end ≈ [R, 0.0, R] atol = 1e-12

    (T_end, N_end, B_end) = end_frame_local(seg)
    @test T_end ≈ [1.0, 0.0, 0.0] atol = 1e-12   # tangent rotated 90° toward +x
    @test is_orthonormal(T_end, N_end, B_end)
end

@testset "BendSegment — end position for half circle" begin
    # T-PHYSICS: half circle of radius R in x-z plane.
    # Start: (0,0,0) tangent z. End: (2R, 0, 0) tangent -z.
    R = 0.05
    seg = BendSegment(R, π, 0.0)
    pos_end = end_position_local(seg)
    @test pos_end ≈ [2R, 0.0, 0.0] atol = 1e-12

    (T_end, _, _) = end_frame_local(seg)
    @test T_end ≈ [0.0, 0.0, -1.0] atol = 1e-12   # reversed
end

@testset "BendSegment — axis_angle rotates bend plane" begin
    # T-PHYSICS: axis_angle = π/2 bends in the y-z plane instead of x-z.
    # Quarter circle end position should be (0, R, R).
    R = 0.1
    seg = BendSegment(R, π / 2, π / 2)
    pos_end = end_position_local(seg)
    @test pos_end ≈ [0.0, R, R] atol = 1e-12
end

@testset "BendSegment — frame orthonormality along arc" begin
    # T-GUARDRAIL: frame must remain orthonormal at all points
    seg = BendSegment(0.08, 2π / 3, π / 6)
    for s in range(0.0, arc_length(seg); length = 9)
        T = tangent_local(seg, s)
        N = normal_local(seg, s)
        B = binormal_local(seg, s)
        @test is_orthonormal(T, N, B)
    end
end

# -----------------------------------------------------------------------
# CatenarySegment
# -----------------------------------------------------------------------

@testset "CatenarySegment — initial tangent along local z" begin
    # T-PHYSICS: catenary vertex has vertical tangent (along z in local frame)
    seg = CatenarySegment(0.5, 1.0)
    @test tangent_local(seg, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-12
end

@testset "CatenarySegment — curvature at vertex" begin
    # T-PHYSICS: κ(0) = 1/a (maximum curvature at vertex)
    a = 0.3
    seg = CatenarySegment(a, 1.0)
    @test curvature(seg, 0.0) ≈ 1.0 / a atol = 1e-12
end

@testset "CatenarySegment — curvature formula" begin
    # T-PHYSICS: κ(s) = a / (a² + s²)
    a = 0.4
    seg = CatenarySegment(a, 2.0)
    for s in [0.0, 0.3, 0.7, 1.2]
        @test curvature(seg, s) ≈ a / (a^2 + s^2) atol = 1e-12
    end
end

@testset "CatenarySegment — zero geometric torsion (planar)" begin
    # T-PHYSICS: catenary is a planar curve, τ_geom = 0
    seg = CatenarySegment(0.2, 1.5)
    @test geometric_torsion(seg, 0.0) == 0.0
    @test geometric_torsion(seg, 0.5) == 0.0
end

@testset "CatenarySegment — frame orthonormality along arc" begin
    # T-GUARDRAIL
    seg = CatenarySegment(0.3, 1.0, π / 4)
    for s in range(0.0, arc_length(seg); length = 9)
        T = tangent_local(seg, s)
        N = normal_local(seg, s)
        B = binormal_local(seg, s)
        @test is_orthonormal(T, N, B)
    end
end

# -----------------------------------------------------------------------
# PathSpec assembly and build
# -----------------------------------------------------------------------

@testset "PathSpec — empty spec fails build" begin
    # T-GUARDRAIL: a path with no segments cannot be built
    spec = PathSpecBuilder()
    @test_throws Exception build(spec)
end

@testset "PathSpec — single straight segment" begin
    # T-PHYSICS: a single straight segment of length L produces a path from
    # (0,0,0) to (0,0,L) with constant tangent (0,0,1).
    spec = PathSpecBuilder()
    straight!(spec; length = 2.0)
    path = build(spec)

    @test path_length(path) ≈ 2.0
    @test start_point(path) ≈ [0.0, 0.0, 0.0] atol = 1e-12
    @test end_point(path)   ≈ [0.0, 0.0, 2.0] atol = 1e-12
    @test start_tangent(path) ≈ [0.0, 0.0, 1.0] atol = 1e-12
    @test end_tangent(path)   ≈ [0.0, 0.0, 1.0] atol = 1e-12
end

@testset "Path — tangent continuity at segment joints" begin
    # T-GUARDRAIL: tangent must be continuous across every segment boundary.
    # Test with a straight → bend → straight sequence.
    spec = PathSpecBuilder()
    straight!(spec; length = 0.5)
    bend!(spec; radius = 0.1, angle = π / 2, axis_angle = 0.0)
    straight!(spec; length = 0.3)
    path = build(spec)

    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        # Tangent just before and just after the joint should agree
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "Path — straight + quarter bend geometry" begin
    # T-PHYSICS: straight segment of length L followed by quarter-circle of
    # radius R. End position should be (R, 0, L+R), end tangent should be (1,0,0).
    L = 1.0; R = 0.2
    spec = PathSpecBuilder()
    straight!(spec; length = L)
    bend!(spec; radius = R, angle = π / 2, axis_angle = 0.0)
    path = build(spec)

    @test end_point(path)   ≈ [R, 0.0, L + R] atol = 1e-10
    @test end_tangent(path) ≈ [1.0, 0.0, 0.0]  atol = 1e-10
end

@testset "Path — full circle returns to start" begin
    # T-PHYSICS: a complete circle of radius R returns to the start position
    # with the same tangent direction.
    R = 0.15
    spec = PathSpecBuilder()
    bend!(spec; radius = R, angle = 2π, axis_angle = 0.0)
    path = build(spec)

    @test start_point(path)   ≈ end_point(path)   atol = 1e-10
    @test start_tangent(path) ≈ end_tangent(path)  atol = 1e-10
    @test path_length(path)   ≈ 2π * R             atol = 1e-10
end

@testset "Path — cartesian_distance vs arc_length" begin
    # T-PHYSICS: for a straight segment, cartesian_distance == arc_length.
    # For a curved path, cartesian_distance < arc_length.
    spec_straight = PathSpecBuilder()
    straight!(spec_straight; length = 3.0)
    path_s = build(spec_straight)
    @test cartesian_distance(path_s, 0.0, 3.0) ≈ 3.0 atol = 1e-12

    spec_bent = PathSpecBuilder()
    bend!(spec_bent; radius = 0.5, angle = π)
    path_b = build(spec_bent)
    L = arc_length(path_b)
    @test cartesian_distance(path_b, 0.0, L) < L   # chord < arc
end

@testset "Path — frame orthonormality along assembled path" begin
    # T-GUARDRAIL: frame must remain orthonormal everywhere in the assembled path
    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.08, angle = π / 3, axis_angle = π / 6)
    catenary!(spec; a = 0.2, length = 0.4)
    path = build(spec)

    ss = range(path.spec.s_start, path.s_end; length = 31)
    for s in ss
        T = tangent(path, s)
        N = normal(path, s)
        B = binormal(path, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "Path — bounding box contains all sampled points" begin
    # T-GUARDRAIL
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.2, angle = π / 2)
    path = build(spec)
    bb = bounding_box(path; n = 256)

    for s in range(path.spec.s_start, path.s_end; length = 64)
        p = position(path, s)
        @test all(p .>= bb.lo .- 1e-10)
        @test all(p .<= bb.hi .+ 1e-10)
    end
end

# -----------------------------------------------------------------------
# Twist meta and material_twist
# -----------------------------------------------------------------------

@testset "Twist — constant rate (Float64) is exact" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 2.0, meta = [Twist(; rate = 1.5, phi_0 = 0.0)])
    path = build(spec)
    @test material_twist(path, 0.0) == 1.5
    @test material_twist(path, 0.7) == 1.5
    @test material_twist(path, 2.0) == 1.5
    @test total_material_twist(path) == 1.5 * 2.0   # exact, no tolerance
end

@testset "Twist — zero outside the run" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    straight!(spec; length = 1.0, meta = [Twist(; rate = 2.0)])
    straight!(spec; length = 1.0)
    path = build(spec)
    # Run is [1.0, 2.0] (until end of next twist anchor — but there is no
    # next anchor, so this run actually extends to s_end = 3.0).
    @test material_twist(path, 0.5) == 0.0
    @test material_twist(path, 1.5) == 2.0
    @test material_twist(path, 2.5) == 2.0   # extends to s_end
end

@testset "Twist — run terminates at next Twist anchor" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0, meta = [Twist(; rate = 1.0)])
    straight!(spec; length = 1.0, meta = [Twist(; rate = 3.0, is_continuous = true)])
    path = build(spec)
    @test length(path.resolved_twists) == 2
    @test path.resolved_twists[1].s_eff_end == 1.0
    @test path.resolved_twists[2].s_eff_start == 1.0
    @test material_twist(path, 0.5) == 1.0
    @test material_twist(path, 1.5) == 3.0
end

@testset "Twist — function rate is invariant under run-local s" begin
    f = s -> sin(s)
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)                                  # no twist
    straight!(spec; length = 2π, meta = [Twist(; rate = f)])
    path = build(spec)
    # At absolute s = 1.0 + 0.7, run-local s_local = 0.7
    @test material_twist(path, 1.7) == f(0.7)
    # Total over the run: ∫₀^{2π} sin(s) ds = 0
    @test isapprox(total_material_twist(path), 0.0; atol = 1e-7)
end

@testset "Twist — oscillatory rate handled by adaptive quadrature" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 2π,
              meta = [Twist(; rate = s -> sin(50 * s))])
    path = build(spec)
    # ∫₀^{2π} sin(50 s) ds = (1 - cos(100π)) / 50 = 0
    @test isapprox(total_material_twist(path), 0.0; atol = 1e-7)
end

@testset "Twist — phi_0 carry-over with is_continuous=true" begin
    spec = PathSpecBuilder()
    L1 = 1.5
    τ1 = 2.0
    straight!(spec; length = L1, meta = [Twist(; rate = τ1, phi_0 = 0.5)])
    straight!(spec; length = 1.0, meta = [Twist(; rate = 1.0, is_continuous = true)])
    path = build(spec)
    @test path.resolved_twists[1].phi_0 == 0.5
    @test isapprox(path.resolved_twists[2].phi_0, 0.5 + τ1 * L1; atol = 1e-12)
end

@testset "Twist — phi_0 carry-over with function rate" begin
    spec = PathSpecBuilder()
    L1 = π
    f1 = s -> cos(s)   # ∫₀^π cos(s) ds = sin(π) - sin(0) = 0
    straight!(spec; length = L1, meta = [Twist(; rate = f1, phi_0 = 0.7)])
    straight!(spec; length = 1.0, meta = [Twist(; rate = 1.0, is_continuous = true)])
    path = build(spec)
    @test isapprox(path.resolved_twists[2].phi_0, 0.7; atol = 1e-8)
end

@testset "Twist — total_material_twist partial interval" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 4.0, meta = [Twist(; rate = 0.5)])
    path = build(spec)
    @test total_material_twist(path; s_start = 1.0, s_end = 3.0) == 0.5 * 2.0
end

@testset "Twist — no anchors → zero everywhere" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    path = build(spec)
    @test path.resolved_twists == ResolvedTwistRate[]
    @test material_twist(path, 0.5) == 0.0
    @test total_material_twist(path) == 0.0
end

@testset "Twist — validation: first Twist with is_continuous=true rejected" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0,
              meta = [Twist(; rate = 1.0, is_continuous = true)])
    @test_throws ArgumentError build(spec)
end

@testset "Twist — validation: two Twists per segment rejected" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0,
              meta = [Twist(; rate = 1.0), Twist(; rate = 2.0)])
    @test_throws ArgumentError build(spec)
end

@testset "Twist — validation: phi_0 with is_continuous rejected at construction" begin
    @test_throws ArgumentError Twist(; rate = 1.0, phi_0 = 0.7, is_continuous = true)
end

@testset "Twist — frame() returns material_twist" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0, meta = [Twist(; rate = 2.5)])
    path = build(spec)
    @test frame(path, 0.4).material_twist == 2.5
end

@testset "Twist — total_frame_rotation = τ_geom + Ω_material" begin
    spec = PathSpecBuilder()
    # straight segment has τ_geom = 0, so total_frame_rotation = ∫τ_mat ds.
    straight!(spec; length = 2.0, meta = [Twist(; rate = 0.5)])
    path = build(spec)
    @test isapprox(total_frame_rotation(path), 1.0; atol = 1e-12)
end

@testset "Twist — path_twist_breakpoints includes run boundaries" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    straight!(spec; length = 1.0, meta = [Twist(; rate = 1.0)])
    straight!(spec; length = 1.0)
    path = build(spec)
    bps = path_twist_breakpoints(path)
    @test 0.0 in bps
    @test 1.0 in bps
    @test 3.0 in bps   # run extends to s_end since no later anchor
end

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

@testset "Path — total_turning_angle of a full circle" begin
    # T-PHYSICS: ∫κ ds over a full circle = 2π
    spec = PathSpecBuilder()
    bend!(spec; radius = 0.1, angle = 2π)
    path = build(spec)
    @test total_turning_angle(path) ≈ 2π atol = 1e-10
end

@testset "Path — total_torsion of straight and bend segments is zero" begin
    # T-PHYSICS: straight and circular-arc segments have zero geometric torsion
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.2, angle = π / 2)
    path = build(spec)
    @test total_torsion(path) == 0.0
end

@testset "Path — writhe of a straight path is zero" begin
    # T-PHYSICS: a straight line has no self-linking, Wr = 0
    spec = PathSpecBuilder()
    straight!(spec; length = 2.0)
    path = build(spec)
    @test abs(writhe(path; n = 64)) < 1e-6
end

# -----------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------

@testset "sample_uniform — returns n frames" begin
    # T-GUARDRAIL
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    bend!(spec; radius = 0.1, angle = π / 2)
    path = build(spec)

    frames = sample_uniform(path; n = 50)
    @test length(frames) == 50
    @test hasproperty(frames[1], :position)
    @test hasproperty(frames[1], :tangent)
    @test hasproperty(frames[1], :curvature)
end

# -----------------------------------------------------------------------
# JumpBy and JumpTo (QuinticConnector)
# -----------------------------------------------------------------------

@testset "JumpBy — endpoint matches delta in local frame" begin
    # T-PHYSICS: JumpBy with delta along the current tangent direction (ẑ)
    # should move the path forward by that amount.
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.0, 0.0, 0.5))
    path = build(spec)

    @test start_point(path) ≈ [0.0, 0.0, 0.0] atol = 1e-10
    @test end_point(path)   ≈ [0.0, 0.0, 0.5] atol = 1e-10
end

@testset "JumpBy — initial tangent is ẑ" begin
    # T-GUARDRAIL: incoming tangent at start of every connector is ẑ
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.1, 0.0, 0.3))
    path = build(spec)
    @test tangent(path, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "JumpBy — outgoing tangent honours tangent_out" begin
    # T-PHYSICS: explicit tangent_out should be the end tangent (in local frame)
    t_out = normalize([1.0, 0.0, 1.0])
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.3, 0.0, 0.3), tangent = t_out)
    path = build(spec)
    @test tangent(path, path.s_end) ≈ t_out atol = 1e-8
end

@testset "JumpBy — frame orthonormality along connector" begin
    # T-GUARDRAIL
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.2, 0.1, 0.4))
    path = build(spec)
    for s in range(path.spec.s_start, path.s_end; length = 11)
        T = tangent(path, s)
        N = normal(path, s)
        B = binormal(path, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "JumpBy — after straight, delta is in rotated local frame" begin
    # T-PHYSICS: after a quarter-circle bend the local ẑ points along global +x.
    # JumpBy with delta=(0,0,d) should therefore move +x in global frame.
    R = 0.1; d = 0.5
    spec = PathSpecBuilder()
    bend!(spec; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpby!(spec; delta = (0.0, 0.0, d))   # local ẑ is now global +x
    path = build(spec)
    # end_point of bend is (R, 0, R); jumpby adds d along global +x
    @test end_point(path) ≈ [R + d, 0.0, R] atol = 1e-8
end

@testset "JumpTo — endpoint matches destination" begin
    # T-PHYSICS: JumpTo should place the path end at the specified global position.
    dest = [1.0, 0.5, 2.0]
    spec = PathSpecBuilder()
    jumpto!(spec; destination = dest)
    path = build(spec)
    @test end_point(path) ≈ dest atol = 1e-10
end

@testset "JumpTo — initial tangent is ẑ" begin
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (0.3, 0.1, 0.8))
    path = build(spec)
    @test tangent(path, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "JumpTo — after bend, destination is in global frame" begin
    # T-PHYSICS: JumpTo destination is always in global frame regardless of prior segments.
    R = 0.1
    dest = [R + 0.5, 0.0, R]   # global position
    spec = PathSpecBuilder()
    bend!(spec; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpto!(spec; destination = dest)
    path = build(spec)
    @test end_point(path) ≈ collect(dest) atol = 1e-8
end

@testset "JumpTo — outgoing tangent honours tangent_out (global)" begin
    # T-PHYSICS: tangent_out for JumpTo is specified in global frame
    t_out_global = normalize([1.0, 0.0, 0.0])
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (1.0, 0.0, 0.5), tangent = t_out_global)
    path = build(spec)
    @test tangent(path, path.s_end) ≈ t_out_global atol = 1e-8
end

@testset "JumpBy/JumpTo — sample_path works on connector" begin
    # T-GUARDRAIL: sample_path must not error on paths containing connectors
    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    jumpby!(spec; delta = (0.1, 0.0, 0.4))
    straight!(spec; length = 0.2)
    path = build(spec)
    ps = sample_path(path, path.spec.s_start, path.s_end)
    @test ps.n >= 4
    @test ps.samples[1].s ≈ path.spec.s_start atol = 1e-12
    @test ps.samples[end].s ≈ path.s_end  atol = 1e-12
end

# -----------------------------------------------------------------------
# Per-segment meta bag
# -----------------------------------------------------------------------

if !isdefined(Main, :Nickname)
    include(joinpath(@__DIR__, "..", "src", "fiber", "fiber-path-meta.jl"))
end
if !isdefined(Main, :Fiber)
    include(joinpath(@__DIR__, "..", "src", "fiber", "fiber-path.jl"))
end
if !isdefined(Main, :modify)
    include(joinpath(@__DIR__, "..", "src", "fiber", "fiber-path-modify.jl"))
end

@testset "per-segment meta — builders forward meta to every segment type" begin
    nick = [Nickname("alpha")]
    mcm  = [MCMadd(:T_K, (:Normal, 0.0, 1.0))]

    spec = PathSpecBuilder()
    straight!(spec; length = 0.1, meta = nick)
    bend!(spec; radius = 0.05, angle = π / 2, meta = mcm)
    helix!(spec; radius = 0.02, pitch = 0.01, turns = 1.0,
           meta = [Nickname("helix"), MCMadd(:T_K, :stub)])
    catenary!(spec; a = 0.04, length = 0.05, meta = [Nickname("cat")])
    jumpby!(spec; delta = (0.0, 0.0, 0.05), meta = [Nickname("jb")])
    jumpto!(spec; destination = (0.0, 0.1, 0.4), meta = [Nickname("jt")])

    segs = spec.segments
    @test segs[1].meta == nick
    @test segs[2].meta == mcm
    @test length(segs[3].meta) == 2
    @test segs[4].meta[1] isa Nickname
    @test segs[5].meta[1] isa Nickname
    @test segs[6].meta[1] isa Nickname

    @test segment_meta(segs[1]) === segs[1].meta
    @test segment_nickname(segs[1]) == "alpha"
    @test isnothing(segment_nickname(segs[2]))  # mcm only
end

@testset "per-segment meta — build() copies jump meta onto QuinticConnector" begin
    spec = PathSpecBuilder()
    straight!(spec; length = 0.1)
    jumpby!(spec; delta = (0.05, 0.0, 0.2),
            meta = [Nickname("connector-1"),
                                MCMadd(:T_K, :stub)])
    path = build(spec)

    hc = path.placed_segments[2].segment
    @test hc isa QuinticConnector
    @test length(hc.meta) == 2
    @test segment_nickname(hc) == "connector-1"
    @test any(m -> m isa MCMadd, hc.meta)
end

@testset "per-segment meta — shrinkage preserves meta" begin
    xs = FiberCrossSection(GermaniaSilicaGlass(0.036), GermaniaSilicaGlass(0.0),
                           8.2e-6, 125e-6)

    spec = PathSpecBuilder()
    straight!(spec; length = 0.1,
              meta = [Nickname("s"), MCMadd(:T_K, 10.0)])
    bend!(spec; radius = 0.05, angle = π / 3,
          meta = [Nickname("b")])
    helix!(spec; radius = 0.02, pitch = 0.01, turns = 1.0,
           meta = [Nickname("h")])
    catenary!(spec; a = 0.03, length = 0.04,
              meta = [Nickname("c")])
    jumpby!(spec; delta = (0.0, 0.0, 0.05),
            meta = [Nickname("j")])
    path = build(spec)
    fiber = Fiber(path; cross_section = xs, T_ref_K = 297.15)
    path2 = modify(fiber)

    @test segment_nickname(path2.placed_segments[1].segment) == "s"
    @test any(m -> m isa MCMadd, path2.placed_segments[1].segment.meta)
    @test segment_nickname(path2.placed_segments[2].segment) == "b"
    @test segment_nickname(path2.placed_segments[3].segment) == "h"
    @test segment_nickname(path2.placed_segments[4].segment) == "c"
    @test segment_nickname(path2.placed_segments[5].segment) == "j"
end

@testset "per-segment meta — segment_meta returns empty default" begin
    seg = StraightSegment(0.1)
    @test segment_meta(seg) == AbstractMeta[]
    @test isnothing(segment_nickname(seg))
end

#=
DISABLED: min_bend_radius enforcement is a stub.

The enforcement code in _build_hermite_connector is commented out (see
path-geometry.jl, inside _build_hermite_connector). Re-enable these testsets
when the scan/bisect block is restored. See commit f06e689.

@testset "JumpBy — min_bend_radius keeps curvature ≤ 1/R_min" begin
    # T-PHYSICS: every point along the connector must have κ ≤ 1/R_min.
    R_min = 0.05
    spec = PathSpecBuilder()
    jumpby!(spec; delta = (0.02, 0.0, 0.03), min_bend_radius = R_min)
    path = build(spec)
    for s in range(path.spec.s_start, path.s_end; length = 51)
        @test curvature(path, s) ≤ 1.0 / R_min
    end
end

@testset "JumpTo — min_bend_radius keeps curvature ≤ 1/R_min" begin
    # T-PHYSICS: min_bend_radius on JumpTo also constrains the connector curvature.
    R_min = 0.04
    spec = PathSpecBuilder()
    jumpto!(spec; destination = (0.05, 0.0, 0.06), min_bend_radius = R_min)
    path = build(spec)
    for s in range(path.spec.s_start, path.s_end; length = 51)
        @test curvature(path, s) ≤ 1.0 / R_min
    end
end

@testset "JumpBy — min_bend_radius does not move endpoint" begin
    # T-GUARDRAIL: the constraint changes handle length only, not the target position.
    delta = (0.03, 0.0, 0.05)
    spec = PathSpecBuilder()
    jumpby!(spec; delta = delta, min_bend_radius = 0.10)
    path = build(spec)
    @test end_point(path) ≈ collect(delta) atol = 1e-10
end

@testset "JumpTo — infeasible min_bend_radius throws ArgumentError" begin
    # T-GUARDRAIL: anti-parallel incoming/outgoing tangents cause peak curvature to
    # increase (not decrease) with handle length — no equal-handle Hermite can satisfy
    # the constraint, so build() must throw ArgumentError.
    spec = PathSpecBuilder()
    straight!(spec; length = 1.0)
    # After straight, incoming tangent is +z. destination=(2,0,1) is 1 m transverse.
    # Outgoing tangent is -z (anti-parallel). κ grows with h → infeasible.
    jumpto!(spec; destination = (2.0, 0.0, 1.0), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 2.0)
    @test_throws ArgumentError build(spec)
end

# JumpTo min_bend_radius — T1/T2/T3 from demo_fiber_path_jumps_min_radius
# Each helper builds the path up to and including the named jump.
# Thresholds are taken from the demo comments.

@testset "JumpTo min_bend_radius — T1: passes at 0.5, fails above 0.5" begin
    # T-PHYSICS: transverse chord (1,0,0), outgoing tangent (0,0,-1) anti-parallel
    # to incoming (0,0,1). Minimum-curvature connector has bend radius ≈ 0.5 m.
    function build_T1(mbr)
        spec = PathSpecBuilder()
        straight!(spec; length = 1)
        jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = mbr)
        build(spec)
    end
    @test_nowarn build_T1(0.49)
    @test_nowarn build_T1(0.5)
    @test_throws ArgumentError build_T1(0.51)
end

@testset "JumpTo min_bend_radius — T2: passes at 0.5, fails above 0.5" begin
    # T-PHYSICS: after T1 the path runs downward (-z). The straight extends
    # that to position (1,0,0). T2 jumps to (2,0,0) with outgoing tangent (0,0,1)
    # — again transverse chord, anti-parallel tangents.
    function build_T2(mbr)
        spec = PathSpecBuilder()
        straight!(spec; length = 1)
        jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = 0.1) # so small it won't error
        straight!(spec; length = 1)
        jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
                min_bend_radius = mbr)
        build(spec)
    end
    @test_nowarn build_T2(0.49)
    @test_nowarn build_T2(0.5)
    @test_throws ArgumentError build_T2(0.51)
end

@testset "JumpTo min_bend_radius — T3: passes at 0.50, fails at 0.51" begin
    # T-PHYSICS: same pattern as T1; threshold is 0.50 m.
    function build_T3(mbr)
        spec = PathSpecBuilder()
        straight!(spec; length = 1)
        jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = 0.1) # so small it won't error
        straight!(spec; length = 1)
        jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
                min_bend_radius = 0.1) # so small it won't error
        straight!(spec; length = 1)
        jumpto!(spec; destination = (3, 0.0, 1), tangent = (0.0, 0.0, -1.0),
                min_bend_radius = mbr)
        build(spec)
    end
    @test_nowarn build_T3(0.49)
    @test_nowarn build_T3(0.50)
    @test_throws ArgumentError build_T3(0.51)
end
=#


# -----------------------------------------------------------------------
# G2 integration: curvature_out propagation
# -----------------------------------------------------------------------

@testset "JumpBy — curvature_out matches sampled κ at end of connector" begin
    """
    G2 outgoing match: when the user specifies curvature_out, the connector's
    sampled scalar curvature at its endpoint must equal ‖curvature_out‖.

    Confirms that the new curvature_out field on JumpBy is plumbed through
    _resolve_at_placement into _build_quintic_connector and survives all the
    way to the curve's terminal κ-sample.
    """
    spec = PathSpecBuilder()
    K1 = (0.0, 2.0, 0.0)   # 2 m⁻¹ in local +y
    jumpby!(spec; delta = (0.5, 0.0, 0.5),
            tangent = (1.0, 0.0, 0.0),
            curvature_out = K1)
    path = build(spec)
    seg = path.placed_segments[1].segment
    L = arc_length(seg)
    @test isapprox(curvature(seg, L), sqrt(K1[1]^2 + K1[2]^2 + K1[3]^2);
                   rtol = 1e-3, atol = 1e-6)
end

@testset "JumpBy — incoming K0 inherited from prior bend (G2 join)" begin
    """
    G2 incoming match: when a JumpBy follows a BendSegment, the connector's
    sampled scalar curvature at its *start* must equal the bend's curvature
    (1/R_bend). This verifies the K_in_global threading in build() and the
    frame-rotation logic in _resolve_at_placement.

    Without this wiring, the connector would start with κ=0 (G1 join) and the
    splice would have a curvature jump.
    """
    R_bend = 0.5
    spec = PathSpecBuilder()
    bend!(spec; radius = R_bend, angle = π/4)
    jumpby!(spec; delta = (0.3, 0.0, 0.3))
    path = build(spec)
    seg = path.placed_segments[2].segment
    @test isapprox(curvature(seg, 0.0), 1.0 / R_bend; rtol = 1e-2, atol = 1e-4)
end

@testset "JumpTo — global-frame curvature_out after bend" begin
    """
    JumpTo with curvature_out specified in the *global* frame, placed after a
    BendSegment that has rotated the local frame. The connector's terminal
    curvature scalar must equal ‖K1_global‖ regardless of frame rotation.

    Probes the global→local rotation of curvature_out in
    _resolve_at_placement(::JumpTo, …).
    """
    spec = PathSpecBuilder()
    bend!(spec; radius = 0.3, angle = π/2)   # rotates local frame
    K1_global = (0.0, 1.0, 0.0)
    jumpto!(spec; destination = (0.6, 0.0, 0.6),
            tangent = (0.0, 0.0, 1.0),
            curvature_out = K1_global)
    path = build(spec)
    seg = path.placed_segments[2].segment
    L = arc_length(seg)
    @test isapprox(curvature(seg, L), 1.0; rtol = 1e-3, atol = 1e-6)
end
