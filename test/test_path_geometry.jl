using Test
using LinearAlgebra
using Bifrost
using Bifrost.PathGeometry: _qc_nominalize


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

is_unit(v)       = abs(norm(v) - 1.0) < 1e-12
is_orthonormal(T, N, B) = is_unit(T) && is_unit(N) && is_unit(B) &&
                           abs(dot(T, N)) < 1e-12 &&
                           abs(dot(T, B)) < 1e-12 &&
                           abs(dot(N, B)) < 1e-12 &&
                           norm(cross(T, N) - B) < 1e-12   # right-handed

# End-of-interior arc-length coordinate. Equals the old `path.s_end` for
# Subpaths whose terminal `jumpto!` sits at the natural endpoint of the
# interior geometry with a matching `incoming_tangent` (so the connector is
# of negligible length).
_s_end_interior(b::SubpathBuilt) = Float64(_qc_nominalize(b.jumpto_placed.s_offset_eff))

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
# SubpathBuilder lifecycle and Subpath construction (T-GUARDRAIL)
# -----------------------------------------------------------------------

@testset "SubpathBuilder — segment before start! is rejected" begin
    sb = SubpathBuilder()
    @test_throws ArgumentError straight!(sb; length = 1.0)
end

@testset "SubpathBuilder — start! after segments is rejected" begin
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 1.0)
    @test_throws ArgumentError start!(sb)
end

@testset "SubpathBuilder — second start! is rejected" begin
    sb = SubpathBuilder()
    start!(sb)
    @test_throws ArgumentError start!(sb)
end

@testset "SubpathBuilder — jumpto! before start! is rejected" begin
    sb = SubpathBuilder()
    @test_throws ArgumentError jumpto!(sb; point = (0,0,1))
end

@testset "SubpathBuilder — second jumpto! is rejected" begin
    sb = SubpathBuilder()
    start!(sb)
    jumpto!(sb; point = (0,0,1))
    @test_throws ArgumentError jumpto!(sb; point = (0,0,2))
end

@testset "SubpathBuilder — segment after jumpto! is rejected" begin
    sb = SubpathBuilder()
    start!(sb)
    jumpto!(sb; point = (0,0,1))
    @test_throws ArgumentError straight!(sb; length = 1.0)
end

@testset "Subpath — missing start! rejected at construction" begin
    sb = SubpathBuilder()
    @test_throws ArgumentError Subpath(sb)
end

@testset "Subpath — missing jumpto! rejected at construction" begin
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 1.0)
    @test_throws ArgumentError Subpath(sb)
end

@testset "Subpath — geometry queries on unbuilt Subpath throw" begin
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 1.0)
    jumpto!(sb; point = (0,0,1.0), incoming_tangent = (0,0,1.0))
    sub = Subpath(sb)
    @test_throws ErrorException arc_length(sub)
    @test_throws ErrorException curvature(sub, 0.5)
    @test_throws ErrorException position(sub, 0.5)
    @test_throws ErrorException tangent(sub, 0.5)
end

# -----------------------------------------------------------------------
# seal! — natural-exit seal (no terminal connector bending)
# -----------------------------------------------------------------------

@testset "seal! — satisfies Subpath construction (no jumpto! needed)" begin
    # T-GUARDRAIL: a builder sealed by seal! is constructible; the seal
    # contract is satisfiable by seal! as well as jumpto!.
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 1.0)
    seal!(sb)
    @test Subpath(sb) isa Subpath
end

@testset "seal! — before start! is rejected" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    @test_throws ArgumentError seal!(sb)
end

@testset "seal! — after jumpto! is rejected (double seal)" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    start!(sb)
    jumpto!(sb; point = (0,0,1))
    @test_throws ArgumentError seal!(sb)
end

@testset "seal! — jumpto! after seal! is rejected (double seal)" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    start!(sb)
    seal!(sb)
    @test_throws ArgumentError jumpto!(sb; point = (0,0,1))
end

@testset "seal! — second seal! is rejected" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    start!(sb)
    seal!(sb)
    @test_throws ArgumentError seal!(sb)
end

@testset "seal! — segment after seal! is rejected" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    start!(sb)
    seal!(sb)
    @test_throws ArgumentError straight!(sb; length = 1.0)
end

@testset "seal! — negative extra is rejected" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    start!(sb)
    @test_throws ArgumentError seal!(sb; extra = -0.1)
end

@testset "seal! — preserves the natural exit (no bending, no added length)" begin
    # T-PHYSICS: with extra = 0 the terminal connector is zero length, so the
    # path ends exactly at the last interior segment's endpoint, with its
    # tangent, and total arc length equals the sum of interior lengths.
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 0.10)
    bend!(sb; radius = 0.05, angle = π / 2)   # quarter circle in x-z plane
    seal!(sb)
    b = build(Subpath(sb))

    L_interior = 0.10 + 0.05 * (π / 2)        # straight + arc length
    s_e = Float64(_qc_nominalize(s_end(b)))
    @test isapprox(s_e, L_interior; atol = 1e-9)

    # Quarter circle from tangent +z about axis_angle=0 ends heading +x at
    # (0.05, 0, 0.10 + 0.05).
    @test isapprox(collect(position(b, s_e)), [0.05, 0.0, 0.15]; atol = 1e-8)
    @test isapprox(collect(tangent(b, s_e)),  [1.0, 0.0, 0.0];  atol = 1e-8)
end

@testset "seal! — extra appends a straight lead-out along the exit tangent" begin
    # T-PHYSICS: extra = L adds exactly L of straight length along the natural
    # exit tangent and nothing else.
    L_extra = 0.07
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 0.10)
    bend!(sb; radius = 0.05, angle = π / 2)
    seal!(sb; extra = L_extra)
    b = build(Subpath(sb))

    L_interior = 0.10 + 0.05 * (π / 2)
    s_e = Float64(_qc_nominalize(s_end(b)))
    @test isapprox(s_e, L_interior + L_extra; atol = 1e-9)

    # Exit was heading +x at (0.05, 0, 0.15); lead-out advances +x by L_extra.
    @test isapprox(collect(position(b, s_e)), [0.05 + L_extra, 0.0, 0.15];
                   atol = 1e-8)
    @test isapprox(collect(tangent(b, s_e)),  [1.0, 0.0, 0.0]; atol = 1e-8)
    # The lead-out is straight: zero curvature near the end.
    @test isapprox(curvature(b, s_e - 1e-4), 0.0; atol = 1e-6)
end

@testset "seal! — naturally-sealed predecessor stitches into a PathBuilt" begin
    # T-GUARDRAIL: conformity check reads the built endpoint of a seal!'d
    # predecessor (jumpto_point === nothing) rather than a global jumpto spec.
    sb1 = SubpathBuilder()
    start!(sb1)
    straight!(sb1; length = 0.2)
    seal!(sb1)                                # natural exit at (0,0,0.2), +z

    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 0.2), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.1)
    seal!(sb2)

    p = build([Subpath(sb1), Subpath(sb2)])
    @test isapprox(s_end(p), 0.3; atol = 1e-9)
end

@testset "seal! — mismatched start after natural seal is rejected" begin
    # T-GUARDRAIL: the built endpoint must still match the next start.
    sb1 = SubpathBuilder()
    start!(sb1)
    straight!(sb1; length = 0.2)
    seal!(sb1)                                # ends at (0,0,0.2)

    sb2 = SubpathBuilder()
    start!(sb2; point = (1.0, 0.0, 0.2), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.1)
    seal!(sb2)

    @test_throws ArgumentError build([Subpath(sb1), Subpath(sb2)])
end

@testset "build(Vector{SubpathBuilder}) — matches explicit Subpath freeze" begin
    # T-GUARDRAIL: build([sb1, sb2]) is a convenience that freezes each builder
    # to a Subpath, so it must agree with build([Subpath(sb1), Subpath(sb2)]).
    make_builders() = begin
        sb1 = SubpathBuilder()
        start!(sb1)
        straight!(sb1; length = 0.2)
        seal!(sb1)                                # natural exit at (0,0,0.2), +z

        sb2 = SubpathBuilder()
        start!(sb2; point = (0.0, 0.0, 0.2), outgoing_tangent = (0.0, 0.0, 1.0))
        straight!(sb2; length = 0.1)
        seal!(sb2)
        (sb1, sb2)
    end

    sb1, sb2 = make_builders()
    p_conv = build([sb1, sb2])
    @test p_conv isa PathBuilt

    sb1e, sb2e = make_builders()
    p_ref = build([Subpath(sb1e), Subpath(sb2e)])

    @test length(p_conv.subpaths) == length(p_ref.subpaths) == 2
    @test isapprox(s_end(p_conv), s_end(p_ref); atol = 1e-12)
    for s in (0.0, 0.15, s_end(p_ref))
        @test isapprox(position(p_conv, s), position(p_ref, s); atol = 1e-10)
    end
end

@testset "build(Vector{SubpathBuilder}) — conformity errors still surface" begin
    # T-GUARDRAIL: a mismatched start between builders must still throw, exactly
    # as it does for the explicit Subpath form.
    sb1 = SubpathBuilder()
    start!(sb1)
    straight!(sb1; length = 0.2)
    seal!(sb1)                                # ends at (0,0,0.2)

    sb2 = SubpathBuilder()
    start!(sb2; point = (1.0, 0.0, 0.2), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.1)
    seal!(sb2)

    @test_throws ArgumentError build([sb1, sb2])
end

# -----------------------------------------------------------------------
# :inherit start-state
# -----------------------------------------------------------------------

# Build the standard transverse-chord predecessor used by several inherit tests:
# a straight up to (0,0,1) sealed by jumpto! landing at (1,0,1) heading -z.
function _inherit_predecessor()
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 1.0)
    jumpto!(sb; point = (1.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.4)
    return sb
end

@testset "inherit — :inherit reproduces hand-loaded coordinates" begin
    # T-SIM-REGRESSION: start!(sb2, :inherit) must yield byte-identical geometry
    # to hand-loading the predecessor endpoint (1,0,1) with tangent -z.
    sb2h = SubpathBuilder()
    start!(sb2h; point = (1.0, 0.0, 1.0), outgoing_tangent = (0.0, 0.0, -1.0))
    straight!(sb2h; length = 1.0)
    seal!(sb2h)
    p_hand = build([Subpath(_inherit_predecessor()), Subpath(sb2h)])

    sb2i = SubpathBuilder()
    start!(sb2i, :inherit)
    straight!(sb2i; length = 1.0)
    seal!(sb2i)
    p_inh = build([Subpath(_inherit_predecessor()), Subpath(sb2i)])

    s_hi = Float64(_qc_nominalize(path_length(p_hand)))
    for s in range(0.0, s_hi; length = 21)
        @test position(p_hand, s) ≈ position(p_inh, s) atol = 1e-12
        @test tangent(p_hand, s)  ≈ tangent(p_inh, s)  atol = 1e-12
    end
end

@testset "inherit — per-field :inherit mixes with explicit values" begin
    # T-SIM-REGRESSION: inherit point only, set tangent explicitly to the same
    # value the predecessor exits with; result equals the all-inherit build.
    sb2m = SubpathBuilder()
    start!(sb2m; point = :inherit, outgoing_tangent = (0.0, 0.0, -1.0))
    straight!(sb2m; length = 1.0)
    seal!(sb2m)
    p_mix = build([Subpath(_inherit_predecessor()), Subpath(sb2m)])

    sb2i = SubpathBuilder()
    start!(sb2i, :inherit)
    straight!(sb2i; length = 1.0)
    seal!(sb2i)
    p_inh = build([Subpath(_inherit_predecessor()), Subpath(sb2i)])

    s_hi = Float64(_qc_nominalize(path_length(p_inh)))
    for s in range(0.0, s_hi; length = 11)
        @test position(p_mix, s) ≈ position(p_inh, s) atol = 1e-12
    end
end

@testset "inherit — tangent from chord-default predecessor" begin
    # T-SIM-REGRESSION: predecessor jumpto! with incoming_tangent=nothing exits
    # along the chord direction; :inherit must resolve that concrete tangent
    # (querying the built geometry), not leave it `nothing`. Check the resolved
    # start tangent directly via the resolver.
    pred = SubpathBuilder()
    start!(pred)
    straight!(pred; length = 1.0)
    jumpto!(pred; point = (0.5, 0.0, 1.5))    # no incoming_tangent → chord dir
    pred_built = build(Subpath(pred))

    sb2 = SubpathBuilder()
    start!(sb2, :inherit)
    straight!(sb2; length = 0.3)
    seal!(sb2)
    resolved = Bifrost.PathGeometry._resolve_inherited_start(Subpath(sb2), pred_built)
    @test !resolved.inherit_start_tangent
    @test collect(resolved.start_outgoing_tangent) ≈ collect(end_tangent(pred_built)) atol = 1e-9
    # And the full join must build without error.
    @test build([Subpath(pred), Subpath(sb2)]) isa PathBuilt
end

@testset "inherit — curvature inherits declared jumpto_incoming_curvature" begin
    # T-SIM-REGRESSION: when the predecessor declares incoming_curvature, an
    # :inherit start picks it up; otherwise it defaults to (0,0,0).
    pred = SubpathBuilder()
    start!(pred)
    bend!(pred; radius = 0.5, angle = π / 2)
    # End of a radius-0.5 quarter bend: declare a matching incoming curvature so
    # the join is G2 and inheritance has a non-zero value to copy.
    jumpto!(pred; point = (0.5, 0.0, 0.5), incoming_tangent = (1.0, 0.0, 0.0),
            incoming_curvature = (0.0, 0.0, -2.0))
    sb2 = SubpathBuilder()
    start!(sb2, :inherit)
    straight!(sb2; length = 0.2)
    seal!(sb2)
    resolved = Bifrost.PathGeometry._resolve_inherited_start(
        Subpath(sb2), build(Subpath(pred)))
    @test resolved.start_outgoing_curvature == (0.0, 0.0, -2.0)
    @test !resolved.inherit_start_curvature
end

@testset "inherit — first Subpath with :inherit is rejected" begin
    # T-GUARDRAIL: there is no predecessor to inherit from.
    sb1 = SubpathBuilder()
    start!(sb1, :inherit)
    straight!(sb1; length = 1.0)
    seal!(sb1)
    sb2 = SubpathBuilder()
    start!(sb2, :inherit)
    straight!(sb2; length = 1.0)
    seal!(sb2)
    @test_throws ArgumentError build([Subpath(sb1), Subpath(sb2)])
end

@testset "inherit — standalone build of an inherit Subpath is rejected" begin
    # T-GUARDRAIL: an unresolved :inherit Subpath cannot be placed alone.
    sb = SubpathBuilder()
    start!(sb, :inherit)
    straight!(sb; length = 1.0)
    seal!(sb)
    @test_throws ArgumentError build(Subpath(sb))
    @test_throws ArgumentError build(sb)
end

@testset "inherit — non-:inherit symbol is rejected" begin
    # T-GUARDRAIL: only :inherit is accepted as a symbol.
    sb = SubpathBuilder()
    @test_throws ArgumentError start!(sb, :nope)
    sb2 = SubpathBuilder()
    @test_throws ArgumentError start!(sb2; point = :nope)
end

@testset "inherit — non-inherit Vector{Subpath} build is unchanged" begin
    # T-SIM-REGRESSION: a hand-loaded multi-Subpath build (no :inherit) joins at
    # the declared coordinates.
    sb2 = SubpathBuilder()
    start!(sb2; point = (1.0, 0.0, 1.0), outgoing_tangent = (0.0, 0.0, -1.0))
    straight!(sb2; length = 1.0)
    seal!(sb2)
    p = build([Subpath(_inherit_predecessor()), Subpath(sb2)])
    # Predecessor ends at (1,0,1) heading -z; a 1 m straight lands at (1,0,0).
    @test end_point(p) ≈ [1.0, 0.0, 0.0] atol = 1e-8
end

# -----------------------------------------------------------------------
# Subpath assembly and build
# -----------------------------------------------------------------------

@testset "Subpath — single straight segment" begin
    # T-PHYSICS: a single straight segment of length L produces a path from
    # (0,0,0) to (0,0,L) with constant tangent (0,0,1).
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 2.0)
    jumpto!(sb; point = (0.0, 0.0, 2.0), incoming_tangent = (0.0, 0.0, 1.0))
    b = build(Subpath(sb))

    s_int = _s_end_interior(b)
    @test s_int ≈ 2.0
    @test position(b, 0.0)   ≈ [0.0, 0.0, 0.0] atol = 1e-12
    @test position(b, s_int) ≈ [0.0, 0.0, 2.0] atol = 1e-10
    @test tangent(b, 0.0)    ≈ [0.0, 0.0, 1.0] atol = 1e-12
    @test tangent(b, s_int)  ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "Path — tangent continuity at segment joints" begin
    # T-GUARDRAIL: tangent must be continuous across every segment boundary.
    # Test with a straight → bend → straight sequence.
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 0.5)
    bend!(sb; radius = 0.1, angle = π / 2, axis_angle = 0.0)
    straight!(sb; length = 0.3)
    # End of the third straight: position (0.1 + 0.3, 0, 0.5 + 0.1) = (0.4, 0, 0.6),
    # tangent (1, 0, 0).
    jumpto!(sb; point = (0.4, 0.0, 0.6), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)

    ps = b.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = ps[i + 1].s_offset_eff
        T_before = tangent(b, s_joint - 1e-9)
        T_after  = tangent(b, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "Path — straight + quarter bend geometry" begin
    # T-PHYSICS: straight segment of length L followed by quarter-circle of
    # radius R. End position should be (R, 0, L+R), end tangent should be (1,0,0).
    L = 1.0; R = 0.2
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = L)
    bend!(sb; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpto!(sb; point = (R, 0.0, L + R), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)

    s_int = _s_end_interior(b)
    @test position(b, s_int) ≈ [R, 0.0, L + R] atol = 1e-10
    @test tangent(b, s_int)  ≈ [1.0, 0.0, 0.0]  atol = 1e-10
end

@testset "Path — full circle returns to start" begin
    # T-PHYSICS: a complete circle of radius R returns to the start position
    # with the same tangent direction.
    R = 0.15
    sb = SubpathBuilder()
    start!(sb)
    bend!(sb; radius = R, angle = 2π, axis_angle = 0.0)
    # After a full circle the natural endpoint is back at the origin with
    # incoming tangent (0,0,1).
    jumpto!(sb; point = (0.0, 0.0, 0.0), incoming_tangent = (0.0, 0.0, 1.0))
    b = build(sb)

    s_int = _s_end_interior(b)
    @test position(b, 0.0)    ≈ position(b, s_int) atol = 1e-10
    @test tangent(b, 0.0)     ≈ tangent(b, s_int)   atol = 1e-10
    @test s_int               ≈ 2π * R              atol = 1e-10
end

@testset "Path — cartesian_distance vs arc_length" begin
    # T-PHYSICS: for a straight segment, cartesian_distance == arc_length.
    # For a curved path, cartesian_distance < arc_length.
    sb_s = SubpathBuilder()
    start!(sb_s)
    straight!(sb_s; length = 3.0)
    jumpto!(sb_s; point = (0.0, 0.0, 3.0), incoming_tangent = (0.0, 0.0, 1.0))
    b_s = build(sb_s)
    s_int_s = _s_end_interior(b_s)
    @test cartesian_distance(b_s, 0.0, s_int_s) ≈ 3.0 atol = 1e-10

    sb_b = SubpathBuilder()
    start!(sb_b)
    bend!(sb_b; radius = 0.5, angle = π)
    # Half-circle: end at (1.0, 0, 0) with incoming tangent (0,0,-1).
    jumpto!(sb_b; point = (1.0, 0.0, 0.0), incoming_tangent = (0.0, 0.0, -1.0))
    b_b = build(sb_b)
    L_int = _s_end_interior(b_b)
    @test cartesian_distance(b_b, 0.0, L_int) < L_int   # chord < arc
end

@testset "Path — frame orthonormality along assembled path" begin
    # T-GUARDRAIL: frame must remain orthonormal everywhere in the assembled path
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 0.3)
    bend!(sb; radius = 0.08, angle = π / 3, axis_angle = π / 6)
    catenary!(sb; a = 0.2, length = 0.4)
    # Compute the natural exit position/tangent by building once with a
    # placeholder jumpto, reading position/tangent from the last interior
    # segment, and re-authoring with the correct jumpto. Simpler alternative:
    # just iterate on s ∈ [0, _s_end_interior(b)] (which is set by the
    # interior segments alone) and accept whatever connector geometry exists.
    jumpto!(sb; point = (0.0, 0.0, 0.0))   # placeholder; orthonormality is
                                           # checked over the *interior* range.
    b = build(sb)

    ss = range(0.0, _s_end_interior(b); length = 31)
    for s in ss
        T = tangent(b, s)
        N = normal(b, s)
        B = binormal(b, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "Path — bounding box contains all sampled points" begin
    # T-GUARDRAIL
    sb = SubpathBuilder()
    start!(sb)
    straight!(sb; length = 1.0)
    bend!(sb; radius = 0.2, angle = π / 2)
    # End of bend: (0.2, 0, 1.0+0.2) with tangent (1,0,0).
    jumpto!(sb; point = (0.2, 0.0, 1.2), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)
    bb = bounding_box(b; n = 256)

    for s in range(0.0, _s_end_interior(b); length = 64)
        p = position(b, s)
        @test all(p .>= bb.lo .- 1e-10)
        @test all(p .<= bb.hi .+ 1e-10)
    end
end

# -----------------------------------------------------------------------
# Spinning meta and spinning_rate
# -----------------------------------------------------------------------

# Helper: terminate a straight-only spec at its natural endpoint with
# incoming_tangent (0,0,1) so the connector is degenerate.
function _seal_at_z(sb::SubpathBuilder, z::Real)
    jumpto!(sb; point = (0.0, 0.0, z), incoming_tangent = (0.0, 0.0, 1.0))
end

@testset "Spinning — constant rate (Float64) is exact" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 2.0, meta = [Spinning(; rate = 1.5, phi_0 = 0.0)])
    _seal_at_z(sb, 2.0)
    b = build(sb)
    @test spinning_rate(b, 0.0) == 1.5
    @test spinning_rate(b, 0.7) == 1.5
    @test spinning_rate(b, 2.0) == 1.5
    # The spinning run extends from s=0 to the end of the Subpath. Over the
    # interior [0, 2.0] the integrated spinning is exactly 1.5*2.0; the
    # connector is degenerate so the full integral is the same to numerical
    # precision.
    @test isapprox(total_spinning(b; s_start = 0.0, s_end = 2.0), 1.5 * 2.0;
                   atol = 1e-12)
end

@testset "Spinning — zero outside the run" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 2.0)])
    straight!(sb; length = 1.0)
    _seal_at_z(sb, 3.0)
    b = build(sb)
    # Run starts at s=1.0 (its anchor segment's offset) and extends to s_end
    # of the Subpath since there is no later anchor.
    @test spinning_rate(b, 0.5) == 0.0
    @test spinning_rate(b, 1.5) == 2.0
    @test spinning_rate(b, 2.5) == 2.0   # extends to s_end
end

@testset "Spinning — run terminates at next Spinning anchor" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 1.0)])
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 3.0, is_continuous = true)])
    _seal_at_z(sb, 2.0)
    b = build(sb)
    @test length(b.resolved_spinning) == 2
    @test b.resolved_spinning[1].s_eff_end == 1.0
    @test b.resolved_spinning[2].s_eff_start == 1.0
    @test spinning_rate(b, 0.5) == 1.0
    @test spinning_rate(b, 1.5) == 3.0
end

@testset "Spinning — function rate is invariant under run-local s" begin
    f = s -> sin(s)
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)                                    # no spinning
    straight!(sb; length = 2π, meta = [Spinning(; rate = f)])
    _seal_at_z(sb, 1.0 + 2π)
    b = build(sb)
    # At absolute s = 1.0 + 0.7, run-local s_local = 0.7
    @test spinning_rate(b, 1.7) == f(0.7)
    # Total over the run (over the interior portion): ∫₀^{2π} sin(s) ds = 0
    @test isapprox(total_spinning(b; s_start = 1.0, s_end = 1.0 + 2π), 0.0;
                   atol = 1e-7)
end

@testset "Spinning — oscillatory rate handled by adaptive quadrature" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 2π,
              meta = [Spinning(; rate = s -> sin(50 * s))])
    _seal_at_z(sb, 2π)
    b = build(sb)
    # ∫₀^{2π} sin(50 s) ds = (1 - cos(100π)) / 50 = 0
    @test isapprox(total_spinning(b; s_start = 0.0, s_end = 2π), 0.0;
                   atol = 1e-7)
end

@testset "Spinning — phi_0 carry-over with is_continuous=true" begin
    sb = SubpathBuilder(); start!(sb)
    L1 = 1.5
    τ1 = 2.0
    straight!(sb; length = L1, meta = [Spinning(; rate = τ1, phi_0 = 0.5)])
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 1.0, is_continuous = true)])
    _seal_at_z(sb, L1 + 1.0)
    b = build(sb)
    @test b.resolved_spinning[1].phi_0 == 0.5
    @test isapprox(b.resolved_spinning[2].phi_0, 0.5 + τ1 * L1; atol = 1e-12)
end

@testset "Spinning — phi_0 carry-over with function rate" begin
    sb = SubpathBuilder(); start!(sb)
    L1 = π
    f1 = s -> cos(s)   # ∫₀^π cos(s) ds = sin(π) - sin(0) = 0
    straight!(sb; length = L1, meta = [Spinning(; rate = f1, phi_0 = 0.7)])
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 1.0, is_continuous = true)])
    _seal_at_z(sb, L1 + 1.0)
    b = build(sb)
    @test isapprox(b.resolved_spinning[2].phi_0, 0.7; atol = 1e-8)
end

@testset "Spinning — total_spinning partial interval" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 4.0, meta = [Spinning(; rate = 0.5)])
    _seal_at_z(sb, 4.0)
    b = build(sb)
    @test total_spinning(b; s_start = 1.0, s_end = 3.0) == 0.5 * 2.0
end

@testset "Spinning — no anchors → zero everywhere" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)
    _seal_at_z(sb, 1.0)
    b = build(sb)
    @test b.resolved_spinning == ResolvedSpinningRate[]
    @test spinning_rate(b, 0.5) == 0.0
    @test total_spinning(b; s_start = 0.0, s_end = 1.0) == 0.0
end

@testset "Spinning — first Spinning with is_continuous=true on Subpath_1 throws at PathBuilt level" begin
    # T-GUARDRAIL: in a single Subpath / Subpath_1 in a PathBuilt, an
    # is_continuous=true first anchor has no prior phase to inherit. The
    # Subpath build leaves it pending; PathBuilt build throws.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0,
              meta = [Spinning(; rate = 1.0, is_continuous = true)])
    _seal_at_z(sb, 1.0)
    b = build(sb)
    @test b.pending_continuous_first_spinning
    @test_throws ArgumentError build([b])
end

@testset "Spinning — validation: two Spinnings per segment rejected" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0,
              meta = [Spinning(; rate = 1.0), Spinning(; rate = 2.0)])
    _seal_at_z(sb, 1.0)
    @test_throws ArgumentError build(sb)
end

@testset "Spinning — validation: phi_0 with is_continuous rejected at construction" begin
    @test_throws ArgumentError Spinning(; rate = 1.0, phi_0 = 0.7, is_continuous = true)
end

@testset "Spinning — frame() returns spinning_rate" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 2.5)])
    _seal_at_z(sb, 1.0)
    b = build(sb)
    @test frame(b, 0.4).spinning_rate == 2.5
end

@testset "Spinning — total_frame_rotation = τ_geom + Ω_spin" begin
    sb = SubpathBuilder(); start!(sb)
    # straight segment has τ_geom = 0, so total_frame_rotation = ∫τ_spin ds.
    straight!(sb; length = 2.0, meta = [Spinning(; rate = 0.5)])
    _seal_at_z(sb, 2.0)
    b = build(sb)
    @test isapprox(total_frame_rotation(b; s_start = 0.0, s_end = 2.0), 1.0; atol = 1e-12)
end

@testset "Spinning — path_spinning_breakpoints includes run boundaries" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 1.0)])
    straight!(sb; length = 1.0)
    _seal_at_z(sb, 3.0)
    b = build(sb)
    bps = path_spinning_breakpoints(b)
    @test 0.0 in bps
    @test 1.0 in bps
    # Run extends to the Subpath's local s_end (interior + connector). The
    # connector is degenerate so this is ≈ 3.0 to ~1e-6 m, but we assert the
    # exact endpoint that the resolver used.
    @test arc_length(b) ∈ bps
end

@testset "Spinning — seal! lead-out carries a spinning run" begin
    # T-PHYSICS: a Spinning on the natural seal anchors a run that begins at the
    # connector (the lead-out) and extends to the Subpath end. Over the interior
    # there is no spinning; over the lead-out the rate is the constant τ, so the
    # integrated spinning equals τ·extra.
    τ = 2.0
    L_int = 1.0
    L_extra = 0.5
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L_int)
    seal!(sb; extra = L_extra, meta = [Spinning(; rate = τ)])
    b = build(sb)
    @test length(b.resolved_spinning) == 1
    @test isapprox(b.resolved_spinning[1].s_eff_start, L_int; atol = 1e-9)
    @test spinning_rate(b, 0.5 * L_int) == 0.0              # interior: no spinning
    @test spinning_rate(b, L_int + 0.5 * L_extra) == τ      # lead-out: spinning τ
    @test isapprox(total_spinning(b; s_start = L_int, s_end = L_int + L_extra),
                   τ * L_extra; atol = 1e-9)
end

@testset "Spinning — jumpto! connector carries a spinning run" begin
    # T-GUARDRAIL: a Spinning on the terminal jumpto! connector anchors a run
    # over the connector, handled exactly like one on an interior segment. The
    # run starts at the connector's arc-length offset and the rate is reported
    # throughout the connector region.
    τ = 1.25
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)
    # Bend the terminal connector toward an off-axis target so it has real length.
    jumpto!(sb; point = (0.3, 0.0, 1.4), incoming_tangent = (1.0, 0.0, 0.0),
            meta = [Spinning(; rate = τ)])
    b = build(sb)
    @test length(b.resolved_spinning) == 1
    s_conn = Float64(_qc_nominalize(b.jumpto_placed.s_offset_eff))
    L = arc_length(b)
    @test L > s_conn                                        # connector has length
    @test isapprox(b.resolved_spinning[1].s_eff_start, s_conn; atol = 1e-9)
    @test spinning_rate(b, 0.5) == 0.0                      # interior: no spinning
    @test spinning_rate(b, 0.5 * (s_conn + L)) == τ         # within the connector
end

@testset "Spinning — interior run ends where the seal's run begins" begin
    # T-GUARDRAIL: an interior Spinning and a seal Spinning produce two runs; the
    # interior run terminates at the connector's start (the next anchor), and the
    # seal run covers the lead-out.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0, meta = [Spinning(; rate = 1.0)])
    seal!(sb; extra = 0.5, meta = [Spinning(; rate = 3.0)])
    b = build(sb)
    @test length(b.resolved_spinning) == 2
    @test isapprox(b.resolved_spinning[1].s_eff_end, 1.0; atol = 1e-9)
    @test isapprox(b.resolved_spinning[2].s_eff_start, 1.0; atol = 1e-9)
    @test spinning_rate(b, 0.5) == 1.0                      # interior run
    @test spinning_rate(b, 1.25) == 3.0                     # seal run over lead-out
end

@testset "Spinning — is_continuous seal run inherits interior phase" begin
    # T-PHYSICS: a seal Spinning with is_continuous=true inherits phi_0 from the
    # interior run's accumulated phase, exactly as a segment-to-segment handoff.
    L1 = 1.5
    τ1 = 2.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L1, meta = [Spinning(; rate = τ1, phi_0 = 0.5)])
    seal!(sb; extra = 1.0, meta = [Spinning(; rate = 1.0, is_continuous = true)])
    b = build(sb)
    @test b.resolved_spinning[1].phi_0 == 0.5
    @test isapprox(b.resolved_spinning[2].phi_0, 0.5 + τ1 * L1; atol = 1e-12)
end

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

@testset "Path — total_turning_angle of a full circle" begin
    # T-PHYSICS: ∫κ ds over a full circle = 2π. After the bend we add a
    # tiny straight to neutralize the terminal curvature so the degenerate
    # terminal connector contributes ≈ 0 to the integrated curvature
    # (K0=K1=0, chord=0, matching tangent → κ = 0 along the connector
    # except at the inflection point, which has measure zero in the
    # integral). The bend itself contributes exactly 2π.
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = 0.1, angle = 2π)
    straight!(sb; length = 1e-6)   # cool off curvature
    # Natural exit: back at origin (~) with tangent (0,0,1).
    jumpto!(sb; point = (0.0, 0.0, 1e-6), incoming_tangent = (0.0, 0.0, 1.0))
    b = build(sb)
    @test total_turning_angle(b) ≈ 2π atol = 1e-3
end

@testset "Path — total_torsion of straight and bend segments is zero" begin
    # T-PHYSICS: straight and circular-arc segments have zero geometric torsion
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)
    bend!(sb; radius = 0.2, angle = π / 2)
    # Natural exit: position (0.2, 0, 1.2) with tangent (1, 0, 0).
    jumpto!(sb; point = (0.2, 0.0, 1.2), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)
    # Geometric torsion of straight + bend is zero. The terminal connector
    # may have nonzero torsion in general, but at chord ≈ 0 it's negligible
    # (the connector polynomial reduces to a near-point curve).
    @test abs(total_torsion(b)) < 1e-3
end

@testset "Path — writhe of a straight path is zero" begin
    # T-PHYSICS: a straight line has no self-linking, Wr = 0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 2.0)
    _seal_at_z(sb, 2.0)
    b = build(sb)
    @test abs(writhe(b; n = 64)) < 1e-6
end

# -----------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------

@testset "sample_uniform — returns n frames" begin
    # T-GUARDRAIL
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 1.0)
    bend!(sb; radius = 0.1, angle = π / 2)
    jumpto!(sb; point = (0.1, 0.0, 1.1), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)

    frames = sample_uniform(b; n = 50)
    @test length(frames) == 50
    @test hasproperty(frames[1], :position)
    @test hasproperty(frames[1], :tangent)
    @test hasproperty(frames[1], :curvature)
end

# -----------------------------------------------------------------------
# JumpBy interior segment
# -----------------------------------------------------------------------

@testset "JumpBy — endpoint matches delta in local frame" begin
    # T-PHYSICS: JumpBy with delta along the current tangent direction (ẑ)
    # should move the path forward by that amount.
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.0, 0.0, 0.5))
    jumpto!(sb; point = (0.0, 0.0, 0.5), incoming_tangent = (0.0, 0.0, 1.0))
    b = build(sb)

    # The end of the JumpBy interior segment sits at the start of the
    # terminal connector (s = jumpto_placed.s_offset_eff).
    @test position(b, 0.0) ≈ [0.0, 0.0, 0.0] atol = 1e-10
    @test position(b, _s_end_interior(b)) ≈ [0.0, 0.0, 0.5] atol = 1e-10
end

@testset "JumpBy — initial tangent is ẑ" begin
    # T-GUARDRAIL: incoming tangent at start of every connector is ẑ
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.1, 0.0, 0.3))
    jumpto!(sb; point = (0.1, 0.0, 0.3))   # default tangent
    b = build(sb)
    @test tangent(b, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "JumpBy — outgoing tangent honours tangent_out" begin
    # T-PHYSICS: explicit tangent_out should be the end tangent (in local frame)
    t_out = normalize([1.0, 0.0, 1.0])
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.3, 0.0, 0.3), tangent = t_out)
    # Terminate with the JumpBy's natural endpoint and matching tangent.
    jumpto!(sb; point = (0.3, 0.0, 0.3), incoming_tangent = t_out)
    b = build(sb)
    @test tangent(b, _s_end_interior(b)) ≈ t_out atol = 1e-8
end

@testset "JumpBy — frame orthonormality along connector" begin
    # T-GUARDRAIL
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.2, 0.1, 0.4))
    jumpto!(sb; point = (0.2, 0.1, 0.4))   # default tangent (chord direction)
    b = build(sb)
    for s in range(0.0, _s_end_interior(b); length = 11)
        T = tangent(b, s)
        N = normal(b, s)
        B = binormal(b, s)
        @test is_orthonormal(T, N, B)
    end
end

@testset "JumpBy — after straight, delta is in rotated local frame" begin
    # T-PHYSICS: after a quarter-circle bend the local ẑ points along global +x.
    # JumpBy with delta=(0,0,d) should therefore move +x in global frame.
    R = 0.1; d = 0.5
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpby!(sb; delta = (0.0, 0.0, d))   # local ẑ is now global +x
    # End of jumpby: (R + d, 0, R) with tangent (1, 0, 0).
    jumpto!(sb; point = (R + d, 0.0, R), incoming_tangent = (1.0, 0.0, 0.0))
    b = build(sb)
    @test position(b, _s_end_interior(b)) ≈ [R + d, 0.0, R] atol = 1e-8
end

# -----------------------------------------------------------------------
# Terminal jumpto behavior
# -----------------------------------------------------------------------

@testset "jumpto — endpoint matches point" begin
    # T-PHYSICS: the terminal connector lands the path at the specified
    # global point.
    dest = (1.0, 0.5, 2.0)
    sb = SubpathBuilder(); start!(sb)
    jumpto!(sb; point = dest)
    b = build(sb)
    # end_point of the SubpathBuilt is at the connector's end = dest.
    @test end_point(b) ≈ collect(dest) atol = 1e-10
end

@testset "jumpto — initial tangent is ẑ" begin
    sb = SubpathBuilder(); start!(sb)
    jumpto!(sb; point = (0.3, 0.1, 0.8))
    b = build(sb)
    @test tangent(b, 0.0) ≈ [0.0, 0.0, 1.0] atol = 1e-10
end

@testset "jumpto — after bend, point is in global frame" begin
    # T-PHYSICS: jumpto.point is always in global frame regardless of prior segments.
    R = 0.1
    dest = (R + 0.5, 0.0, R)
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = R, angle = π / 2, axis_angle = 0.0)
    jumpto!(sb; point = dest)
    b = build(sb)
    @test end_point(b) ≈ collect(dest) atol = 1e-8
end

@testset "jumpto — incoming_tangent honoured (global frame)" begin
    # T-PHYSICS: incoming_tangent for jumpto is specified in global frame
    t_in_global = normalize([1.0, 0.0, 0.0])
    sb = SubpathBuilder(); start!(sb)
    jumpto!(sb; point = (1.0, 0.0, 0.5), incoming_tangent = t_in_global)
    b = build(sb)
    @test tangent(b, Float64(_qc_nominalize(arc_length(b)))) ≈ t_in_global atol = 1e-8
end

@testset "JumpBy / jumpto — sample_path works on connectors" begin
    # T-GUARDRAIL: sample_path must not error on paths containing connectors
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.3)
    jumpby!(sb; delta = (0.1, 0.0, 0.4))
    straight!(sb; length = 0.2)
    # End of last straight: (0.1, 0, 0.9) with tangent (0,0,1).
    jumpto!(sb; point = (0.1, 0.0, 0.9), incoming_tangent = (0.0, 0.0, 1.0))
    b = build(sb)
    s_lo = 0.0
    s_hi = Float64(_qc_nominalize(arc_length(b)))
    ps = sample_path(b, s_lo, s_hi)
    @test ps.n >= 4
    @test ps.samples[1].s   ≈ s_lo atol = 1e-12
    @test ps.samples[end].s ≈ s_hi atol = 1e-12
end

# -----------------------------------------------------------------------
# PathBuilt assembly
# -----------------------------------------------------------------------

@testset "PathBuilt — builds from Vector{Subpath}" begin
    sub1 = let sb = SubpathBuilder()
        start!(sb)
        straight!(sb; length = 1.0)
        jumpto!(sb; point = (0.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, 1.0))
        Subpath(sb)
    end
    sub2 = let sb = SubpathBuilder()
        start!(sb; point = (0.0, 0.0, 1.0), outgoing_tangent = (0.0, 0.0, 1.0))
        straight!(sb; length = 0.5)
        jumpto!(sb; point = (0.0, 0.0, 1.5), incoming_tangent = (0.0, 0.0, 1.0))
        Subpath(sb)
    end
    p = build([sub1, sub2])
    @test length(p.subpaths) == 2
    @test isapprox(s_end(p), 1.5; atol = 1e-3)
    @test position(p, 0.0)       ≈ [0.0, 0.0, 0.0] atol = 1e-8
    @test isapprox(end_point(p), [0.0, 0.0, 1.5]; atol = 1e-3)
end

@testset "PathBuilt — Vector{SubpathBuilt} stitching matches Vector{Subpath}" begin
    # T-GUARDRAIL: build(::Vector{SubpathBuilt}) and build(::Vector{Subpath})
    # produce equivalent PathBuilts.
    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = 1.0)
    jumpto!(sb1; point = (0.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, 1.0))
    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 1.0), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.7)
    jumpto!(sb2; point = (0.0, 0.0, 1.7), incoming_tangent = (0.0, 0.0, 1.0))

    sub1 = Subpath(sb1); sub2 = Subpath(sb2)
    pa = build([sub1, sub2])
    pb = build([build(sub1), build(sub2)])
    @test length(pa.subpaths) == length(pb.subpaths) == 2
    @test isapprox(s_end(pa), s_end(pb); atol = 1e-12)
    @test position(pa, 0.5) ≈ position(pb, 0.5) atol = 1e-10
end

@testset "PathBuilt — single SubpathBuilt convenience" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.4)
    _seal_at_z(sb, 0.4)
    p = build(build(sb))   # build(SubpathBuilt) → PathBuilt
    @test p isa PathBuilt
    @test length(p.subpaths) == 1
end

@testset "PathBuilt — conformity check rejects mismatched start_point" begin
    # T-GUARDRAIL: Subpath_2 start_point must equal Subpath_1 jumpto_point.
    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = 1.0)
    jumpto!(sb1; point = (0.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, 1.0))

    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 2.0))   # mismatch
    straight!(sb2; length = 0.5)
    jumpto!(sb2; point = (0.0, 0.0, 2.5), incoming_tangent = (0.0, 0.0, 1.0))

    @test_throws ArgumentError build([Subpath(sb1), Subpath(sb2)])
end

@testset "PathBuilt — conformity check rejects mismatched tangent" begin
    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = 1.0)
    jumpto!(sb1; point = (0.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, 1.0))

    sb2 = SubpathBuilder()
    # Mismatched outgoing tangent: prev declared (0,0,1) coming in.
    start!(sb2; point = (0.0, 0.0, 1.0), outgoing_tangent = (1.0, 0.0, 0.0))
    straight!(sb2; length = 0.5)
    jumpto!(sb2; point = (0.5, 0.0, 1.0), incoming_tangent = (1.0, 0.0, 0.0))

    @test_throws ArgumentError build([Subpath(sb1), Subpath(sb2)])
end

@testset "PathBuilt — Spinning(is_continuous) inherits from prior Subpath" begin
    # T-PHYSICS: the second Subpath's is_continuous=true first anchor inherits
    # phi_0 from the prior Subpath's terminal spinning phase.
    L1 = 2.0; rate1 = 0.5; phi0_1 = 0.3
    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = L1, meta = [Spinning(; rate = rate1, phi_0 = phi0_1)])
    jumpto!(sb1; point = (0.0, 0.0, L1), incoming_tangent = (0.0, 0.0, 1.0))

    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, L1), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 1.0,
              meta = [Spinning(; rate = 1.0, is_continuous = true)])
    jumpto!(sb2; point = (0.0, 0.0, L1 + 1.0), incoming_tangent = (0.0, 0.0, 1.0))

    p = build([Subpath(sb1), Subpath(sb2)])
    # Subpath_1's terminal phase = phi0_1 + rate1 * (run length up to s_end of subpath1).
    # The run extends from s=0 to the Subpath_1's local s_end.
    sub1_s_end = arc_length(p.subpaths[1])
    expected_phi0 = phi0_1 + rate1 * Float64(_qc_nominalize(sub1_s_end))
    @test isapprox(p.subpaths[2].resolved_spinning[1].phi_0, expected_phi0;
                   atol = 1e-8)
end

# -----------------------------------------------------------------------
# Per-segment meta bag
# -----------------------------------------------------------------------

@testset "per-segment meta — builders forward meta to every segment type" begin
    nick = [Nickname("alpha")]
    mcm  = [MCMadd(:T_K, (:Normal, 0.0, 1.0))]

    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.1, meta = nick)
    bend!(sb; radius = 0.05, angle = π / 2, meta = mcm)
    helix!(sb; radius = 0.02, pitch = 0.01, turns = 1.0,
           meta = [Nickname("helix"), MCMadd(:T_K, :stub)])
    catenary!(sb; a = 0.04, length = 0.05, meta = [Nickname("cat")])
    jumpby!(sb; delta = (0.0, 0.0, 0.05), meta = [Nickname("jb")])
    jumpto!(sb; point = (0.0, 0.1, 0.4))

    segs = sb.segments
    @test segs[1].meta == nick
    @test segs[2].meta == mcm
    @test length(segs[3].meta) == 2
    @test segs[4].meta[1] isa Nickname
    @test segs[5].meta[1] isa Nickname

    @test segment_meta(segs[1]) === segs[1].meta
    @test segment_nickname(segs[1]) == "alpha"
    @test isnothing(segment_nickname(segs[2]))  # mcm only
end

@testset "per-segment meta — build() copies jumpby meta onto QuinticConnector" begin
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.1)
    jumpby!(sb; delta = (0.05, 0.0, 0.2),
            meta = [Nickname("connector-1"),
                    MCMadd(:T_K, :stub)])
    jumpto!(sb; point = (0.05, 0.0, 0.3))   # placeholder
    b = build(sb)

    hc = b.placed_segments[2].segment
    @test hc isa QuinticConnector
    @test length(hc.meta) == 2
    @test segment_nickname(hc) == "connector-1"
    @test any(m -> m isa MCMadd, hc.meta)
end

@testset "per-segment meta — segment_meta returns empty default" begin
    seg = StraightSegment(0.1)
    @test segment_meta(seg) == AbstractMeta[]
    @test isnothing(segment_nickname(seg))
end

# -----------------------------------------------------------------------
# G2 integration: curvature_out propagation
# -----------------------------------------------------------------------

@testset "JumpBy — curvature_out matches sampled κ at end of connector" begin
    # T-PHYSICS: G2 outgoing match. The connector's sampled scalar curvature
    # at its endpoint must equal ‖curvature_out‖.
    sb = SubpathBuilder(); start!(sb)
    K1 = (0.0, 2.0, 0.0)   # 2 m⁻¹ in local +y
    jumpby!(sb; delta = (0.5, 0.0, 0.5),
            tangent = (1.0, 0.0, 0.0),
            curvature_out = K1)
    jumpto!(sb; point = (0.5, 0.0, 0.5), incoming_tangent = (1.0, 0.0, 0.0),
            incoming_curvature = (0.0, 2.0, 0.0))
    b = build(sb)
    seg = b.placed_segments[1].segment
    L = arc_length(seg)
    @test isapprox(curvature(seg, L), sqrt(K1[1]^2 + K1[2]^2 + K1[3]^2);
                   rtol = 1e-3, atol = 1e-6)
end

@testset "JumpBy — incoming K0 inherited from prior bend (G2 join)" begin
    # T-PHYSICS: G2 incoming match. After a bend the connector's start curvature
    # must equal 1/R_bend.
    R_bend = 0.5
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = R_bend, angle = π/4)
    jumpby!(sb; delta = (0.3, 0.0, 0.3))
    jumpto!(sb; point = (0.3 + R_bend*sin(π/4), 0.0, R_bend*(1-cos(π/4)) + 0.3),
            incoming_tangent = (sin(π/4) + 0.3 / 0.3, 0.0, cos(π/4)))   # placeholder; ok if connector has small length
    b = build(sb)
    seg = b.placed_segments[2].segment
    @test isapprox(curvature(seg, 0.0), 1.0 / R_bend; rtol = 1e-2, atol = 1e-4)
end

@testset "jumpto — global-frame incoming_curvature after bend" begin
    # T-PHYSICS: jumpto with incoming_curvature specified in the *global* frame,
    # placed after a BendSegment that has rotated the local frame.
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = 0.3, angle = π/2)   # rotates local frame
    K1_global = (0.0, 1.0, 0.0)
    jumpto!(sb; point = (0.6, 0.0, 0.6),
            incoming_tangent = (0.0, 0.0, 1.0),
            incoming_curvature = K1_global)
    b = build(sb)
    # The terminal connector is stored on jumpto_quintic_connector / jumpto_placed.
    seg = b.jumpto_quintic_connector
    L = arc_length(seg)
    @test isapprox(curvature(seg, L), 1.0; rtol = 1e-3, atol = 1e-6)
end

# -----------------------------------------------------------------------
# start_outgoing_curvature non-default (T-PHYSICS)
# -----------------------------------------------------------------------

@testset "start_outgoing_curvature — non-zero curvature flows into first JumpBy" begin
    # T-PHYSICS: a Subpath with non-zero start_outgoing_curvature followed by
    # a JumpBy interior segment yields a connector whose incoming κ matches
    # the start curvature (G2 join at the Subpath start).
    κ_start = 1.5
    sb = SubpathBuilder()
    start!(sb; point = (0.0, 0.0, 0.0),
              outgoing_tangent = (0.0, 0.0, 1.0),
              outgoing_curvature = (κ_start, 0.0, 0.0))   # along +x normal
    jumpby!(sb; delta = (0.2, 0.0, 0.3))
    jumpto!(sb; point = (0.2, 0.0, 0.3))
    b = build(sb)
    seg = b.placed_segments[1].segment
    @test seg isa QuinticConnector
    @test isapprox(curvature(seg, 0.0), κ_start; rtol = 1e-2, atol = 1e-4)
end
