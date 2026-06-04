using Bifrost
using Bifrost.PathGeometry: _qc_nominalize
# Fiber assembly tests for `fiber-path.jl`: low-level binding of a built path,
# and the builder-accepting constructors that apply perturbation meta during the
# single build (thermal `:T_K` via the cladding CTE, field-level MCM, and the
# terminal-connector thermal expansion of issue #33).

using Test
using LinearAlgebra


const _XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036), SilicaGermaniaGlass(0.0),
    8.2e-6, 125e-6,
)
const _T_REF = 297.15
const _α_LIN = cte(_XS.cladding_material, _T_REF)

# ----------------------------
# Helpers
# ----------------------------

# Single-Subpath builder with a straight-line interior; seal at (0,0,L).
function _trivial_subpath(L::Float64; meta = AbstractMeta[])
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = meta)
    jumpto!(sb; point = (0.0, 0.0, L), incoming_tangent = (0.0, 0.0, 1.0))
    return Subpath(sb)
end

# ----------------------------
# 3a. Fiber construction (low-level: bind an already-built path as-is)
# ----------------------------

@testset "Fiber — construction from SubpathBuilt" begin
    # T-GUARDRAIL
    b = build(_trivial_subpath(2.0))
    f = Fiber(b; cross_section = _XS, T_ref_K = _T_REF)
    @test f.s_start == 0.0
    @test isapprox(f.s_end, arc_length(b); atol = 1e-12)
    @test f.path === b
end

@testset "Fiber — construction from PathBuilt" begin
    # T-GUARDRAIL
    sub1 = _trivial_subpath(1.0)
    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 1.0), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.5)
    jumpto!(sb2; point = (0.0, 0.0, 1.5), incoming_tangent = (0.0, 0.0, 1.0))
    p = build([sub1, Subpath(sb2)])
    f = Fiber(p; cross_section = _XS, T_ref_K = _T_REF)
    @test f.s_start == 0.0
    @test isapprox(f.s_end, arc_length(p); atol = 1e-12)
    @test f.path === p
end

@testset "Fiber — generator closures are callable" begin
    # T-GUARDRAIL
    b = build(_trivial_subpath(0.5))
    f = Fiber(b; cross_section = _XS, T_ref_K = _T_REF)
    K = generator_K(f, 1550e-9)
    Kω = generator_Kω(f, 1550e-9)
    M  = K(0.25)
    Mω = Kω(0.25)
    @test size(M) == (2, 2)
    @test size(Mω) == (2, 2)
    @test eltype(M) == ComplexF64
end

# ----------------------------
# 3b. Fiber(builder) — thermal :T_K scaling on a straight segment
# ----------------------------

@testset "Fiber — :T_K scales StraightSegment length by (1 + α·ΔT)" begin
    # T-PHYSICS
    L = 0.5
    ΔT = 100.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = [MCMadd(:T_K, ΔT)])
    jumpto!(sb; point = (0.0, 0.0, L), incoming_tangent = (0.0, 0.0, 1.0))
    f_nominal = Fiber(build(sb); cross_section = _XS, T_ref_K = _T_REF)  # bind as-is
    f         = Fiber(sb; cross_section = _XS, T_ref_K = _T_REF)         # thermal applied
    expected = L * (1 + _α_LIN * ΔT)
    @test isapprox(f.path.placed_segments[1].segment.length, expected; atol = 1e-12)
    # Binding the built path as-is leaves :T_K inert (length unchanged).
    @test f_nominal.path.placed_segments[1].segment.length == L
end

# ----------------------------
# 3c. Fiber(builder) — bend radius MCMadd (field-level perturbation)
# ----------------------------

@testset "Fiber — MCMadd(:radius) on BendSegment shifts radius only" begin
    # T-PHYSICS
    R0   = 0.10
    Δr   = 0.02
    θ    = π / 3
    sb = SubpathBuilder(); start!(sb)
    bend!(sb; radius = R0, angle = θ, axis_angle = 0.0,
          meta = [MCMadd(:radius, Δr)])
    # Pin endpoint at the natural exit of an *unperturbed* bend; the build
    # uses the same authored jumpto_point regardless.
    jumpto!(sb; point = (R0 * (1 - cos(θ)), 0.0, R0 * sin(θ)),
            incoming_tangent = (sin(θ), 0.0, cos(θ)))
    f = Fiber(sb; cross_section = _XS, T_ref_K = _T_REF)
    seg = f.path.placed_segments[1].segment
    @test isapprox(seg.radius, R0 + Δr; atol = 1e-12)
    @test seg.angle == θ            # untouched
    @test seg.axis_angle == 0.0      # untouched
end

# ----------------------------
# 3d. Fiber(builder) — issue #33: terminal connector thermal expansion
# ----------------------------

@testset "Fiber — :T_K on jumpto! expands the connector by τ, endpoint fixed" begin
    # T-PHYSICS (issue #33): a jumpto! seal carrying :T_K thermally expands the
    # terminal connector — its arc length scales by τ = 1 + α·ΔT — while still
    # landing at the fixed jumpto_point.
    L  = 0.5
    ΔT = 100.0
    P  = (0.1, 0.0, L)
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L)
    jumpto!(sb; point = P, incoming_tangent = (1.0, 0.0, 0.0), meta = [MCMadd(:T_K, ΔT)])
    L0 = Float64(_qc_nominalize(arc_length(build(sb).jumpto_quintic_connector)))

    f = Fiber(sb; cross_section = _XS, T_ref_K = _T_REF)
    τ = 1 + _α_LIN * ΔT
    L_conn = Float64(_qc_nominalize(arc_length(f.path.jumpto_quintic_connector)))
    @test isapprox(L_conn, τ * L0; rtol = 1e-5)

    # Endpoint preserved at the fixed jumpto_point.
    s_end = Float64(_qc_nominalize(arc_length(f.path)))
    @test isapprox(collect(position(f.path, s_end)), collect(P); atol = 1e-6)
end

@testset "Fiber — jumpto! without :T_K leaves the connector unconstrained" begin
    # T-GUARDRAIL: thermal connector expansion is opt-in via :T_K on the seal.
    # Without it, an expanding interior grows the total and the connector solves
    # naturally to the fixed endpoint.
    L  = 0.5
    ΔT = 100.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = [MCMadd(:T_K, ΔT)])   # interior expands
    jumpto!(sb; point = (0.1, 0.0, L), incoming_tangent = (1.0, 0.0, 0.0))  # no seal :T_K
    L_baseline = Float64(_qc_nominalize(arc_length(build(sb))))
    f = Fiber(sb; cross_section = _XS, T_ref_K = _T_REF)
    L_modified = Float64(_qc_nominalize(arc_length(f.path)))
    @test L_modified > L_baseline + 1e-9
    @test isapprox(f.path.placed_segments[1].segment.length, L * (1 + _α_LIN * ΔT);
                   atol = 1e-12)
end

# ----------------------------
# 3e. Fiber(builder) — PathBuilt across Subpaths
# ----------------------------

@testset "Fiber — PathBuilt builds each Subpath independently" begin
    # T-GUARDRAIL: a PathBuilt of two Subpaths with :T_K only on Subpath 2
    # leaves Subpath 1 unchanged.
    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = 0.4)
    jumpto!(sb1; point = (0.0, 0.0, 0.4), incoming_tangent = (0.0, 0.0, 1.0))

    sb2 = SubpathBuilder()
    start!(sb2; point = (0.0, 0.0, 0.4), outgoing_tangent = (0.0, 0.0, 1.0))
    straight!(sb2; length = 0.6, meta = [MCMadd(:T_K, 50.0)])
    jumpto!(sb2; point = (0.0, 0.0, 1.0), incoming_tangent = (0.0, 0.0, 1.0))

    subs = [Subpath(sb1), Subpath(sb2)]
    f_nominal = Fiber(build(subs); cross_section = _XS, T_ref_K = _T_REF)  # bind as-is
    f         = Fiber(subs; cross_section = _XS, T_ref_K = _T_REF)         # thermal applied
    @test f.path isa PathBuilt
    @test length(f.path.subpaths) == 2

    # Subpath 1 unchanged.
    seg1_orig = f_nominal.path.subpaths[1].placed_segments[1].segment
    seg1_mod  = f.path.subpaths[1].placed_segments[1].segment
    @test seg1_orig.length == seg1_mod.length

    # Subpath 2's straight expanded.
    expected = 0.6 * (1 + _α_LIN * 50.0)
    seg2_mod = f.path.subpaths[2].placed_segments[1].segment
    @test isapprox(seg2_mod.length, expected; atol = 1e-12)
end

# ----------------------------
# 3f. JumpBy interior segment — :T_K passes through (not modeled)
# ----------------------------

@testset "Fiber — :T_K on a JumpBy is inert (connector not scaled)" begin
    # T-GUARDRAIL: a JumpBy is resolved to a connector at placement time and
    # passes through the thermal transform unchanged, so :T_K on a JumpBy does
    # not scale its connector. (Thermal expansion of jump connectors is not
    # modeled; express such jumps with thermally-scaled segments instead.)
    ΔT = 200.0
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.0, 0.0, 0.4), meta = [MCMadd(:T_K, ΔT)])
    jumpto!(sb; point = (0.0, 0.0, 0.4), incoming_tangent = (0.0, 0.0, 1.0))
    L_baseline = Float64(_qc_nominalize(arc_length(build(sb).placed_segments[1].segment)))

    f = Fiber(sb; cross_section = _XS, T_ref_K = _T_REF)
    L_after = Float64(_qc_nominalize(arc_length(f.path.placed_segments[1].segment)))
    @test isapprox(L_after, L_baseline; rtol = 1e-12)
end
