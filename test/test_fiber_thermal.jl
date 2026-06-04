using Test
using LinearAlgebra
using Bifrost
using Bifrost.PathGeometry: _qc_nominalize
using MonteCarloMeasurements

# Fiber-layer thermal interpretation: Fiber(builder/Subpath; …) bakes :T_K into
# an isotropic length scaling (1 + α_lin·ΔT) using the cladding CTE, then builds
# once. α_lin is consulted lazily — only when a :T_K segment is present.

const _FT_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),     # pure silica cladding → α_lin = SILICA_CTE
    8.2e-6,
    125e-6,
)
const _FT_T_REF = 297.15
const _FT_ALPHA = cte(_FT_XS.cladding_material, _FT_T_REF)

# ΔT that produces a given total-length factor α = 1 + α_lin·ΔT.
_ft_ΔT_for(α) = (α - 1) / _FT_ALPHA
_ft_mcm(α)    = [MCMadd(:T_K, _ft_ΔT_for(α))]

function _ft_subpath(f::Function; spin_rate = nothing)
    sb = SubpathBuilder(); start!(sb; spin_rate = spin_rate)
    f(sb)
    isnothing(sb.jumpto_point) && seal!(sb)
    return Subpath(sb)
end

_ft_baseline(f; spin_rate = nothing) = build(_ft_subpath(f; spin_rate = spin_rate))
_ft_scaled(f; spin_rate = nothing)   =
    Fiber(_ft_subpath(f; spin_rate = spin_rate);
          cross_section = _FT_XS, T_ref_K = _FT_T_REF).path

# -----------------------------------------------------------------------
# Uniform thermal scaling via :T_K
# -----------------------------------------------------------------------

@testset "Fiber :T_K — uniform scales interior arc length by α" begin
    # T-PHYSICS: uniform α on every segment multiplies total interior arc length
    # by α.
    α = 1.05
    spec = sb -> begin
        straight!(sb; length = 1.0, meta = _ft_mcm(α))
        bend!(sb; radius = 0.1, angle = π / 2, meta = _ft_mcm(α))
    end
    ib = sum(arc_length(ps.segment) for ps in _ft_baseline(spec).placed_segments)
    is = sum(arc_length(ps.segment) for ps in _ft_scaled(spec).placed_segments)
    @test is ≈ α * ib atol = 1e-9
end

@testset "Fiber :T_K — preserves joint tangent continuity" begin
    # T-GUARDRAIL
    spec = sb -> begin
        straight!(sb; length = 0.5, meta = _ft_mcm(0.98))
        bend!(sb; radius = 0.1, angle = π / 3, meta = _ft_mcm(0.98))
    end
    path = _ft_scaled(spec)
    ps = path.placed_segments
    for i in 1:(length(ps) - 1)
        s_joint = Float64(_qc_nominalize(ps[i + 1].s_offset_eff))
        T_before = tangent(path, s_joint - 1e-9)
        T_after  = tangent(path, s_joint + 1e-9)
        @test norm(T_before - T_after) < 1e-6
    end
end

@testset "Fiber :T_K — BendSegment preserves angle, scales radius" begin
    # T-PHYSICS: α scales R_eff = α·R; swept angle preserved.
    spec = sb -> bend!(sb; radius = 0.1, angle = π / 3, axis_angle = 0.0, meta = _ft_mcm(1.1))
    seg1 = _ft_baseline(spec).placed_segments[1].segment
    seg2 = _ft_scaled(spec).placed_segments[1].segment
    @test seg2.radius ≈ 1.1 * seg1.radius
    @test seg2.angle == seg1.angle
    @test arc_length(seg2) ≈ 1.1 * arc_length(seg1)
end

@testset "Fiber :T_K — CatenarySegment scales arc length and parameter a" begin
    # T-PHYSICS
    spec = sb -> catenary!(sb; a = 0.2, length = 1.0, meta = _ft_mcm(0.95))
    seg  = _ft_baseline(spec).placed_segments[1].segment
    segs = _ft_scaled(spec).placed_segments[1].segment
    @test segs.a ≈ 0.95 * seg.a
    @test segs.length ≈ 0.95 * seg.length
    @test arc_length(segs) ≈ 0.95 * arc_length(seg)
end

@testset "Fiber :T_K — HelixSegment scales arc length but preserves turns" begin
    spec = sb -> helix!(sb; radius = 0.03, pitch = 0.01, turns = 2.0, meta = _ft_mcm(0.9))
    seg  = _ft_baseline(spec).placed_segments[1].segment
    segs = _ft_scaled(spec).placed_segments[1].segment
    @test segs.turns == seg.turns
    @test segs.radius ≈ 0.9 * seg.radius
    @test segs.pitch ≈ 0.9 * seg.pitch
    @test arc_length(segs) ≈ 0.9 * arc_length(seg)
end

# -----------------------------------------------------------------------
# Per-segment independence
# -----------------------------------------------------------------------

@testset "Fiber :T_K — per-segment scaling applies independently" begin
    path = _ft_scaled() do sb
        straight!(sb; length = 1.0, meta = _ft_mcm(1.0))   # α = 1 via ΔT = 0
        straight!(sb; length = 1.0, meta = _ft_mcm(0.5))
    end
    @test arc_length(path.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path.placed_segments[2].segment) ≈ 0.5
end

@testset "Fiber :T_K — segments without MCM annotations default to α = 1.0" begin
    path = _ft_scaled() do sb
        straight!(sb; length = 1.0)                     # no meta → α = 1
        straight!(sb; length = 2.0, meta = _ft_mcm(0.5))
    end
    @test arc_length(path.placed_segments[1].segment) ≈ 1.0
    @test arc_length(path.placed_segments[2].segment) ≈ 1.0
end

# -----------------------------------------------------------------------
# Composition: :T_K then field-level MCM
# -----------------------------------------------------------------------

@testset "Fiber :T_K — combined :T_K and direct :length on StraightSegment" begin
    # Thermal length scaling (via CTE) first, then the direct :length additive
    # offset on top.
    seg = _ft_scaled() do sb
        straight!(sb; length = 1.0, meta = [MCMadd(:T_K, 10.0), MCMadd(:length, 0.001)])
    end.placed_segments[1].segment
    expected = 1.0 * (1 + _FT_ALPHA * 10.0) + 0.001
    @test seg.length ≈ expected atol = 1e-14
end

# -----------------------------------------------------------------------
# MCM: Particles flow through via a Particles-valued ΔT
# -----------------------------------------------------------------------

@testset "Fiber :T_K — MCMadd Particles ΔT scales arc length" begin
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        # Particles-valued ΔT with zero mean → nominal interior arc length
        # preserved; output segment length lifts to Particles.
        ΔT = 0.0 ± (0.01 / _FT_ALPHA)    # σ_α = 0.01 around α = 1
        seg = _ft_scaled() do sb
            bend!(sb; radius = 0.05, angle = π / 2, meta = [MCMadd(:T_K, ΔT)])
        end.placed_segments[1].segment
        @test arc_length(seg) isa Particles
        @test pmean(arc_length(seg)) ≈ 0.05 * (π / 2) rtol = 1e-3
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end

# -----------------------------------------------------------------------
# Lazy α_lin: cte is consulted only when :T_K is present
# -----------------------------------------------------------------------

@testset "Fiber :T_K — cte-less cladding errors only when :T_K is present" begin
    # T-GUARDRAIL: a fluorine-doped cladding has no defined CTE. A non-thermal
    # fiber builds (cte never consulted); a :T_K fiber errors from cte.
    xs_f = StepIndexCrossSection(SilicaGermaniaGlass(0.036), SilicaFluorinatedGlass(0.01),
                             8.2e-6, 125e-6)
    plain = Fiber(_ft_subpath(sb -> straight!(sb; length = 1.0));
                  cross_section = xs_f, T_ref_K = _FT_T_REF)
    @test plain isa Fiber
    @test_throws ArgumentError Fiber(
        _ft_subpath(sb -> straight!(sb; length = 1.0, meta = [MCMadd(:T_K, 10.0)]));
        cross_section = xs_f, T_ref_K = _FT_T_REF)
end

# -----------------------------------------------------------------------
# Issue #33 — terminal jumpto! connector thermal expansion
# -----------------------------------------------------------------------

@testset "Fiber :T_K — jumpto! seal expands the connector by τ, endpoint fixed" begin
    # T-PHYSICS: a :T_K on the jumpto! seal scales the terminal connector's arc
    # length by τ (issue #33), still landing at the fixed jumpto_point.
    ΔT = 100.0
    P  = (0.1, 0.0, 0.5)
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.5)
    jumpto!(sb; point = P, incoming_tangent = (1.0, 0.0, 0.0), meta = [MCMadd(:T_K, ΔT)])
    L0 = Float64(_qc_nominalize(arc_length(build(sb).jumpto_quintic_connector)))

    f = Fiber(sb; cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    τ = 1 + _FT_ALPHA * ΔT
    @test isapprox(Float64(_qc_nominalize(arc_length(f.path.jumpto_quintic_connector))),
                   τ * L0; rtol = 1e-5)
    s_end = Float64(_qc_nominalize(arc_length(f.path)))
    @test isapprox(collect(position(f.path, s_end)), collect(P); atol = 1e-6)
end

@testset "Fiber Vector{SubpathBuilder} — routes through thermal jumpto! handling" begin
    # T-GUARDRAIL: Fiber([sb1, sb2]) must freeze builders to Subpaths and reuse
    # the Vector{Subpath} thermal path, NOT build([sb1, sb2]) first (which would
    # bypass per-subpath jumpto_target_length). A terminal jumpto! :T_K on the
    # last subpath proves the connector expansion (#33) still applies.
    ΔT = 100.0
    P1 = (0.0, 0.0, 0.5)
    P2 = (0.1, 0.0, 1.0)
    make_builders() = begin
        sb1 = SubpathBuilder(); start!(sb1)
        straight!(sb1; length = 0.5)
        jumpto!(sb1; point = P1, incoming_tangent = (0.0, 0.0, 1.0))

        sb2 = SubpathBuilder()
        start!(sb2; point = P1, outgoing_tangent = (0.0, 0.0, 1.0))
        straight!(sb2; length = 0.4)
        jumpto!(sb2; point = P2, incoming_tangent = (1.0, 0.0, 0.0),
                meta = [MCMadd(:T_K, ΔT)])
        (sb1, sb2)
    end

    sb1, sb2 = make_builders()
    f_conv = Fiber([sb1, sb2]; cross_section = _FT_XS, T_ref_K = _FT_T_REF)

    sb1e, sb2e = make_builders()
    f_ref = Fiber([Subpath(sb1e), Subpath(sb2e)];
                  cross_section = _FT_XS, T_ref_K = _FT_T_REF)

    @test length(f_conv.path.subpaths) == length(f_ref.path.subpaths) == 2

    s_conv = Float64(_qc_nominalize(arc_length(f_conv.path)))
    s_ref  = Float64(_qc_nominalize(arc_length(f_ref.path)))
    @test isapprox(s_conv, s_ref; atol = 1e-9)
    @test isapprox(collect(position(f_conv.path, s_conv)),
                   collect(position(f_ref.path, s_ref)); atol = 1e-9)

    # The terminal connector of the last subpath thermally expanded by τ — and
    # it is identical between the convenience and explicit forms.
    conn_conv = arc_length(f_conv.path.subpaths[end].jumpto_quintic_connector)
    conn_ref  = arc_length(f_ref.path.subpaths[end].jumpto_quintic_connector)
    @test isapprox(Float64(_qc_nominalize(conn_conv)),
                   Float64(_qc_nominalize(conn_ref)); rtol = 1e-9)
end

# -----------------------------------------------------------------------
# Seal meta validation: only MCMadd(:T_K, …) is supported on jumpto!
# -----------------------------------------------------------------------

@testset "Fiber — field MCM on a jumpto! seal errors; :T_K / Nickname are fine" begin
    # T-GUARDRAIL: the terminal connector supports only MCMadd(:T_K, …) (#33).
    # A field-level MCMadd/MCMmul (or a multiplicative :T_K) on the seal is
    # rejected at Fiber construction rather than silently ignored.
    mk(meta) = begin
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 0.5)
        jumpto!(sb; point = (0.1, 0.0, 0.5), incoming_tangent = (1.0, 0.0, 0.0), meta = meta)
        sb
    end
    @test_throws ArgumentError Fiber(mk([MCMadd(:length, 0.01)]);
                                     cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    @test_throws ArgumentError Fiber(mk([MCMmul(:T_K, 1.1)]);
                                     cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    # Supported: thermal :T_K, and a plain Nickname.
    @test Fiber(mk([MCMadd(:T_K, 10.0)]); cross_section = _FT_XS, T_ref_K = _FT_T_REF) isa Fiber
    @test Fiber(mk([Nickname("seal")]);  cross_section = _FT_XS, T_ref_K = _FT_T_REF) isa Fiber
end

# -----------------------------------------------------------------------
# min_bend_radius is still honored on jumpto! and jumpby!
# -----------------------------------------------------------------------

@testset "Fiber — jumpto! min_bend_radius honored (no :T_K)" begin
    # T-GUARDRAIL: with a generous radius limit the build succeeds and respects it.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.2)
    jumpto!(sb; point = (0.1, 0.0, 0.4), incoming_tangent = (0.0, 0.0, 1.0),
            min_bend_radius = 0.02)
    f = Fiber(sb; cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    conn = f.path.jumpto_quintic_connector
    L = arc_length(conn)
    κ_peak = maximum(curvature(conn, s) for s in range(0.0, Float64(_qc_nominalize(L)); length = 64))
    @test κ_peak <= 1 / 0.02 + 1e-6
end

@testset "Fiber — jumpto! min_bend_radius honored under :T_K (#33)" begin
    # T-GUARDRAIL: a seal carrying :T_K thermally expands the connector
    # (target-length solve) yet still respects a feasible min_bend_radius — the
    # solver validates peak curvature against the limit post-hoc.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.2)
    jumpto!(sb; point = (0.1, 0.0, 0.4), incoming_tangent = (0.0, 0.0, 1.0),
            min_bend_radius = 0.02, meta = [MCMadd(:T_K, 50.0)])
    f = Fiber(sb; cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    conn = f.path.jumpto_quintic_connector
    L = Float64(_qc_nominalize(arc_length(conn)))
    κ_peak = maximum(curvature(conn, s) for s in range(0.0, L; length = 64))
    @test κ_peak <= 1 / 0.02 + 1e-6
end

@testset "Fiber — jumpby! min_bend_radius honored" begin
    # T-GUARDRAIL: JumpBy's min_bend_radius is honored at placement-time
    # resolution (unaffected by the thermal/perturb path).
    sb = SubpathBuilder(); start!(sb)
    jumpby!(sb; delta = (0.1, 0.0, 0.3), tangent = (0.0, 0.0, 1.0), min_bend_radius = 0.02)
    seal!(sb)
    f = Fiber(sb; cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    conn = f.path.placed_segments[1].segment
    L = arc_length(conn)
    κ_peak = maximum(curvature(conn, s) for s in range(0.0, Float64(_qc_nominalize(L)); length = 64))
    @test κ_peak <= 1 / 0.02 + 1e-6
end

# -----------------------------------------------------------------------
# Material spin under thermal expansion
# -----------------------------------------------------------------------

@testset "Fiber :T_K — preserves constant spin rate; total scales with length" begin
    # T-PHYSICS: thermal expansion scales arc length by α but leaves the spin
    # rate τ unchanged, so the integrated spin ∫τ ds scales by α.
    α = 1.05
    τ = 1.5
    spec = sb -> straight!(sb; length = 2.0, meta = _ft_mcm(α))

    base = _ft_baseline(spec; spin_rate = τ)
    scal = _ft_scaled(spec; spin_rate = τ)

    @test base.spin_rate == τ
    @test scal.spin_rate == τ          # rate preserved, not scaled

    Lb = Float64(_qc_nominalize(arc_length(base)))
    Ls = Float64(_qc_nominalize(arc_length(scal)))
    @test Ls ≈ α * Lb atol = 1e-9

    Ωb = total_spin(base; s_start = 0.0, s_end = Lb)
    Ωs = total_spin(scal; s_start = 0.0, s_end = Ls)
    @test Ωb ≈ τ * Lb atol = 1e-9
    @test Ωs ≈ α * Ωb rtol = 1e-9
end
