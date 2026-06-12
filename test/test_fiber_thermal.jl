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
    spec = sb -> bend!(sb; radius = 0.1, angle = π / 3, axis_angle = 0.0,
                        meta = _ft_mcm(1.1))
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

@testset "Fiber :T_K — temperature derived on demand from segment meta" begin
    # T-GUARDRAIL: `temperature(f, s)` recovers each segment's `:T_K` excursion at
    # query time via `local_segment(f.path, s)` — no thermal state on `Fiber`.
    ΔT = 12.0
    sub = _ft_subpath() do sb
        straight!(sb; length = 1.0, meta = [MCMadd(:T_K, ΔT)])
        straight!(sb; length = 1.0)
    end
    f = Fiber(sub; cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    ps = f.path.placed_segments
    s1 = 0.5 * arc_length(ps[1].segment)
    s2 = Float64(_qc_nominalize(ps[2].s_offset_eff)) +
         0.5 * Float64(_qc_nominalize(arc_length(ps[2].segment)))

    @test temperature(f, s1) ≈ _FT_T_REF + ΔT
    @test temperature(f, s2) ≈ _FT_T_REF
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
# Terminal jumpto! connector thermal expansion
# -----------------------------------------------------------------------

@testset "Fiber :T_K — jumpto! seal expands the connector by τ, endpoint fixed" begin
    # T-PHYSICS: a :T_K on the jumpto! seal scales the terminal connector's arc
    # length by τ, still landing at the fixed jumpto_point.
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
    # last subpath proves the connector expansion still applies.
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
    # T-GUARDRAIL: the terminal connector supports only MCMadd(:T_K, …).
    # A field-level MCMadd/MCMmul (or a multiplicative :T_K) on the seal is
    # rejected at Fiber construction rather than silently ignored.
    mk(meta) = begin
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 0.5)
        jumpto!(sb; point = (0.1, 0.0, 0.5),
                incoming_tangent = (1.0, 0.0, 0.0), meta = meta)
        sb
    end
    @test_throws ArgumentError Fiber(mk([MCMadd(:length, 0.01)]);
                                     cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    @test_throws ArgumentError Fiber(mk([MCMmul(:T_K, 1.1)]);
                                     cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    # Supported: thermal :T_K, and a plain Nickname.
    @test Fiber(mk([MCMadd(:T_K, 10.0)]);
                cross_section = _FT_XS, T_ref_K = _FT_T_REF) isa Fiber
    @test Fiber(mk([Nickname("seal")]);
                cross_section = _FT_XS, T_ref_K = _FT_T_REF) isa Fiber
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
    ss = range(0.0, Float64(_qc_nominalize(L)); length = 64)
    κ_peak = maximum(curvature(conn, s) for s in ss)
    @test κ_peak <= 1 / 0.02 + 1e-6
end

@testset "Fiber — jumpto! min_bend_radius honored under :T_K" begin
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
    ss = range(0.0, Float64(_qc_nominalize(L)); length = 64)
    κ_peak = maximum(curvature(conn, s) for s in ss)
    @test κ_peak <= 1 / 0.02 + 1e-6
end

# -----------------------------------------------------------------------
# Material spin under thermal expansion
# -----------------------------------------------------------------------

@testset "Fiber :T_K — spin rotations per Subpath fixed regardless of temperature" begin
    # T-PHYSICS (issue #32): the number of material-spin rotations per Subpath
    # is fixed regardless of temperature. The rotations are baked into the
    # material at draw time and an isotropic thermal expansion by α only
    # stretches the existing pattern over a longer length — no rotations are
    # added or removed. As a consequence the spinning period (1/τ) is
    # temperature-dependent: τ divides by α so the period scales by α.
    # This mirrors how mechanical twist conserves total turns under thermal
    # expansion (see the twist tests below).
    α = 1.05
    τ = 1.5
    spec = sb -> straight!(sb; length = 2.0, meta = _ft_mcm(α))

    base = _ft_baseline(spec; spin_rate = τ)
    scal = _ft_scaled(spec; spin_rate = τ)

    @test base.spin_rate == τ
    @test scal.spin_rate ≈ τ / α rtol = 1e-12   # rate divides, total turns fixed

    Lb = Float64(_qc_nominalize(arc_length(base)))
    Ls = Float64(_qc_nominalize(arc_length(scal)))
    @test Ls ≈ α * Lb atol = 1e-9

    # Total spin rotations ∫τ ds are conserved across the thermal expansion.
    Ωb = total_spin(base; s_start = 0.0, s_end = Lb)
    Ωs = total_spin(scal; s_start = 0.0, s_end = Ls)
    @test Ωb ≈ τ * Lb atol = 1e-9
    @test Ωs ≈ Ωb rtol = 1e-9

    # The spinning period (length per rotation) is temperature-dependent.
    period_base = 2π / spin_rate(base, 0.5 * Lb)
    period_scal = 2π / spin_rate(scal, 0.5 * Ls)
    @test period_scal ≈ α * period_base rtol = 1e-12
end

@testset "Fiber :T_K — function-valued spin rate conserves turns (reparametrized)" begin
    # T-PHYSICS (issue #32, follow-on): a function rate τ(s_local) is
    # reparametrized onto the stretched local arc length, g(s) = τ(s/α)/α, so
    # ∫g over the elongated Subpath equals ∫τ over the original — turns
    # conserved for non-constant spin too.
    α    = 1.1
    rate = s -> 1.0 + s          # rad/m, linear in Subpath-local arc length
    spec = sb -> straight!(sb; length = 2.0, meta = _ft_mcm(α))

    base = _ft_baseline(spec; spin_rate = rate)
    scal = _ft_scaled(spec; spin_rate = rate)

    Lb = Float64(_qc_nominalize(arc_length(base)))
    Ls = Float64(_qc_nominalize(arc_length(scal)))
    @test total_spin(base; s_start = 0.0, s_end = Lb) ≈
          total_spin(scal; s_start = 0.0, s_end = Ls) rtol = 1e-9
end

# -----------------------------------------------------------------------
# Mechanical twist under thermal expansion (inverse-length scaling)
# -----------------------------------------------------------------------

@testset "Fiber :T_K — divides twist rate by α; conserves total turns" begin
    # T-PHYSICS: mechanical twist τ_m is an inverse-length rate (rad/m). Thermal
    # expansion scales arc length by α, so the *rate* divides by α and the total
    # accumulated twist Φ = ∫τ_m ds is conserved (the frozen-in turns just stretch
    # over a longer length — no turns are added). Contrast material spinning above,
    # whose rate is preserved and whose total therefore scales with length.
    α   = 1.05
    τm0 = 3.0
    spec = sb -> straight!(sb; length = 2.0, twist = τm0, meta = _ft_mcm(α))

    base = _ft_baseline(spec)
    scal = _ft_scaled(spec)

    Lb = Float64(_qc_nominalize(arc_length(base)))
    Ls = Float64(_qc_nominalize(arc_length(scal)))
    @test Ls ≈ α * Lb atol = 1e-9

    # Rate divides by α.
    @test twist_rate(base, 0.5 * Lb) ≈ τm0
    @test twist_rate(scal, 0.5 * Ls) ≈ τm0 / α rtol = 1e-12

    # Total twist Φ = ∫τ_m ds is conserved across the expansion.
    Φb = twist_phase(base, Lb)
    Φs = twist_phase(scal, Ls)
    @test Φb ≈ τm0 * Lb atol = 1e-9
    @test Φs ≈ Φb rtol = 1e-9
end

@testset "Fiber :T_K — function-valued twist conserves turns (reparametrized)" begin
    # T-PHYSICS: a function rate τ_m(s_local) is reparametrized onto the stretched
    # local arc length, g(s) = τ_m(s/α)/α, so ∫g over the elongated segment equals
    # ∫τ_m over the original — turns conserved for non-constant twist too.
    α    = 1.1
    rate = s -> 1.0 + s          # rad/m, linear in segment-local arc length
    spec = sb -> straight!(sb; length = 2.0, twist = rate, meta = _ft_mcm(α))

    base = _ft_baseline(spec)
    scal = _ft_scaled(spec)
    Lb = Float64(_qc_nominalize(arc_length(base)))
    Ls = Float64(_qc_nominalize(arc_length(scal)))

    @test twist_phase(base, Lb) ≈ twist_phase(scal, Ls) rtol = 1e-9
end

@testset "Fiber — non-thermal fiber leaves twist rate untouched" begin
    # T-GUARDRAIL: with no :T_K/:tension meta the scaling factor is never applied,
    # so the twist rate is bit-for-bit the authored value.
    τm0 = 2.5
    f = Fiber(_ft_subpath(sb -> straight!(sb; length = 1.0, twist = τm0));
              cross_section = _FT_XS, T_ref_K = _FT_T_REF)
    L = Float64(_qc_nominalize(arc_length(f.path)))
    @test twist_rate(f.path, 0.5 * L) == τm0
    @test twist_phase(f.path, L) ≈ τm0 * L
end

@testset "Fiber :T_K — twist rate divides under Particles ΔT" begin
    # MCM: a Particles-valued ΔT lifts the twist rate to Particles, dividing by the
    # Particles factor α = 1 + α_lin·ΔT (no ::Real slot, no coercion).
    MonteCarloMeasurements.unsafe_comparisons(true)
    try
        τm0 = 4.0
        ΔT  = 0.0 ± (0.01 / _FT_ALPHA)        # α = 1 ± 0.01
        f = Fiber(_ft_subpath(sb -> straight!(sb; length = 1.0, twist = τm0,
                                              meta = [MCMadd(:T_K, ΔT)]));
                  cross_section = _FT_XS, T_ref_K = _FT_T_REF)
        L = Float64(_qc_nominalize(arc_length(f.path)))
        r = twist_rate(f.path, 0.5 * L)
        @test r isa Particles
        @test pmean(r) ≈ τm0 rtol = 1e-3      # 1/α with zero-mean α centers on τm0
    finally
        MonteCarloMeasurements.unsafe_comparisons(false)
    end
end
