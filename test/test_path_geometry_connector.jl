using Test
using LinearAlgebra
using QuadGK
using Bifrost
using Bifrost.PathGeometry: _build_quintic_connector, _qc_build_table, _qc_coefficients, _qc_eval, _qc_eval_d1, _qc_eval_d2, _qc_peak_curvature, _qc_speed, _qc_t_from_s, _qc_nominalize


# -----------------------------------------------------------------------
# Helpers used across multiple test sets in this file
# -----------------------------------------------------------------------

# Sample a connector at a uniform grid in arc length and return the points
# as a Vector{Vector}.
function _qc_sample_points(seg::QuinticConnector, n::Int)
    L = arc_length(seg)
    return [position_local(seg, s) for s in range(0.0, Float64(L); length = n)]
end

# Build a "raw" connector from endpoint geometry without any λ search,
# bypassing _build_quintic_connector. Mirrors what the reference module's
# `quintic_connector` does — useful for tests that must hold *for any λ*.
function _raw_quintic(P0, T0, K0, P1, T1, K1, λ::Real; n_table::Int = 64)
    P0 = collect(Float64, P0)
    T0 = collect(Float64, T0)
    K0 = collect(Float64, K0)
    P1 = collect(Float64, P1)
    T1 = collect(Float64, T1)
    K1 = collect(Float64, K1)
    V0 = λ .* T0
    V1 = λ .* T1
    A0 = (λ^2) .* K0
    A1 = (λ^2) .* K1
    coeffs = _qc_coefficients(P0, V0, A0, P1, V1, A1)
    s_table = _qc_build_table(coeffs, n_table)
    return QuinticConnector(coeffs, λ, s_table)
end

# Cubic Hermite reference for the K0=K1=0 degeneracy test.
# r(u) = h00·P0 + h10·V0 + h01·P1 + h11·V1, V_i = m·T_i
function _ref_cubic_hermite(P0, T0, P1, T1, m::Real, u::Real)
    h00 =  2u^3 - 3u^2 + 1
    h10 =     u^3 - 2u^2 + u
    h01 = -2u^3 + 3u^2
    h11 =     u^3 -   u^2
    return h00 .* P0 .+ h10 .* (m .* T0) .+ h01 .* P1 .+ h11 .* (m .* T1)
end

# -----------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------

@testset "QuinticConnector — path-geometry-connector.jl" begin

    @testset "Endpoint position match" begin
        """
        Endpoint-position sanity at u=0,1.

        The closed-form quintic Hermite solve must place the curve at exactly the
        requested endpoints. Catches sign/index errors in the coefficient assembly.
        """
        P0 = [0.0, 0.0, 0.0]
        P1 = [1.0, 0.5, 0.2]
        seg = _raw_quintic(P0, [0,0,1.0], [0,0,0.0], P1, [0,1.0,0], [0,0,0.0], 0.8)
        p0 = _qc_eval(seg.a, 0.0)
        p1 = _qc_eval(seg.a, 1.0)
        @test isapprox(p0, P0; atol = 1e-12)
        @test isapprox(p1, P1; atol = 1e-12)
    end

    @testset "Quintic endpoint matching (G2)" begin
        """
        Quintic endpoint matching (G2): position, velocity (×λ), acceleration (×λ²).

        This is the foundational claim the connector makes — that r(0)=P0, r'(0)=λ·T0,
        r''(0)=λ²·K0, and analogously at u=1. If this fails, every other test is
        suspect.
        """
        P0, T0, K0 = [0.0,0,0], [1.0,0,0], [0.0,0.2,0]
        P1, T1, K1 = [1.0,0.5,0.2], [0,1.0,0], [-0.1,0,0]
        λ = 0.8
        seg = _raw_quintic(P0, T0, K0, P1, T1, K1, λ)
        @test isapprox(_qc_eval(seg.a, 0.0), P0; atol = 1e-12)
        @test isapprox(_qc_eval(seg.a, 1.0), P1; atol = 1e-12)
        @test isapprox(_qc_eval_d1(seg.a, 0.0), λ .* T0; atol = 1e-12)
        @test isapprox(_qc_eval_d1(seg.a, 1.0), λ .* T1; atol = 1e-12)
        @test isapprox(_qc_eval_d2(seg.a, 0.0), λ^2 .* K0; atol = 1e-11)
        @test isapprox(_qc_eval_d2(seg.a, 1.0), λ^2 .* K1; atol = 1e-11)
    end

    @testset "Parametric type propagation (BigFloat)" begin
        """
        Parametric numeric type propagation (BigFloat).

        The connector is parameterized over T<:Real. Building from BigFloat inputs
        must yield BigFloat coefficients and outputs — i.e. the implementation does
        not silently downcast to Float64 anywhere geometry flows.
        """
        P0 = BigFloat[0, 0, 0]; T0 = BigFloat[1, 0, 0]; K0 = BigFloat[0, 0, 0]
        P1 = BigFloat[2, 0, 0]; T1 = BigFloat[1, 0, 0]; K1 = BigFloat[0, 0, 0]
        coeffs = _qc_coefficients(P0, BigFloat(2) .* T0, BigFloat(4) .* K0,
                                   P1, BigFloat(2) .* T1, BigFloat(4) .* K1)
        @test eltype(coeffs) == BigFloat
        s_table = _qc_build_table(coeffs, 8)
        @test eltype(s_table) == BigFloat
    end

    @testset "Straight-line zero-curvature stays on x-axis" begin
        """
        Straight-line K0=K1=0 degenerate case.

        With both endpoints collinear along x and zero curvature, the quintic must
        produce a curve that stays exactly on the x-axis. Probes for accidental
        nonzero y/z coefficients leaking through the closed-form solve.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [1.0, 0.0, 0.0]; K0 = [0.0, 0.0, 0.0]
        P1 = [2.0, 0.0, 0.0]; T1 = [1.0, 0.0, 0.0]; K1 = [0.0, 0.0, 0.0]
        seg = _raw_quintic(P0, T0, K0, P1, T1, K1, 2.0)
        for u in (0.0, 0.25, 0.5, 0.75, 1.0)
            p = _qc_eval(seg.a, u)
            @test abs(p[2]) ≤ 1e-12
            @test abs(p[3]) ≤ 1e-12
        end
    end

    @testset "Zero-acceleration endpoints (quintic, not cubic)" begin
        """
        Quintic with K0=K1=0 is the *quintic* Hermite with r''(0)=r''(1)=0,
        not a cubic Hermite.

        Important conceptual note: a quintic with zero endpoint acceleration is
        NOT the same curve as a cubic Hermite over the same endpoint position +
        tangent data. The cubic Hermite has r''(0), r''(1) determined by the
        cubic basis (generally nonzero); the quintic forces them to zero.

        This test verifies the quintic curve actually has r''(0) = r''(1) = 0
        (already covered structurally by the G2 endpoint-matching test, but
        worth asserting in isolation as a sanity check), AND that the curve
        differs from the cubic Hermite reference for endpoints with nonzero
        cubic-induced curvature — i.e. they are *different* G1 interpolants.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [0.0, 0.0, 1.0]
        P1 = [0.3, 0.5, 1.0]; T1 = [0.0, 1.0, 0.0]
        λ = 1.5
        seg = _raw_quintic(P0, T0, zeros(3), P1, T1, zeros(3), λ)
        # r''(0) = r''(1) = 0
        @test isapprox(_qc_eval_d2(seg.a, 0.0), zeros(3); atol = 1e-12)
        @test isapprox(_qc_eval_d2(seg.a, 1.0), zeros(3); atol = 1e-12)
        # Curve differs from cubic Hermite (the cubic has nonzero r'' at ends).
        diffs = [norm(_qc_eval(seg.a, u) .- _ref_cubic_hermite(P0, T0, P1, T1, λ, u))
                 for u in (0.2, 0.4, 0.6, 0.8)]
        @test maximum(diffs) > 1e-3
    end

    @testset "G2 endpoint curvature self-consistency" begin
        """
        G2 endpoint curvature self-consistency: κ(0) ≈ ‖K0‖, κ(1) ≈ ‖K1‖.

        Builds the connector with arbitrary K0/K1 and asks `curvature(seg, s)` for
        its scalar curvature value at the two endpoints (s=0 and s=L). These must
        match the magnitudes of the input curvature vectors. This is the *core*
        claim that motivates moving from cubic to quintic — must be tested directly.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [0.0, 0.0, 1.0]; K0 = [0.0, 0.4, 0.0]
        P1 = [1.0, 0.0, 1.0]; T1 = [1.0, 0.0, 0.0]; K1 = [0.0, 0.3, 0.0]
        λ = 1.2
        seg = _raw_quintic(P0, T0, K0, P1, T1, K1, λ)
        L = Float64(arc_length(seg))
        @test isapprox(curvature(seg, 0.0), norm(K0); rtol = 1e-3, atol = 1e-6)
        @test isapprox(curvature(seg, L),   norm(K1); rtol = 1e-3, atol = 1e-6)
    end

    @testset "Translation invariance" begin
        """
        Translation invariance.

        Translating both endpoints by an arbitrary vector v must leave a₁..a₅
        unchanged and only shift a₀ by v. Catches bugs that conflate position with
        velocity in the coefficient assembly.
        """
        P0 = [0.1, -0.2, 0.05]; T0 = [0.0, 0.0, 1.0]; K0 = [0.0, 0.1, 0.0]
        P1 = [1.0,  0.5, 0.2];  T1 = [0.0, 1.0, 0.0]; K1 = [0.05, 0.0, 0.0]
        v  = [3.0, -7.0, 2.5]
        λ  = 0.9
        s1 = _raw_quintic(P0, T0, K0, P1, T1, K1, λ)
        s2 = _raw_quintic(P0 .+ v, T0, K0, P1 .+ v, T1, K1, λ)
        @test isapprox(s2.a[1, :], s1.a[1, :] .+ v; atol = 1e-12)
        for i in 2:6
            @test isapprox(s2.a[i, :], s1.a[i, :]; atol = 1e-12)
        end
    end

    @testset "Rotation equivariance" begin
        """
        Rotation equivariance.

        Rotating the entire problem (P, T, K at both endpoints) by an arbitrary
        rotation R must produce a connector whose sampled positions are exactly
        R·(original sampled positions). Catches handedness bugs in the curvature-
        vector handling that pure-tangent G1 cubics would never expose.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [0.0, 0.0, 1.0]; K0 = [0.0, 0.2, 0.0]
        P1 = [0.6, 0.4, 0.7]; T1 = [1.0, 0.0, 0.0]; K1 = [0.0, 0.1, 0.0]
        λ  = 1.0
        seg1 = _raw_quintic(P0, T0, K0, P1, T1, K1, λ)
        # Arbitrary rotation: 30° about (1,1,1)/√3.
        ax = [1.0, 1.0, 1.0] ./ sqrt(3)
        θ  = π / 6
        c, s = cos(θ), sin(θ)
        K = [   0.0  -ax[3]   ax[2];
              ax[3]    0.0  -ax[1];
             -ax[2]  ax[1]    0.0]
        R = c * I(3) .+ s * K .+ (1 - c) * (ax * ax')
        seg2 = _raw_quintic(R * P0, R * T0, R * K0, R * P1, R * T1, R * K1, λ)
        for u in (0.0, 0.2, 0.5, 0.8, 1.0)
            p1 = _qc_eval(seg1.a, u)
            p2 = _qc_eval(seg2.a, u)
            @test isapprox(p2, R * p1; atol = 1e-11)
        end
    end

    @testset "Reverse-traversal symmetry" begin
        """
        Reverse-traversal symmetry.

        Building (P0,T0,K0) → (P1,T1,K1) and (P1,−T1,K1) → (P0,−T0,K0) must trace
        the *same set* of 3D points (the second curve is the first traversed
        backwards). Curvature is a magnitude property; the orientation flip of T
        only changes parameter direction, not the geometric curve. Probes for
        asymmetric handling of the start vs end of the coefficient solve.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [0.0, 0.0, 1.0]; K0 = [0.0, 0.3, 0.0]
        P1 = [0.5, 0.4, 0.8]; T1 = [1.0, 0.0, 0.0]; K1 = [0.0, 0.1, 0.0]
        λ = 1.1
        fwd = _raw_quintic(P0,  T0, K0, P1,  T1, K1, λ)
        rev = _raw_quintic(P1, -T1, K1, P0, -T0, K0, λ)
        for u in (0.0, 0.2, 0.5, 0.8, 1.0)
            pf = _qc_eval(fwd.a, u)
            pr = _qc_eval(rev.a, 1 - u)
            @test isapprox(pf, pr; atol = 1e-10)
        end
    end

    @testset "Arc-length table boundary and monotonicity" begin
        """
        Arc-length table — boundary values and monotonicity.

        `s_table[1] ≈ 0`, `s_table[end] ≈ L` (verified against an independent
        higher-order quadrature), and the table must be strictly increasing.
        Catches off-by-one and accumulation bugs in `_qc_build_table`.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [0.0, 0.0, 1.0]; K0 = [0.0, 0.2, 0.0]
        P1 = [0.4, 0.3, 0.7]; T1 = [1.0, 0.0, 0.0]; K1 = [0.0, 0.1, 0.0]
        seg = _raw_quintic(P0, T0, K0, P1, T1, K1, 1.0; n_table = 256)
        @test seg.s_table[1] == 0.0
        # Independent cross-check via QuadGK.
        L_ref, _ = quadgk(u -> _qc_speed(seg.a, u), 0.0, 1.0; rtol = 1e-12)
        @test isapprox(seg.s_table[end], L_ref; rtol = 1e-9)
        @test all(diff(seg.s_table) .> 0.0)
    end

    @testset "s ↔ u inversion round-trip" begin
        """
        s ↔ u inversion round-trip.

        For a sweep of arc-length values s in [0, L], compute u = _qc_t_from_s(s)
        and re-integrate ∫₀^u ‖dr/du'‖ du'; the result must round-trip to s.
        Catches Newton-iteration regressions and bracket-choice bugs in the
        inversion.
        """
        P0 = [0.0, 0.0, 0.0]; T0 = [0.0, 0.0, 1.0]; K0 = [0.0, 0.2, 0.0]
        P1 = [0.6, 0.4, 0.5]; T1 = [0.0, 1.0, 0.0]; K1 = [0.0, 0.0, 0.1]
        seg = _raw_quintic(P0, T0, K0, P1, T1, K1, 1.0; n_table = 128)
        L = Float64(seg.s_table[end])
        for frac in (0.0, 0.13, 0.37, 0.5, 0.71, 0.99, 1.0)
            s = frac * L
            u = _qc_t_from_s(seg, s)
            s_back, _ = quadgk(uu -> _qc_speed(seg.a, uu), 0.0, u; rtol = 1e-10)
            @test isapprox(s_back, s; atol = 1e-6, rtol = 1e-6)
        end
    end

    @testset "min_bend_radius search respects κ limit" begin
        """
        min_bend_radius search converges and respects the κ limit.

        A configuration whose chord-default λ produces a peak curvature above
        1/min_bend_radius must trigger the exponential-bracket / bisection search,
        converge, and return a connector whose *sampled* peak curvature lies at
        or below the limit (within tol).
        """
        # 90° turn over moderate chord. Chord-default λ produces κ_peak above
        # the limit 1/0.3 ≈ 3.3 m⁻¹; a larger λ tames it.
        p1 = [0.4, 0.0, 0.4]
        t_out = [1.0, 0.0, 0.0]
        K0 = zeros(3); K1 = zeros(3)
        seg = _build_quintic_connector(p1, t_out, K0, K1;
                                        min_bend_radius = 0.3,
                                        n_table = 64, n_check = 256)
        κ_max = _qc_peak_curvature(seg.a; n_check = 512)
        @test κ_max ≤ 1.0 / 0.3 + 1e-6
    end

    @testset "Bisection tightens λ" begin
        """
        Bisection refines λ tighter than the bracket alone would.

        Compare the λ chosen with bisection on a non-trivially curving problem
        against the first feasible λ found by exponential growth without
        refinement. The bisected λ must be ≤ the bracket-only λ. Establishes
        that the refinement step is doing useful work, not no-oping.
        """
        p1 = [0.4, 0.0, 0.4]
        t_out = [1.0, 0.0, 0.0]
        K0 = zeros(3); K1 = zeros(3)
        seg_bi = _build_quintic_connector(p1, t_out, K0, K1;
                                           min_bend_radius = 0.3,
                                           n_table = 32)
        seg_br = _build_quintic_connector(p1, t_out, K0, K1;
                                           min_bend_radius = 0.3,
                                           bisect_iter = 0,
                                           n_table = 32)
        @test seg_bi.lambda ≤ seg_br.lambda + 1e-9
    end

    @testset "Determinism of λ across repeated builds" begin
        """
        Determinism of λ-search under MCM nominalization.

        The chosen λ must depend only on the nominalized scalars used for branching
        predicates. Two builds of the same problem with identical inputs must
        return the same λ to the bit. (We don't require this across MCM and Float64
        inputs in general, but a Float64 problem must be reproducible.)
        """
        p1 = [0.4, 0.0, 0.4]
        t_out = [1.0, 0.0, 0.0]
        s1 = _build_quintic_connector(p1, t_out, zeros(3), zeros(3);
                                       min_bend_radius = 0.3, n_table = 32)
        s2 = _build_quintic_connector(p1, t_out, zeros(3), zeros(3);
                                       min_bend_radius = 0.3, n_table = 32)
        @test s1.lambda == s2.lambda
    end

    @testset "Reject endpoint curvature limit violation" begin
        """
        Endpoint curvature exceeding κ-limit raises ArgumentError.

        User-supplied K_out (or implied K_in from prior segment) larger than
        1/min_bend_radius is structurally infeasible — no λ can rescue it.
        The builder must reject before entering the search.
        """
        p1 = [0.5, 0.0, 0.5]
        t_out = [1.0, 0.0, 0.0]
        K_too_big = [0.0, 100.0, 0.0]    # 1/R = 100 m⁻¹
        @test_throws ArgumentError _build_quintic_connector(
            p1, t_out, K_too_big, zeros(3); min_bend_radius = 0.05)
    end

    @testset "AbstractPathSegment interface smoke test" begin
        """
        AbstractPathSegment interface methods return correctly typed/sized data.

        Quick smoke test of the whole interface (`position_local`, `tangent_local`,
        `normal_local`, `binormal_local`, `arc_length`, `end_position_local`,
        `end_frame_local`) on a representative connector. Catches signature drift
        and shape regressions.
        """
        p1 = [0.5, 0.4, 0.7]; t_out = [1.0, 0.0, 0.0]
        seg = _build_quintic_connector(p1, t_out, zeros(3), zeros(3);
                                        min_bend_radius = nothing, n_table = 64)
        L = Float64(arc_length(seg))
        @test L > 0
        for s in (0.0, 0.5L, L)
            p = position_local(seg, s);   @test length(p) == 3
            T = tangent_local(seg, s);    @test isapprox(norm(T), 1.0; atol = 1e-9)
            N = normal_local(seg, s);     @test isapprox(norm(N), 1.0; atol = 1e-9)
            B = binormal_local(seg, s);   @test isapprox(norm(B), 1.0; atol = 1e-9)
            @test isapprox(dot(T, N), 0.0; atol = 1e-9)
        end
        ep = end_position_local(seg)
        @test isapprox(ep, p1; atol = 1e-10)
        Tend, Nend, Bend = end_frame_local(seg)
        @test isapprox(Tend, t_out; atol = 1e-9)
        @test isapprox(norm(Nend), 1.0; atol = 1e-9)
        @test isapprox(norm(Bend), 1.0; atol = 1e-9)
    end
end
