using Test
using LinearAlgebra
using Bifrost
using Bifrost.PathGeometry: _qc_nominalize

# -----------------------------------------------------------------------
# Transported (Bishop) frame: physics and guardrail tests
#
# The geometry layer exposes `normal`/`binormal` as the parallel-transported
# (relatively-parallel) pair (e1, e2): zero twist about the tangent, continuous
# along the whole path, anchored at s = 0 by the static lab-frame Gram–Schmidt
# convention of `_initial_frame_from_tangent`. The fiber layer projects the
# curvature vector onto that pair to orient bend birefringence. These tests pin
# the defining properties of the transport and the Jones-level consequences
# that the old ∫τ_geom (Frenet) gauge got wrong (issues #88, #89).
# -----------------------------------------------------------------------

const BISHOP_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like (bishop tests)",
)

const BISHOP_XS_ELLIPTICAL = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like elliptical (bishop tests)",
    ellipticity_axis_ratio = 1.05,
)

const BISHOP_λ = 1550e-9

# Mixed path exercising every authored segment type plus both connector kinds.
function _bishop_mixed_path()
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.2)
    bend!(sb; radius = 0.05, angle = π / 2, axis_angle = π / 6)
    catenary!(sb; a = 0.3, length = 0.2)
    helix!(sb; radius = 0.04, pitch = 0.03, turns = 1.5)
    jumpby!(sb; delta = (0.05, 0.03, 0.1))
    jumpto!(sb; point = (0.0, 0.0, 0.9), min_bend_radius = 0.01)
    return build(sb)
end

# Smallest absolute angular difference modulo 2π.
_wrap_angle(x) = mod(x + π, 2π) - π

@testset "Bishop — (T, e1, e2) orthonormal over a mixed path" begin
    # T-GUARDRAIL: the transported frame must stay right-handed orthonormal on
    # every segment type, including both connector kinds.
    b = _bishop_mixed_path()
    L = Float64(_qc_nominalize(s_end(b)))
    for s in range(1e-6, L - 1e-6; length = 401)
        T = tangent(b, s); e1 = normal(b, s); e2 = binormal(b, s)
        @test abs(norm(T) - 1) < 1e-12
        @test abs(norm(e1) - 1) < 1e-9
        @test abs(dot(T, e1)) < 1e-9
        @test norm(cross(T, e1) - e2) < 1e-9
    end
end

@testset "Bishop — zero twist about the tangent" begin
    # T-PHYSICS: the defining property of a relatively-parallel field is
    # de1/ds ∥ T, i.e. ⟨de1/ds, e2⟩ = 0 — the frame never rotates about the
    # tangent. Verified by central differences inside every segment (segment
    # joints are breakpoints, where e1 is continuous but de1/ds may jump).
    b = _bishop_mixed_path()
    L = Float64(_qc_nominalize(s_end(b)))
    bps = breakpoints(b)
    h = 1e-6
    for s in range(1e-4, L - 1e-4; length = 801)
        minimum(abs.(bps .- s)) > 2h || continue
        e2 = binormal(b, s)
        de1 = (normal(b, s + h) - normal(b, s - h)) ./ (2h)
        @test abs(dot(de1, e2)) < 1e-6   # rad/m, FD-limited
    end
end

@testset "Bishop — static lab-frame anchor at s = 0" begin
    # T-PHYSICS: e1(0) is the Gram–Schmidt projection of the world axis least
    # aligned with the start tangent — a static lab convention, independent of
    # any curvature that follows. Derivation: ref = argmin_axis |T·axis|;
    # e1(0) = (ref − (ref·T) T) / ‖·‖.
    cases = (
        ((0.0, 0.0, 1.0), [1.0, 0.0, 0.0]),   # T = ẑ → ref = x̂ → e1 = x̂
        ((1.0, 0.0, 0.0), [0.0, 1.0, 0.0]),   # T = x̂ → ref = ŷ → e1 = ŷ
        ((0.0, 1.0, 0.0), [1.0, 0.0, 0.0]),   # T = ŷ → ref = x̂ → e1 = x̂
    )
    for (tdir, e1_expected) in cases
        sb = SubpathBuilder(); start!(sb; outgoing_tangent = tdir)
        straight!(sb; length = 0.3)
        # A bend AFTER the start must not affect the anchor.
        bend!(sb; radius = 0.05, angle = π / 3, axis_angle = 0.4)
        seal!(sb)
        b = build(sb)
        @test normal(b, 0.0) ≈ e1_expected atol = 1e-12
    end
    # Oblique tangent: e1 = (x̂ − (x̂·T)T)/‖·‖ with T = (1,1,1)/√3.
    T = [1.0, 1.0, 1.0] ./ sqrt(3.0)
    e1_expected = normalize([1.0, 0.0, 0.0] .- dot([1.0, 0.0, 0.0], T) .* T)
    sb = SubpathBuilder(); start!(sb; outgoing_tangent = Tuple(T))
    straight!(sb; length = 0.3)
    seal!(sb)
    @test normal(build(sb), 0.0) ≈ e1_expected atol = 1e-12
end

@testset "Bishop — perpendicular corner has distinct Jones axes (#88)" begin
    # T-PHYSICS: two equal 90° bends of radius R. In the transported gauge the
    # first bend's axis is θ_b = 0 and the second's is θ_b = axis_angle (the
    # transport carries e1 through the joint unchanged; the curvature direction
    # jumps physically).
    #
    # Same plane (axis_angle = 0): one retarder, K = (iΔβ/2)σ3, total phase
    #   Δβ·2ℓ  →  J = diag(e^{iΔβℓ}, e^{−iΔβℓ}),  ℓ = Rπ/2.
    # Perpendicular planes (axis_angle = π/2): cos2φ flips sign, so the second
    # generator is exactly −K and the bend contributions cancel:
    #   J = e^{ℓ(−K)} e^{ℓK} = 𝟙.
    R = 0.05
    ℓ = R * π / 2
    function corner_J(axis2)
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 0.1)
        bend!(sb; radius = R, angle = π / 2)
        bend!(sb; radius = R, angle = π / 2, axis_angle = axis2)
        straight!(sb; length = 0.1)
        seal!(sb)
        f = Fiber(build(sb); cross_section = BISHOP_XS, T_ref_K = 297.15)
        J, _ = propagate_fiber(f; λ_m = BISHOP_λ, verbose = false)
        return J, f
    end
    J_same, f_same = corner_J(0.0)
    J_perp, f_perp = corner_J(π / 2)

    Δβ = bending_birefringence(BISHOP_XS, BISHOP_λ, 297.15; bend_radius_m = R)
    @test J_same[1, 1] ≈ exp(1im * Δβ * ℓ) atol = 1e-8
    @test abs(J_same[1, 2]) < 1e-8
    @test J_perp ≈ Matrix{ComplexF64}(I, 2, 2) atol = 1e-8
    @test maximum(abs.(J_same .- J_perp)) > 1e-3   # physically distinct

    # The projected curvature components rotate with the bend plane.
    s_b1 = 0.1 + 0.5ℓ
    s_b2 = 0.1 + ℓ + 0.5ℓ
    bc1 = bend_components(f_perp.path, s_b1)
    bc2 = bend_components(f_perp.path, s_b2)
    @test bc1.kx ≈ 1 / R atol = 1e-9
    @test abs(bc1.ky) < 1e-9
    @test abs(bc2.kx) < 1e-9
    @test bc2.ky ≈ 1 / R atol = 1e-9
end

@testset "Bishop — S-bend equals one continuous bend (mod-π axis)" begin
    # T-PHYSICS: bend(θ, axis 0) then bend(θ, axis π): the second curvature
    # vector points along −e1, i.e. axis angle π. A linear retarder depends on
    # its axis only mod π (cos2φ, sin2φ are π-periodic), so the S-bend's Jones
    # matrix equals that of ONE continuous bend of angle 2θ — even though the
    # geometry differs completely. (The old gauge already treated these as the
    # same axis but for the wrong reason — no torsion anywhere; in the
    # transported gauge the equality is the physically derived statement.)
    R = 0.05
    θ = π / 2
    sb1 = SubpathBuilder(); start!(sb1)
    bend!(sb1; radius = R, angle = θ)
    bend!(sb1; radius = R, angle = θ, axis_angle = π)
    seal!(sb1)
    f_s = Fiber(build(sb1); cross_section = BISHOP_XS, T_ref_K = 297.15)
    J_s, _ = propagate_fiber(f_s; λ_m = BISHOP_λ, verbose = false)

    sb2 = SubpathBuilder(); start!(sb2)
    bend!(sb2; radius = R, angle = 2θ)
    seal!(sb2)
    f_c = Fiber(build(sb2); cross_section = BISHOP_XS, T_ref_K = 297.15)
    J_c, _ = propagate_fiber(f_c; λ_m = BISHOP_λ, verbose = false)

    @test J_s ≈ J_c atol = 1e-8

    # The S-joint axis flip is visible in the curvature projections: (κ, 0)
    # before the joint, (−κ, 0) after.
    ℓ = R * θ
    bc_a = bend_components(f_s.path, 0.5ℓ)
    bc_b = bend_components(f_s.path, 1.5ℓ)
    @test bc_a.kx ≈ 1 / R atol = 1e-9
    @test bc_b.kx ≈ -1 / R atol = 1e-9
end

@testset "Bishop — helix bend axis advances at the torsion rate" begin
    # T-PHYSICS: inside a helix the Frenet normal (and with it the curvature
    # vector k⃗ = κN̂) rotates about the tangent at the geometric torsion rate
    # τ = h/(R² + h²) relative to any relatively-parallel field, so
    #   θ_b(s) − θ_b(0⁺) = τ·s.
    # Across a helix → bend(axis_angle = 0) joint the curvature direction is
    # continuous (the bend inherits the helix's end normal), so θ_b must not
    # jump — the one case the old hybrid gauge got right, kept as regression.
    R = 0.05
    pitch = 0.02
    turns = 3.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.3)
    helix!(sb; radius = R, pitch = pitch, turns = turns)
    bend!(sb; radius = R, angle = π / 2, axis_angle = 0.0)
    seal!(sb)
    f = Fiber(build(sb); cross_section = BISHOP_XS, T_ref_K = 297.15)

    h = pitch / (2π)
    τ = h / (R^2 + h^2)
    L_helix = turns * 2π * sqrt(R^2 + h^2)
    s0 = 0.3 + 1e-9
    θ0 = bend_geometry(f, s0).theta_b
    for frac in (0.2, 0.5, 0.8)
        s = 0.3 + frac * L_helix
        θ = bend_geometry(f, s).theta_b
        @test abs(_wrap_angle(θ - θ0 - τ * frac * L_helix)) < 1e-6
    end
    # Joint continuity helix → bend(axis_angle = 0).
    s_joint = 0.3 + L_helix
    θ_before = bend_geometry(f, s_joint - 1e-9).theta_b
    θ_after  = bend_geometry(f, s_joint + 1e-9).theta_b
    @test abs(_wrap_angle(θ_after - θ_before)) < 1e-6
end

@testset "Bishop — pre-refactor regression anchors" begin
    # T-SIM-REGRESSION: values captured on the pre-Bishop code (commit
    # 507c58c) for two reference fibers at 1550 nm, circular core.
    #
    # (a) Planar straight(0.5) → bend(R=0.05, π/2) → straight(0.5): the old
    #     gauge is exactly correct on planar single-subpath paths (torsion
    #     phase ≡ 0 and the anchor coincides), so J must match bit-for-bit
    #     (tolerance 1e-12) and DGD to well below 1 fs.
    # (b) straight(0.3) → helix(R=0.05, pitch=0.02, turns=3) →
    #     bend(R=0.05, π/2, axis 0): old and new gauges differ by a constant
    #     conjugation J ↦ R J Rᵀ from the helix entry onward, so J is NOT
    #     comparable — but DGD is gauge-invariant and must match.
    sbB = SubpathBuilder(); start!(sbB)
    straight!(sbB; length = 0.5)
    bend!(sbB; radius = 0.05, angle = π / 2)
    straight!(sbB; length = 0.5)
    seal!(sbB)
    fibB = Fiber(build(sbB); cross_section = BISHOP_XS, T_ref_K = 297.15)
    JB, GB, _ = propagate_fiber_sensitivity(fibB; λ_m = BISHOP_λ, verbose = false)
    @test JB[1, 1] ≈ +9.994511167114496e-01 - 3.312801388910522e-02im atol = 1e-12
    @test JB[2, 2] ≈ +9.994511167114496e-01 + 3.312801388910522e-02im atol = 1e-12
    @test abs(JB[1, 2]) < 1e-12
    @test abs(JB[2, 1]) < 1e-12
    @test output_dgd_2x2(JB, GB) ≈ 5.661025681660720e-17 atol = 1e-24

    sbA = SubpathBuilder(); start!(sbA)
    straight!(sbA; length = 0.3)
    helix!(sbA; radius = 0.05, pitch = 0.02, turns = 3.0)
    bend!(sbA; radius = 0.05, angle = π / 2, axis_angle = 0.0)
    seal!(sbA)
    fibA = Fiber(build(sbA); cross_section = BISHOP_XS, T_ref_K = 297.15)
    JA, GA, _ = propagate_fiber_sensitivity(fibA; λ_m = BISHOP_λ, verbose = false)
    @test output_dgd_2x2(JA, GA) ≈ 5.530465145008127e-16 atol = 1e-24
end

@testset "Bishop — one vs two Subpaths give identical optics (#89)" begin
    # T-PHYSICS: splitting a path into Subpaths must not change the optics.
    # The transported gauge is made continuous across the boundary by
    # `_resolve_bishop_gauge`, so the frame field and the propagated Jones
    # matrix (elliptical core + twist — sensitive to the relative
    # bend-vs-intrinsic axis angle) must equal the unsplit build's. Exercised
    # for the `:inherit` start and a hand-loaded start.
    #
    # The split suffix uses only frame-independent authoring (straight
    # segments): an `axis_angle`-bearing segment after a boundary would change
    # the authored 3D shape itself, because the successor's construction frame
    # is re-derived from the tangent rather than continued — a pre-existing
    # authoring-layer limitation (subpath concatenation, issues #51/#32),
    # separate from the optical gauge tested here.
    R = 0.05
    L1 = 0.2
    τm = 6.0   # rad/m mechanical twist on the suffix
    sb_one = SubpathBuilder(); start!(sb_one)
    straight!(sb_one; length = L1)
    bend!(sb_one; radius = R, angle = π / 2, axis_angle = π / 6)
    straight!(sb_one; length = 0.25, twist = τm)
    straight!(sb_one; length = 0.15)
    seal!(sb_one)
    p_one = build([Subpath(sb_one)])

    function split_path(handload::Bool)
        sb1 = SubpathBuilder(); start!(sb1)
        straight!(sb1; length = L1)
        bend!(sb1; radius = R, angle = π / 2, axis_angle = π / 6)
        seal!(sb1)
        b1 = build(sb1)
        sb2 = SubpathBuilder()
        if handload
            ep = end_point(b1); et = end_tangent(b1)
            start!(sb2; point = Tuple(Float64.(ep)),
                   outgoing_tangent = Tuple(Float64.(et)))
        else
            start!(sb2, :inherit)
        end
        straight!(sb2; length = 0.25, twist = τm)
        straight!(sb2; length = 0.15)
        seal!(sb2)
        return build([Subpath(sb1), Subpath(sb2)])
    end

    L = Float64(_qc_nominalize(s_end(p_one)))
    for handload in (false, true)
        p_two = split_path(handload)
        @test Float64(_qc_nominalize(s_end(p_two))) ≈ L atol = 1e-9
        # The split is only a meaningful test if the boundary gauge correction
        # is nontrivial (the successor's lab anchor differs from the continued
        # frame after the oblique bend).
        @test abs(p_two.subpaths[2]._bishop_gauge_at_s0) > 0.1
        # Identical geometry...
        for s in range(1e-6, L - 1e-6; length = 51)
            @test position(p_two, s) ≈ position(p_one, s) atol = 1e-9
        end
        # ...identical transported frame field (gauge-continuous across the
        # boundary)...
        for s in range(1e-6, L - 1e-6; length = 101)
            @test normal(p_two, s) ≈ normal(p_one, s) atol = 1e-8
        end
        # ...identical optics with axis-angle-sensitive sources.
        f_one = Fiber(p_one; cross_section = BISHOP_XS_ELLIPTICAL,
                      T_ref_K = 297.15)
        f_two = Fiber(p_two; cross_section = BISHOP_XS_ELLIPTICAL,
                      T_ref_K = 297.15)
        J1, _ = propagate_fiber(f_one; λ_m = BISHOP_λ, verbose = false)
        J2, _ = propagate_fiber(f_two; λ_m = BISHOP_λ, verbose = false)
        @test J1 ≈ J2 atol = 1e-8
    end
end

@testset "Bishop — no torsion spike through a near-straight connector" begin
    # T-PHYSICS: the transported frame rotates at rate |k⃗·e1| ≤ κ — never
    # faster than the curvature. The old ∫τ_geom gauge diverged here: Frenet
    # torsion carries a ~κ² denominator, so a nearly-straight quintic connector
    # injected enormous spurious axis rotation. Bound: ‖Δe1/Δs‖ ≤ 2·max κ.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.1)
    jumpby!(sb; delta = (1e-4, -2e-5, 0.3))   # almost straight jump
    straight!(sb; length = 0.1)
    seal!(sb)
    b = build(sb)
    L = Float64(_qc_nominalize(s_end(b)))
    κ_max = maximum(curvature(b, s) for s in range(0.0, L; length = 2001))
    @test κ_max < 0.05   # the connector really is nearly straight
    h = 1e-5
    rate_max = maximum(norm(normal(b, s + h) - normal(b, s - h)) / (2h)
                       for s in range(0.11, 0.39; length = 501))
    @test rate_max <= 2 * κ_max + 1e-6
end

@testset "Bishop — pure mechanical twist is untouched" begin
    # T-GUARDRAIL: mechanical twist is material, not geometric. On a straight
    # twisted fiber the transported frame is constant, twist_phase still
    # accumulates ∫τ_m ds = τ·L, and the Jones output is the pure circular-
    # birefringence rotation J = R(Δβc·L/2) (real rotation matrix).
    τm = 4π   # rad/m
    Lf = 0.5
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = Lf, twist = τm)
    seal!(sb)
    b = build(sb)
    @test twist_phase(b, Lf) ≈ τm * Lf atol = 1e-12
    @test normal(b, Lf - 1e-6) ≈ normal(b, 1e-6) atol = 1e-12

    f = Fiber(b; cross_section = BISHOP_XS, T_ref_K = 297.15)
    J, _ = propagate_fiber(f; λ_m = BISHOP_λ, verbose = false)
    Δβc = twisting_birefringence(BISHOP_XS, BISHOP_λ, 297.15;
                                 twist_rate_rad_per_m = τm)
    φ = 0.5 * Δβc * Lf
    @test J ≈ [cos(φ) -sin(φ); sin(φ) cos(φ)] atol = 1e-8
end
