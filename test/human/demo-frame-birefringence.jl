# =====================================================================
# demo-frame-birefringence.jl
#
# Visual companion to the Bishop-frame (parallel-transport) refactor
# (issues #88, #89; see docs/src/frame-and-gauge.md). One function,
# one HTML artifact:
#
#     demo_frame_birefringence() → output/frame-birefringence.html
#
# Sections inside the single HTML:
#   1. Frame primer — Frenet normal vs transported e1 on an S-bend.
#   2. All birefringence sources — axis angle and magnitude per source,
#      extracted from the actual fiber generators K(s).
#   3. Pathological cases in the old ∫τ_geom gauge, each shown alongside
#      the transported-frame resolution.
#   4. Where the old hybrid was right — helix → bend(axis_angle = 0).
#   5. The static lab-frame anchor convention.
#
# Run:
#
#     include("test/human/demo-frame-birefringence.jl")
#     demo_frame_birefringence()
# =====================================================================

using Bifrost
using Bifrost.PathGeometry: _qc_nominalize
using LinearAlgebra
using QuadGK
using Printf
include("demo-index-helpers.jl")

const _FB_OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "output")

const _FB_λ = 1550e-9
const _FB_T = 297.15

const _FB_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like (frame demo)",
)

const _FB_XS_ELL = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like elliptical (frame demo)",
    ellipticity_axis_ratio = 1.05,
)

# ---------------------------------------------------------------------
# Presentation helpers (demo-local; JSON-ish serialization for Plotly)
# ---------------------------------------------------------------------

_fb_num(x::Real) = isnan(x) ? "null" : isfinite(x) ? string(Float64(x)) : "null"
_fb_arr(xs) = "[" * join((_fb_num(Float64(x)) for x in xs), ",") * "]"

# One Plotly panel: returns the <div> + <script> fragment. `traces` are JSON
# trace strings; `layout` is a JSON layout string.
function _fb_panel(divid::AbstractString, traces::Vector{String}, layout::AbstractString;
                   height::Int = 420)
    return """
    <div id="$(divid)" style="width:100%;height:$(height)px;"></div>
    <script>
      Plotly.newPlot("$(divid)", [$(join(traces, ","))], $(layout),
                     {responsive: true});
    </script>
    """
end

function _fb_line_trace(x, y; name::AbstractString, dash::AbstractString = "solid",
                        width::Int = 2, yaxis::AbstractString = "y")
    return """{type:"scatter",mode:"lines",x:$(_fb_arr(x)),y:$(_fb_arr(y)),
               name:"$(name)",yaxis:"$(yaxis)",
               line:{dash:"$(dash)",width:$(width)}}"""
end

function _fb_layout_s(title::AbstractString, ylabel::AbstractString;
                      y2label::Union{Nothing, String} = nothing)
    y2 = isnothing(y2label) ? "" :
        """,yaxis2:{title:"$(y2label)",overlaying:"y",side:"right",
                    gridcolor:"#333"}"""
    return """{title:{text:"$(title)"},paper_bgcolor:"#111",plot_bgcolor:"#181818",
               font:{color:"#ddd"},
               xaxis:{title:"s (m)",gridcolor:"#333"},
               yaxis:{title:"$(ylabel)",gridcolor:"#333"},
               legend:{orientation:"h"}$(y2)}"""
end

function _fb_scatter3d_trace(x, y, z; name::AbstractString,
                             mode::AbstractString = "lines",
                             color::AbstractString = "#4db87a", width::Int = 4)
    return """{type:"scatter3d",mode:"$(mode)",x:$(_fb_arr(x)),y:$(_fb_arr(y)),
               z:$(_fb_arr(z)),name:"$(name)",
               line:{color:"$(color)",width:$(width)}}"""
end

const _FB_LAYOUT3D = """{paper_bgcolor:"#111",font:{color:"#ddd"},
    scene:{aspectmode:"data",
           xaxis:{title:"x (m)",gridcolor:"#333",backgroundcolor:"#181818"},
           yaxis:{title:"y (m)",gridcolor:"#333",backgroundcolor:"#181818"},
           zaxis:{title:"z (m)",gridcolor:"#333",backgroundcolor:"#181818"}},
    legend:{orientation:"h"}}"""

# Vector-glyph trace: short segments r → r + scale·v with NaN separators.
function _fb_glyphs(rs, vs; scale::Float64, name::AbstractString,
                    color::AbstractString)
    gx = Float64[]; gy = Float64[]; gz = Float64[]
    for (r, v) in zip(rs, vs)
        any(isnan, v) && continue
        append!(gx, (r[1], r[1] + scale * v[1], NaN))
        append!(gy, (r[2], r[2] + scale * v[2], NaN))
        append!(gz, (r[3], r[3] + scale * v[3], NaN))
    end
    return _fb_scatter3d_trace(gx, gy, gz; name, color, width = 2)
end

# ---------------------------------------------------------------------
# Physics extraction helpers (demo-local ground truth)
# ---------------------------------------------------------------------

# Decompose a lossless Jones generator K = (iΔβ_lin/2)(c2φ σ3 + s2φ σ1)
#                                          + (Δβ_circ/2)(iσ2-form rotation)
# into (Δβ_lin, φ, Δβ_circ): Δβ_lin = 2·hypot(Im K11, Im K12),
# φ = atan(Im K12, Im K11)/2, Δβ_circ = 2·Re K21.
function _fb_K_decompose(K::AbstractMatrix)
    Δβ_lin = 2 * hypot(imag(K[1, 1]), imag(K[1, 2]))
    φ = 0.5 * atan(imag(K[1, 2]), imag(K[1, 1]))
    Δβ_circ = 2 * real(K[2, 1])
    return (Δβ_lin = Δβ_lin, φ = φ, Δβ_circ = Δβ_circ)
end

# Frenet normal by central differences of the path tangent: N = T′/‖T′‖.
# Returns NaN vector where the path is locally straight (‖T′‖ below tol) —
# exactly the undefinedness the transported frame does not suffer from.
function _fb_fs_normal_fd(path, s; h = 1e-6, tol = 1e-3)
    dT = (tangent(path, s + h) .- tangent(path, s - h)) ./ (2h)
    n = norm(dT)
    n < tol && return [NaN, NaN, NaN]
    return dT ./ n
end

# Reconstruction of the deleted ∫τ_geom gauge (old `torsion_phase`) by
# quadrature of the retained `geometric_torsion` diagnostic.
function _fb_old_torsion_phase(path, s)
    val, _ = QuadGK.quadgk(u -> geometric_torsion(path, u), 0.0, s;
                           rtol = 1e-7, atol = 1e-10)
    return val
end

# Fixed-step in-demo Jones propagation for a generator closure (used to
# reconstruct what the OLD gauge would have produced; the library propagator
# only knows the new gauge).
function _fb_propagate(K, s0, s1; n = 4001)
    J = Matrix{ComplexF64}(I, 2, 2)
    h = (s1 - s0) / n
    for i in 1:n
        J = exp(h .* K(s0 + (i - 0.5) * h)) * J
    end
    return J
end

_fb_fmt_J(J) = @sprintf("[%.4f%+.4fim  %.4f%+.4fim; %.4f%+.4fim  %.4f%+.4fim]",
                        real(J[1, 1]), imag(J[1, 1]), real(J[1, 2]), imag(J[1, 2]),
                        real(J[2, 1]), imag(J[2, 1]), real(J[2, 2]), imag(J[2, 2]))

# ---------------------------------------------------------------------
# Section 1 — frame primer: Frenet vs transported on an S-bend
# ---------------------------------------------------------------------

function _fb_section_primer()
    PG = PathGeometry
    R = 0.05
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.06, meta = [PG.Nickname("lead-in")])
    bend!(sb; radius = R, angle = π / 2, meta = [PG.Nickname("90° bend")])
    bend!(sb; radius = R, angle = π / 2, axis_angle = π,
          meta = [PG.Nickname("counter-bend (S)")])
    straight!(sb; length = 0.06, meta = [PG.Nickname("lead-out")])
    seal!(sb)
    b = build(sb)
    L = Float64(_qc_nominalize(s_end(b)))

    ss = range(0.0, L; length = 401)
    xs = [Float64(position(b, s)[1]) for s in ss]
    ys = [Float64(position(b, s)[2]) for s in ss]
    zs = [Float64(position(b, s)[3]) for s in ss]

    sg = range(1e-4, L - 1e-4; length = 61)
    rs = [Float64.(position(b, s)) for s in sg]
    e1s = [Float64.(bishop_e1(b, s)) for s in sg]
    fss = [_fb_fs_normal_fd(b, s) for s in sg]

    traces = String[
        _fb_scatter3d_trace(xs, ys, zs; name = "centerline", color = "#4db87a"),
        _fb_glyphs(rs, e1s; scale = 0.018, name = "transported e1 (continuous)",
                   color = "#5aa0ff"),
        _fb_glyphs(rs, fss; scale = 0.012, name = "Frenet normal (flips, undefined)",
                   color = "#ff6a5a"),
    ]
    panel = _fb_panel("fb-primer", traces, _FB_LAYOUT3D; height = 540)

    return """
    <h2>1. Frame primer — Frenet–Serret vs the transported (Bishop) frame</h2>
    <p>An S-bend: a 90° bend followed by an equal counter-bend
    (<code>axis_angle = π</code>). The <b>Frenet normal</b> (red) always points
    toward the local center of curvature, so it flips by 180° at the S-joint,
    and on the straight lead-in/lead-out it is undefined (no glyphs). The
    <b>transported frame</b> e1 (blue) — defined by zero twist about the
    tangent — is continuous through the joint and perfectly defined on
    straights. Polarization physically follows the transported frame (Rytov's
    law), which is why it is the propagation gauge; the curvature
    <i>direction</i> jumping relative to it at the joint is real physics, and
    enters the generators through the curvature-vector projection.</p>
    $(panel)
    """
end

# ---------------------------------------------------------------------
# Section 2 — all birefringence sources, straight from the generators
# ---------------------------------------------------------------------

function _fb_source_panel(divid, title, fiber, λ; smax = nothing, note = "")
    L = isnothing(smax) ? Float64(_qc_nominalize(s_end(fiber_path(fiber)))) : smax
    K = generator_K(fiber, λ)
    ss = collect(range(1e-6, L - 1e-6; length = 601))
    dec = [_fb_K_decompose(K(s)) for s in ss]
    Δβl = [d.Δβ_lin for d in dec]
    φπ = [d.Δβ_lin > 1e-9 ? d.φ / π : NaN for d in dec]
    Δβc = [d.Δβ_circ for d in dec]

    traces = String[
        _fb_line_trace(ss, Δβl; name = "|Δβ| linear (rad/m)"),
        _fb_line_trace(ss, φπ; name = "axis angle φ/π", yaxis = "y2",
                       dash = "dot"),
    ]
    any(x -> abs(x) > 1e-12, Δβc) &&
        push!(traces, _fb_line_trace(ss, Δβc; name = "Δβ circular (rad/m)",
                                     dash = "dash"))
    layout = _fb_layout_s(title, "Δβ (rad/m)"; y2label = "axis φ/π")
    return """
    <h3>$(title)</h3>
    <p>$(note)</p>
    $(_fb_panel(divid, traces, layout))
    """
end

function _fb_section_sources()
    PG = PathGeometry
    R = 0.06

    # Bend: straight → bend → straight, circular core.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.05, meta = [PG.Nickname("lead-in")])
    bend!(sb; radius = R, angle = π / 2, meta = [PG.Nickname("90° bend")])
    straight!(sb; length = 0.05, meta = [PG.Nickname("lead-out")])
    seal!(sb)
    f_bend = Fiber(build(sb); cross_section = _FB_XS, T_ref_K = _FB_T)

    # Tension: same geometry, axial tension meta on the bend.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.05, meta = [PG.Nickname("lead-in")])
    bend!(sb; radius = R, angle = π / 2,
          meta = [PG.Nickname("tensioned bend"), MCMadd(:tension, 1.0)])
    straight!(sb; length = 0.05, meta = [PG.Nickname("lead-out")])
    seal!(sb)
    f_tens = Fiber(build(sb); cross_section = _FB_XS, T_ref_K = _FB_T)

    # Ellipticity: straight path, elliptical core.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.2, meta = [PG.Nickname("straight")])
    seal!(sb)
    f_ell = Fiber(build(sb); cross_section = _FB_XS_ELL, T_ref_K = _FB_T)

    # Spin: same elliptical core, spun.
    sb = SubpathBuilder(); start!(sb; spin_rate = 8π)
    straight!(sb; length = 0.2, meta = [PG.Nickname("spun straight")])
    seal!(sb)
    f_spin = Fiber(build(sb); cross_section = _FB_XS_ELL, T_ref_K = _FB_T)

    # Twist: elliptical core, mechanically twisted straight.
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.2, twist = 8π, meta = [PG.Nickname("twisted straight")])
    seal!(sb)
    f_twist = Fiber(build(sb); cross_section = _FB_XS_ELL, T_ref_K = _FB_T)

    return """
    <h2>2. The birefringence sources in the transported gauge</h2>
    <p>Each panel evaluates the <i>actual</i> fiber generator K(s) for a fiber
    carrying a single dominant source and decomposes it into a linear
    birefringence magnitude |Δβ|, its eigen-axis angle φ in the (e1, e2)
    frame, and a circular component. Free propagation contributes
    <b>nothing</b>: where the fiber is straight, untwisted, and circular,
    K ≡ 0 — there is no geometric-torsion term in this gauge.</p>
    $(_fb_source_panel("fb-src-bend", "Bend (curvature) birefringence", f_bend, _FB_λ;
        note = "Linear, axis = the curvature direction projected on (e1, e2): " *
               "φ = 0 inside this planar bend, |Δβ| ∝ 1/R², zero on the straights."))
    $(_fb_source_panel("fb-src-tension", "Axial-tension birefringence", f_tens, _FB_λ;
        note = "Shares the bend eigen-axis (∝ 1/R): present only where the " *
               "tensioned segment is also bent."))
    $(_fb_source_panel("fb-src-ell", "Core-ellipticity (intrinsic) birefringence",
        f_ell, _FB_λ;
        note = "Linear, frozen into the glass: axis fixed at the ellipse " *
               "angle (here 0), magnitude independent of s."))
    $(_fb_source_panel("fb-src-spin", "Spun intrinsic birefringence", f_spin, _FB_λ;
        note = "Spin (frozen-in rotation of the glass during draw) carries the " *
               "intrinsic axes with it: φ(s) advances at the spin rate while " *
               "|Δβ| is unchanged. Spin is a material rotation — it does not " *
               "touch the transported frame."))
    $(_fb_source_panel("fb-src-twist", "Mechanically twisted fiber", f_twist, _FB_λ;
        note = "Elastic twist does two things: photoelastic <i>circular</i> " *
               "birefringence (dashed) and co-rotation of the intrinsic axes " *
               "(dotted axis trace) at the twist rate."))
    """
end

# ---------------------------------------------------------------------
# Section 3 — pathological cases of the old gauge, and their resolution
# ---------------------------------------------------------------------

function _fb_patho_corner()
    PG = PathGeometry
    R = 0.05
    ℓ = R * π / 2
    function corner(axis2)
        sb = SubpathBuilder(); start!(sb)
        straight!(sb; length = 0.05, meta = [PG.Nickname("lead-in")])
        bend!(sb; radius = R, angle = π / 2, meta = [PG.Nickname("bend 1")])
        bend!(sb; radius = R, angle = π / 2, axis_angle = axis2,
              meta = [PG.Nickname("bend 2")])
        straight!(sb; length = 0.05, meta = [PG.Nickname("lead-out")])
        seal!(sb)
        return Fiber(build(sb); cross_section = _FB_XS, T_ref_K = _FB_T)
    end
    f_same = corner(0.0)
    f_perp = corner(π / 2)
    L = Float64(_qc_nominalize(s_end(fiber_path(f_perp))))
    ss = collect(range(1e-6, L - 1e-6; length = 601))

    θ_new = [bend_geometry(f_perp, s).theta_b / π for s in ss]
    # Old gauge: θ_b = ∫τ_geom ds ≡ 0 on this all-planar path, for both fibers.
    θ_old = [_fb_old_torsion_phase(fiber_path(f_perp), s) / π for s in ss]
    κs = [curvature(fiber_path(f_perp), s) for s in ss]

    traces = String[
        _fb_line_trace(ss, θ_new; name = "new: θ_b/π (projected k⃗)"),
        _fb_line_trace(ss, θ_old; name = "old: θ_b/π (∫τ_geom — blind to the corner)",
                       dash = "dash"),
        _fb_line_trace(ss, κs ./ maximum(κs); name = "κ/κ_max (segment marker)",
                       dash = "dot", yaxis = "y2"),
    ]
    layout = _fb_layout_s("Perpendicular corner — bend axis angle", "θ_b/π";
                          y2label = "κ/κ_max")
    panel = _fb_panel("fb-patho-corner", traces, layout)

    J_perp_new, _ = propagate_fiber(f_perp; λ_m = _FB_λ, verbose = false)
    J_same_new, _ = propagate_fiber(f_same; λ_m = _FB_λ, verbose = false)
    # Old-gauge reconstruction: both corners produce the same single-axis
    # retarder because the axis never moves.
    Δβ = bending_birefringence(_FB_XS, _FB_λ, _FB_T; bend_radius_m = R)
    K_old = s -> begin
        κ = curvature(fiber_path(f_perp), s)
        Float64(κ) < 1e-9 ? zeros(ComplexF64, 2, 2) :
            ComplexF64[0.5im*Δβ 0.0; 0.0 -0.5im*Δβ]
    end
    J_old = _fb_propagate(K_old, 0.0, L)

    return """
    <h3>3a. Perpendicular-plane corner (issue #88)</h3>
    <p>Two consecutive 90° bends, second one in the <i>perpendicular</i> plane
    (<code>axis_angle = π/2</code>). The old gauge oriented the bend axis by
    ∫τ_geom, which is identically zero on planar segments — the corner was
    invisible, and the perpendicular corner produced <i>exactly</i> the same
    Jones matrix as two bends in the same plane. The projected curvature
    vector sees the corner as a clean π/2 axis jump on an existing breakpoint.
    With equal bend strengths the two perpendicular retarders cancel:</p>
    <pre>
old gauge:  J(perpendicular) = J(same-plane) = $(_fb_fmt_J(J_old))
new gauge:  J(same-plane)    = $(_fb_fmt_J(J_same_new))
            J(perpendicular) = $(_fb_fmt_J(J_perp_new))   ≈ identity (retarders cancel)
    </pre>
    $(panel)
    """
end

function _fb_patho_connector()
    PG = PathGeometry
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.1, meta = [PG.Nickname("lead-in")])
    jumpby!(sb; delta = (1e-4, -2e-5, 0.3), meta = [PG.Nickname("near-straight jump")])
    straight!(sb; length = 0.1, meta = [PG.Nickname("lead-out")])
    seal!(sb)
    b = build(sb)
    L = Float64(_qc_nominalize(s_end(b)))
    ss = collect(range(1e-6, L - 1e-6; length = 801))

    τfs = [geometric_torsion(b, s) for s in ss]
    # Cumulative trapezoid of τ_FS — reconstructs the old torsion_phase.
    ϕ_old = similar(τfs)
    ϕ_old[firstindex(τfs)] = 0.0
    for i in eachindex(τfs)
        i == firstindex(τfs) && continue
        ϕ_old[i] = ϕ_old[i-1] + 0.5 * (τfs[i] + τfs[i-1]) * (ss[i] - ss[i-1])
    end
    f = Fiber(b; cross_section = _FB_XS, T_ref_K = _FB_T)
    θ_new = [bend_geometry(f, s).theta_b for s in ss]
    κmax = maximum(curvature(b, s) for s in ss)

    traces = String[
        _fb_line_trace(ss, ϕ_old; name = "old: ∫τ_FS ds (spikes near κ→0)"),
        _fb_line_trace(ss, θ_new; name = "new: θ_b (bounded by transport)"),
        _fb_line_trace(ss, τfs; name = "Frenet torsion τ_FS (rad/m)",
                       dash = "dot", yaxis = "y2"),
    ]
    layout = _fb_layout_s("Near-straight quintic connector", "phase / angle (rad)";
                          y2label = "τ_FS (rad/m)")
    return """
    <h3>3b. Torsion spike in a near-straight connector</h3>
    <p>A <code>jumpby!</code> connector that is almost straight
    (transverse offset 10⁻⁴ m over 0.3 m). Frenet torsion carries a ~κ²
    denominator, so where the connector's curvature passes near zero, τ_FS
    spikes by orders of magnitude and the old ∫τ_FS bend-axis phase jumps by
    O(1) rad through a region that is optically almost inert
    (max κ = $(@sprintf("%.4f", κmax)) m⁻¹). The transported frame cannot
    rotate faster than κ, so the new axis angle stays bounded and smooth.</p>
    $(_fb_panel("fb-patho-connector", traces, layout))
    """
end

function _fb_patho_boundary()
    PG = PathGeometry
    R = 0.05
    L1 = 0.2
    # One Subpath vs two (split after the oblique bend; straight suffix keeps
    # the authored 3D geometry identical — see test_bishop_frame.jl).
    sb_one = SubpathBuilder(); start!(sb_one)
    straight!(sb_one; length = L1, meta = [PG.Nickname("lead-in")])
    bend!(sb_one; radius = R, angle = π / 2, axis_angle = π / 6,
          meta = [PG.Nickname("oblique bend")])
    straight!(sb_one; length = 0.4, meta = [PG.Nickname("suffix")])
    seal!(sb_one)
    p_one = build([Subpath(sb_one)])

    sb1 = SubpathBuilder(); start!(sb1)
    straight!(sb1; length = L1, meta = [PG.Nickname("lead-in")])
    bend!(sb1; radius = R, angle = π / 2, axis_angle = π / 6,
          meta = [PG.Nickname("oblique bend")])
    seal!(sb1)
    sb2 = SubpathBuilder(); start!(sb2, :inherit)
    straight!(sb2; length = 0.4, meta = [PG.Nickname("suffix")])
    seal!(sb2)
    p_two = build([Subpath(sb1), Subpath(sb2)])

    L = Float64(_qc_nominalize(s_end(p_one)))
    s_bnd = L1 + R * π / 2
    ss = collect(range(1e-6, L - 1e-6; length = 601))

    # Angle of each build's e1 relative to the single-subpath reference.
    function rel_angle(p, s)
        e1r = bishop_e1(p_one, s); e2r = bishop_e2(p_one, s)
        e1 = bishop_e1(p, s)
        return atan(dot(e1, e2r), dot(e1, e1r))
    end
    θ_two = [rel_angle(p_two, s) for s in ss]
    # Naive (old) behavior: each Subpath re-anchored from its own tangent.
    # Emulate by querying the second Subpath standalone past the boundary.
    sb2n = SubpathBuilder()
    et = Float64.(end_tangent(build(Subpath(sb1))))
    ep = Float64.(end_point(build(Subpath(sb1))))
    start!(sb2n; point = Tuple(ep), outgoing_tangent = Tuple(et))
    straight!(sb2n; length = 0.4)
    seal!(sb2n)
    b2_alone = build(sb2n)
    θ_naive = [s <= s_bnd ? 0.0 :
               (e1 = bishop_e1(b2_alone, s - s_bnd);
                atan(dot(e1, bishop_e2(p_one, s)), dot(e1, bishop_e1(p_one, s))))
               for s in ss]

    traces = String[
        _fb_line_trace(ss, θ_naive ./ π;
                       name = "old: per-Subpath re-anchored e1 (gauge jump)"),
        _fb_line_trace(ss, θ_two ./ π;
                       name = "new: gauge-resolved PathBuilt e1 (continuous)"),
    ]
    layout = _fb_layout_s("Subpath boundary — e1 angle relative to unsplit build",
                          "angle/π")

    f_one = Fiber(p_one; cross_section = _FB_XS_ELL, T_ref_K = _FB_T)
    f_two = Fiber(p_two; cross_section = _FB_XS_ELL, T_ref_K = _FB_T)
    J1, _ = propagate_fiber(f_one; λ_m = _FB_λ, verbose = false)
    J2, _ = propagate_fiber(f_two; λ_m = _FB_λ, verbose = false)

    return """
    <h3>3c. Subpath boundary gauge break (issue #89)</h3>
    <p>The same physical fiber authored as one Subpath or as two. Each Subpath
    builds with its own static lab anchor; without correction, the frame —
    and with it every birefringence axis — jumps by a constant angle at the
    boundary (here $(@sprintf("%.3f", abs(p_two.subpaths[2]._bishop_gauge_at_s0)))
    rad after the oblique bend). <code>build(::Vector)</code> now resolves one
    constant gauge rotation per Subpath (exact: parallel transport commutes
    with constant transverse rotations), so the split build is optically
    identical to the unsplit one:</p>
    <pre>
J(one Subpath)  = $(_fb_fmt_J(J1))
J(two Subpaths) = $(_fb_fmt_J(J2))    max|ΔJ| = $(@sprintf("%.2e", maximum(abs.(J1 .- J2))))
    </pre>
    <p>Caveat (authoring layer, unchanged by this refactor): segments whose
    parameters reference the construction frame (<code>axis_angle</code>,
    <code>jumpby!</code> deltas) are interpreted in a re-derived frame after a
    boundary, so frame-dependent authoring after a split can change the 3D
    shape itself — the subpath-concatenation epic (#51/#32) owns that.</p>
    $(_fb_panel("fb-patho-boundary", traces, layout))
    """
end

function _fb_patho_framerate()
    PG = PathGeometry
    R = 0.05
    pitch = 0.02
    sb = SubpathBuilder(); start!(sb; spin_rate = 2π)
    straight!(sb; length = 0.1, meta = [PG.Nickname("lead-in")])
    helix!(sb; radius = R, pitch = pitch, turns = 2.0, meta = [PG.Nickname("helix")])
    straight!(sb; length = 0.1, meta = [PG.Nickname("lead-out")])
    seal!(sb)
    b = build(sb)
    L = Float64(_qc_nominalize(s_end(b)))
    ss = collect(range(1e-6, L - 1e-6; length = 601))

    old_rate = [geometric_torsion(b, s) + spin_rate(b, s) for s in ss]
    new_rate = [spin_rate(b, s) for s in ss]
    traces = String[
        _fb_line_trace(ss, old_rate; name = "old display: τ_geom + spin (jumps)"),
        _fb_line_trace(ss, new_rate; name = "new display: spin_rate (material only)"),
    ]
    layout = _fb_layout_s("Spun straight → helix → straight", "rate (rad/m)")
    return """
    <h3>3d. The \"frame rate\" jump at helix boundaries (issue #24)</h3>
    <p>The old diagnostics displayed a combined frame rotation rate
    τ_geom + spin, which jumps discontinuously at helix entry/exit — the
    visual artifact reported in #24 (<code>helix-mcm-spinning</code>). In the
    transported gauge the frame itself never rotates about the tangent, so the
    only meaningful rotation rate is the material spin — continuous across
    the helix. (∫τ_geom survives only as the shape diagnostic
    <code>total_torsion</code>, with no optical role.)</p>
    $(_fb_panel("fb-patho-framerate", traces, layout))
    """
end

function _fb_section_pathologies()
    return """
    <h2>3. Pathological cases of the old gauge — and their resolution</h2>
    <p>Each panel reconstructs the old ∫τ_geom behavior in-demo (the quantity
    itself is no longer part of the optics) and overlays the transported-gauge
    result. A further item in this family, the S-shape conformity rejection
    (issue #62: Frenet-flip artifacts in endpoint curvature checks), is
    analyzed in <code>docs/src/frame-and-gauge.md</code> — its fix lives in
    the authoring conformity check, not in the optical gauge.</p>
    $(_fb_patho_corner())
    $(_fb_patho_connector())
    $(_fb_patho_boundary())
    $(_fb_patho_framerate())
    """
end

# ---------------------------------------------------------------------
# Section 4 — where the old hybrid was right: helix → bend(axis_angle 0)
# ---------------------------------------------------------------------

function _fb_section_hybrid()
    PG = PathGeometry
    R = 0.05
    pitch = 0.02
    turns = 3.0
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = 0.3, meta = [PG.Nickname("lead-in")])
    helix!(sb; radius = R, pitch = pitch, turns = turns, meta = [PG.Nickname("helix")])
    bend!(sb; radius = R, angle = π / 2, meta = [PG.Nickname("exit bend")])
    seal!(sb)
    f = Fiber(build(sb); cross_section = _FB_XS, T_ref_K = _FB_T)
    b = fiber_path(f)
    L = Float64(_qc_nominalize(s_end(b)))
    ss = collect(range(1e-6, L - 1e-6; length = 801))

    h = pitch / (2π)
    τ = h / (R^2 + h^2)
    L_helix = turns * 2π * sqrt(R^2 + h^2)
    θ_new = [bend_geometry(f, s).theta_b for s in ss]
    # Analytic reference inside the helix: θ_b(0⁺) + τ·(s − 0.3), wrapped.
    θ0 = bend_geometry(f, 0.3 + 1e-9).theta_b
    θ_ref = [(s > 0.3 && s < 0.3 + L_helix) ?
             mod(θ0 + τ * (s - 0.3) + π, 2π) - π : NaN for s in ss]

    J, G, _ = propagate_fiber_sensitivity(f; λ_m = _FB_λ, verbose = false)
    dgd_new = output_dgd_2x2(J, G)

    traces = String[
        _fb_line_trace(ss, [mod(θ + π, 2π) - π for θ in θ_new];
                       name = "θ_b (projected, wrapped)"),
        _fb_line_trace(ss, θ_ref; name = "analytic: θ_b(0⁺) + τ·s", dash = "dot"),
    ]
    layout = _fb_layout_s("Helix → bend(axis_angle = 0)", "θ_b (rad)")
    return """
    <h2>4. Where the old hybrid was right — helix → bend(axis_angle = 0)</h2>
    <p>Inside a helix the curvature direction rotates at exactly the torsion
    rate τ relative to the transported frame, so the old ∫τ_geom phase tracked
    the bend axis <i>correctly</i> there; and because the placement loop chains
    each segment's Frenet end frame, a following bend with
    <code>axis_angle = 0</code> continued the axis smoothly. Both gauges agree
    on every gauge-invariant observable on this geometry — kept as a
    regression so the fix cannot silently break the case the old code got
    right. DGD here: $(@sprintf("%.6e", dgd_new)) s vs pre-refactor
    5.530465e-16 s (gauge-invariant; equality enforced to 10⁻²⁴ s in
    <code>test/test_bishop_frame.jl</code>). The Jones matrices differ between
    gauges only by a fixed conjugation R·J·Rᵀ (a constant axis offset from the
    helix entry), which no polarization observable distinguishes.</p>
    $(_fb_panel("fb-hybrid", traces, layout))
    """
end

# ---------------------------------------------------------------------
# Section 5 — the anchor convention
# ---------------------------------------------------------------------

function _fb_section_anchor()
    PG = PathGeometry
    cases = [
        ((0.0, 0.0, 1.0), "launch ∥ ẑ"),
        ((1.0, 0.0, 0.0), "launch ∥ x̂"),
        ((1.0, 1.0, 1.0) ./ sqrt(3.0), "launch ∥ (1,1,1)/√3"),
    ]
    rows = String[]
    traces = String[]
    colors = ("#4db87a", "#5aa0ff", "#ff6a5a")
    for (i, (tdir, label)) in enumerate(cases)
        sb = SubpathBuilder(); start!(sb; outgoing_tangent = Tuple(Float64.(tdir)))
        straight!(sb; length = 0.05, meta = [PG.Nickname("lead-in")])
        bend!(sb; radius = 0.05, angle = π / 2, meta = [PG.Nickname("bend")])
        straight!(sb; length = 0.05, meta = [PG.Nickname("lead-out")])
        seal!(sb)
        f = Fiber(build(sb); cross_section = _FB_XS, T_ref_K = _FB_T)
        b = fiber_path(f)
        e10 = Float64.(bishop_e1(b, 0.0))
        J, G, _ = propagate_fiber_sensitivity(f; λ_m = _FB_λ, verbose = false)
        dgd = output_dgd_2x2(J, G)
        push!(rows, @sprintf("%-22s e1(0) = (%+.3f, %+.3f, %+.3f)   DGD = %.6e s",
                             label, e10[1], e10[2], e10[3], dgd))
        L = Float64(_qc_nominalize(s_end(b)))
        ss = range(0.0, L; length = 201)
        push!(traces, _fb_scatter3d_trace(
            [Float64(position(b, s)[1]) for s in ss],
            [Float64(position(b, s)[2]) for s in ss],
            [Float64(position(b, s)[3]) for s in ss];
            name = label, color = colors[i]))
        rs = [Float64.(position(b, s)) for s in range(0.0, L; length = 13)]
        e1s = [Float64.(bishop_e1(b, s)) for s in range(0.0, L; length = 13)]
        push!(traces, _fb_glyphs(rs, e1s; scale = 0.015,
                                 name = "e1 — $(label)", color = colors[i]))
    end
    return """
    <h2>5. The anchor convention — static, lab-frame, curvature-blind</h2>
    <p>e1(0) is fixed by one static rule: Gram–Schmidt of the world axis least
    aligned with the launch tangent (x̂ when launching along ẑ). It depends
    only on the launch direction — never on the curvature that follows — so it
    is trivial to reason about and reproducible by hand. The same bend
    launched along three different directions gets three different anchors,
    and every polarization observable is unchanged (DGD identical below);
    only the gauge in which J is reported rotates. To express J on specific
    lab axes, conjugate by the constant rotation between (e1, e2) and those
    axes at the fiber ends.</p>
    <pre>$(join(rows, "\n"))</pre>
    $(_fb_panel("fb-anchor", traces, _FB_LAYOUT3D; height = 540))
    """
end

# ---------------------------------------------------------------------
# Assembly — one function, one HTML file
# ---------------------------------------------------------------------

"""
    demo_frame_birefringence(; output) -> String

Write the multi-section visual companion to the Bishop-frame refactor to
`output/frame-birefringence.html` and return the path. See the file header
for the section list.
"""
function demo_frame_birefringence(;
    output::AbstractString = joinpath(_FB_OUTPUT_DIR, "frame-birefringence.html"),
)
    isdir(_FB_OUTPUT_DIR) || mkpath(_FB_OUTPUT_DIR)
    html = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Transported (Bishop) frame and birefringence</title>
  <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
  <style>
    body { font-family: sans-serif; max-width: 1100px; margin: 2em auto;
           background: #111; color: #ddd; }
    h1   { font-size: 1.6em; border-bottom: 1px solid #444;
           padding-bottom: 0.3em; }
    h2   { font-size: 1.25em; margin-top: 2.0em; color: #4db87a; }
    h3   { font-size: 1.05em; margin-top: 1.5em; color: #ddd; }
    p    { line-height: 1.45; }
    pre  { background: #181818; color: #bbb; padding: 0.8em;
           overflow-x: auto; font-size: 0.85em; }
    code { color: #bbb; }
  </style>
</head>
<body>
  <h1>Transported (Bishop) frame and the birefringence gauge</h1>
  <p>Companion artifact to the Bishop-frame refactor (issues #88, #89) and to
  <code>docs/src/frame-and-gauge.md</code>. All "old gauge" curves are
  reconstructed in-demo from retained diagnostics; everything else queries the
  library as shipped.</p>
  $(_fb_section_primer())
  $(_fb_section_sources())
  $(_fb_section_pathologies())
  $(_fb_section_hybrid())
  $(_fb_section_anchor())
</body>
</html>
"""
    open(output, "w") do io
        write(io, html)
    end
    println("Wrote ", output)
    return String(output)
end

# Demo-index registration helper: one entry, one artifact.
function demo_frame_birefringence_entries()
    println("[ demo ] demo_frame_birefringence")
    path = demo_frame_birefringence()
    return Tuple{String, String, String, String}[(
        "frame",
        basename(path),
        path,
        "Frenet vs transported (Bishop) frame: primer, all birefringence " *
        "sources, old-gauge pathologies with resolutions, the helix " *
        "regression case, and the lab-frame anchor convention.",
    )]
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_frame_birefringence()
end
