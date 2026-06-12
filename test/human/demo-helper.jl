# =====================================================================
# demo-helper.jl — plotting and logistics support for bifrost-demos.ipynb
# =====================================================================
#
# Everything visual lives here so the notebook cells can stay focused on the
# Bifrost API. The renderers implement the visual-technique catalog in
# demo-intent.md (§V1–V8): the interactive 3D path inspector, side-by-side
# variant rows, baseline-vs-modified overlays, ensemble scatters, the adaptive
# step-doubling panels, and the benchmark chart/table. Helpers here may be as
# abstract as convenient; readability lives in the notebook.
#
# Usage from the notebook:
#
#     include("demo-helper.jl")
#     using .DemoHelper

module DemoHelper

using LinearAlgebra
using Printf
using Markdown
using PlotlyJS

using Bifrost
using Bifrost.PathGeometry: SubpathBuilt, PathBuilt, PlacedSegment, Nickname,
    sample_path, s_end, s_offsets, path_length, arc_length, position, total_spin,
    _qc_nominalize, _all_placed_segs

export dh_path_inspector,
    dh_variant_row, dh_variant_row_2d, dh_overlay_compare, dh_try_build,
    dh_jones_to_stokes, dh_stokes_ensemble, dh_ensemble_scatter, dh_poincare_equatorial,
    dh_adaptive_panels, dh_benchmark_chart, dh_benchmark_table,
    dh_temperature_ptf_row, dh_temperature_scatter_row

# ---------------------------------------------------------------------
# Palette and layout conventions (demo-intent.md “cross-cutting conventions”)
# ---------------------------------------------------------------------

const COLOR_SUBJECT = "#e15555"   # red: the element under study
const COLOR_CONTEXT = "#66cc66"   # green: fixed context segments
const COLOR_FAINT   = "#555555"   # gray: incidental terminal connectors
const COLOR_BASE    = "#111111"   # black: baseline curve in overlays

const _DARK_BG  = "#111111"
const _DARK_AX  = attr(gridcolor = "#333", zerolinecolor = "#333", color = "#aaa")
# 3D scene axes also need their background panes darkened (the default template
# draws light panes that override scene.bgcolor).
const _DARK_AX3 = attr(gridcolor = "#333", zerolinecolor = "#333", color = "#aaa",
                       showbackground = true, backgroundcolor = _DARK_BG)

_nom(x) = Float64(_qc_nominalize(x))

"""
    _dark3d(title) -> Layout

Shared dark 3D scene layout: equal aspect, muted axes, legend styled for dark.
"""
function _dark3d(title::AbstractString; height::Int = 560, ranges = nothing,
                 camera_eye = nothing)
    xax = merge(_DARK_AX3, attr(title = "x (m)"))
    yax = merge(_DARK_AX3, attr(title = "y (m)"))
    zax = merge(_DARK_AX3, attr(title = "z (m)"))
    if ranges === nothing
        scene = attr(bgcolor = _DARK_BG, xaxis = xax, yaxis = yax, zaxis = zax,
                     aspectmode = "data")
    else
        # Fixed bounding box: axes never re-range as slider/cursor traces move.
        xr, yr, zr = ranges
        spans = (xr[2] - xr[1], yr[2] - yr[1], zr[2] - zr[1])
        m = max(spans...)
        scene = attr(bgcolor = _DARK_BG,
            xaxis = merge(xax, attr(range = [xr...], autorange = false)),
            yaxis = merge(yax, attr(range = [yr...], autorange = false)),
            zaxis = merge(zax, attr(range = [zr...], autorange = false)),
            aspectmode = "manual",
            aspectratio = attr(x = spans[1] / m, y = spans[2] / m, z = spans[3] / m))
    end
    # A larger camera eye distance zooms the initial view out. `camera_eye`
    # overrides Plotly's default eye of (1.25, 1.25, 1.25) when supplied.
    camera_eye === nothing ||
        (scene[:camera] = attr(eye = attr(x = camera_eye[1], y = camera_eye[2],
                                          z = camera_eye[3])))
    return Layout(
        title = attr(text = title, font = attr(color = "#eee", size = 15)),
        paper_bgcolor = _DARK_BG,
        scene = scene,
        height = height,
        margin = attr(l = 0, r = 0, t = 42, b = 0),
        legend = attr(font = attr(color = "#ddd")),
    )
end

function _dark2d(title::AbstractString; height::Int = 470)
    return Layout(
        title = attr(text = title, font = attr(color = "#eee", size = 15)),
        paper_bgcolor = _DARK_BG, plot_bgcolor = _DARK_BG,
        xaxis = merge(_DARK_AX, attr(title = "x (m)")),
        yaxis = merge(_DARK_AX, attr(title = "z (m)", scaleanchor = "x",
                                     scaleratio = 1)),
        height = height,
        margin = attr(l = 55, r = 15, t = 42, b = 45),
        legend = attr(font = attr(color = "#ddd")),
    )
end

function _light2d(title::AbstractString; height::Int = 470)
    return Layout(
        title = attr(text = title, font = attr(size = 15)),
        paper_bgcolor = "#fafafa", plot_bgcolor = "#fafafa",
        xaxis = attr(title = "x (m)", gridcolor = "#e5e5e5"),
        yaxis = attr(title = "z (m)", gridcolor = "#e5e5e5", scaleanchor = "x",
                     scaleratio = 1),
        height = height,
        margin = attr(l = 55, r = 15, t = 42, b = 45),
    )
end

# ---------------------------------------------------------------------
# Path access shims (uniform across SubpathBuilt / PathBuilt)
# ---------------------------------------------------------------------

"""
    _placed(path) -> Vector{NamedTuple}

All placed segments of a built path in global order. Each entry carries the
segment, its global arc-length window `(s0, s1)`, its `Nickname` (or `""`),
and whether it is a terminal connector of its subpath.
"""
function _placed(path::SubpathBuilt)
    out = NamedTuple[]
    segs = _all_placed_segs(path)
    for (i, ps) in enumerate(segs)
        s0 = _nom(ps.s_offset_eff)
        push!(out, (segment = ps.segment, s0 = s0,
                    s1 = s0 + _nom(arc_length(ps.segment)),
                    nickname = _nickname(ps.segment),
                    terminal = i == length(segs)))
    end
    return out
end

function _placed(path::PathBuilt)
    out = NamedTuple[]
    offs = s_offsets(path)
    for (k, sp) in enumerate(path.subpaths)
        segs = _all_placed_segs(sp)
        for (i, ps) in enumerate(segs)
            s0 = offs[k] + _nom(ps.s_offset_eff)
            push!(out, (segment = ps.segment, s0 = s0,
                        s1 = s0 + _nom(arc_length(ps.segment)),
                        nickname = _nickname(ps.segment),
                        terminal = i == length(segs)))
        end
    end
    return out
end

function _nickname(seg)
    if hasfield(typeof(seg), :meta)
        for m in seg.meta
            m isa Nickname && return String(m.label)
        end
    end
    return ""
end

"""
    _segment_curve(path, s0, s1; n = 96) -> (x, y, z)

Sample the global-frame centerline of one segment window.
"""
function _segment_curve(path, s0::Float64, s1::Float64; n::Int = 96)
    ss = range(s0, min(s1, _nom(s_end(path))); length = n)
    xs = Vector{Float64}(undef, n); ys = similar(xs); zs = similar(xs)
    for (i, s) in enumerate(ss)
        p = position(path, s)
        xs[i] = _nom(p[1]); ys[i] = _nom(p[2]); zs[i] = _nom(p[3])
    end
    return xs, ys, zs
end

# ---------------------------------------------------------------------
# V1 — interactive 3D path inspector
# ---------------------------------------------------------------------

"""
    dh_path_inspector(path; title = "", fidelity = 1.0, nsteps = 50) -> Plot

Interactive inspector for a built path (`SubpathBuilt` or `PathBuilt`):

- faint centerline with dot markers color-graded by arc length,
- open circles at segment boundaries, green start / red end dots,
- `Nickname` labels floated at segment midpoints,
- an arc-length slider that scrubs a cursor along the path, carrying a
  translucent normal–binormal plane square, a T̂/N̂/B̂ axis triad
  (orange/blue/green), a red in-plane spin arrow at the accumulated spin
  phase, and a live readout (s; x,y,z; κ; τ_geom; τ_spin; ∫τ_spin ds).

The frame shown is the parallel-transport frame reported by `sample_path`.
"""
function dh_path_inspector(path; title::AbstractString = "", fidelity::Float64 = 1.0,
                        nsteps::Int = 50)
    hi = _nom(s_end(path))
    ps = sample_path(path, 0.0, hi; fidelity = fidelity)
    n  = ps.n
    sv   = [_nom(sm.s) for sm in ps.samples]
    xv   = [_nom(sm.position[1]) for sm in ps.samples]
    yv   = [_nom(sm.position[2]) for sm in ps.samples]
    zv   = [_nom(sm.position[3]) for sm in ps.samples]
    tv   = [Float64[_nom(sm.tangent[i])  for i in 1:3] for sm in ps.samples]
    nv   = [Float64[_nom(sm.normal[i])   for i in 1:3] for sm in ps.samples]
    bv   = [Float64[_nom(sm.binormal[i]) for i in 1:3] for sm in ps.samples]
    κv   = [_nom(sm.curvature) for sm in ps.samples]
    τgv  = [_nom(sm.geometric_torsion) for sm in ps.samples]
    τsv  = [_nom(sm.spin_rate) for sm in ps.samples]
    # Cumulative spin phase by trapezoid over the dense samples.
    ϕv = zeros(n)
    for i in 2:n
        ϕv[i] = ϕv[i-1] + 0.5 * (τsv[i] + τsv[i-1]) * (sv[i] - sv[i-1])
    end

    diam = max(maximum(xv) - minimum(xv), maximum(yv) - minimum(yv),
               maximum(zv) - minimum(zv), 1e-9)
    halfw = 0.10 * diam      # transverse square half-width
    axlen = 0.14 * diam      # T/N/B axis length

    traces = GenericTrace[]
    push!(traces, scatter3d(x = xv, y = yv, z = zv, mode = "lines",
        line = attr(color = "#bbbbbb", width = 3), name = "path",
        hoverinfo = "skip", showlegend = false))
    push!(traces, scatter3d(x = xv, y = yv, z = zv, mode = "markers",
        marker = attr(size = 2.4, color = sv, colorscale = "Reds",
                      showscale = false),
        name = "arc length", hoverinfo = "skip", showlegend = false))

    segs = _placed(path)
    bx = Float64[]; by = Float64[]; bz = Float64[]
    lx = Float64[]; ly = Float64[]; lz = Float64[]; ltxt = String[]
    for sg in segs
        p = position(path, min(sg.s0, hi))
        push!(bx, _nom(p[1])); push!(by, _nom(p[2])); push!(bz, _nom(p[3]))
        if !isempty(sg.nickname)
            pm = position(path, min(0.5 * (sg.s0 + sg.s1), hi))
            push!(lx, _nom(pm[1])); push!(ly, _nom(pm[2]))
            push!(lz, _nom(pm[3]) + 0.05 * diam)
            push!(ltxt, sg.nickname)
        end
    end
    push!(traces, scatter3d(x = bx, y = by, z = bz, mode = "markers",
        marker = attr(size = 4.5, color = "rgba(0,0,0,0)",
                      line = attr(color = "#888", width = 2), symbol = "circle"),
        name = "segment joins", hoverinfo = "skip", showlegend = false))
    push!(traces, scatter3d(x = [xv[1]], y = [yv[1]], z = [zv[1]],
        mode = "markers", marker = attr(size = 5, color = "#2ca02c"),
        name = "start", showlegend = false))
    push!(traces, scatter3d(x = [xv[end]], y = [yv[end]], z = [zv[end]],
        mode = "markers", marker = attr(size = 5, color = "#d62728"),
        name = "end", showlegend = false))
    isempty(ltxt) || push!(traces, scatter3d(x = lx, y = ly, z = lz,
        mode = "text", text = ltxt, textfont = attr(color = "#ddd", size = 12),
        hoverinfo = "skip", showlegend = false))

    # --- cursor traces (updated by the slider); start at index 1 ---
    cursor_geo(i) = begin
        p = [xv[i], yv[i], zv[i]]
        nh = nv[i]; bh = bv[i]; th = tv[i]
        corners = [p .+ halfw .* (sx .* nh .+ sy .* bh)
                   for (sx, sy) in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
        spin = cos(ϕv[i]) .* nh .+ sin(ϕv[i]) .* bh
        (p = p, t = th, nh = nh, bh = bh, corners = corners, spin = spin)
    end
    g = cursor_geo(1)
    plane = mesh3d(
        x = [c[1] for c in g.corners], y = [c[2] for c in g.corners],
        z = [c[3] for c in g.corners], i = [0, 0], j = [1, 2], k = [2, 3],
        color = "#7da7e8", opacity = 0.30, name = "N–B plane",
        hoverinfo = "skip", showlegend = false)
    axis_tr(v, col, nm) = scatter3d(
        x = [g.p[1], g.p[1] + axlen * v[1]], y = [g.p[2], g.p[2] + axlen * v[2]],
        z = [g.p[3], g.p[3] + axlen * v[3]], mode = "lines",
        line = attr(color = col, width = 7), name = nm, showlegend = true,
        hoverinfo = "skip")
    tT = axis_tr(g.t, "#ff7f0e", "T̂")
    tN = axis_tr(g.nh, "#1f77b4", "N̂")
    tB = axis_tr(g.bh, "#2ca02c", "B̂")
    tS = axis_tr(g.spin, "#d62728", "∫τ_spin")
    dot = scatter3d(x = [g.p[1]], y = [g.p[2]], z = [g.p[3]], mode = "markers",
        marker = attr(size = 4, color = "#000",
                      line = attr(color = "#fff", width = 1)),
        name = "cursor", showlegend = false, hoverinfo = "skip")
    cursor_start = length(traces)            # 0-based index of `plane`
    append!(traces, [plane, tT, tN, tB, tS, dot])
    cursor_idx = collect(cursor_start:(cursor_start + 5))

    readout(i) = @sprintf(
        "s = %.4f m | x,y,z = %.4f, %.4f, %.4f | κ = %.3g 1/m | τ_geom = %.3g | τ_spin = %.3g rad/m | ∫τ_spin = %.3f rad (%.1f°)",
        sv[i], xv[i], yv[i], zv[i], κv[i], τgv[i], τsv[i], ϕv[i], rad2deg(ϕv[i]))

    steps = []
    for k in 1:nsteps
        # Step targets are uniform in arc length (sample density is not).
        s_target = hi * (k - 1) / (nsteps - 1)
        i = clamp(searchsortedfirst(sv, s_target), 1, n)
        gi = cursor_geo(i)
        seg_x(v) = [gi.p[1], gi.p[1] + axlen * v[1]]
        seg_y(v) = [gi.p[2], gi.p[2] + axlen * v[2]]
        seg_z(v) = [gi.p[3], gi.p[3] + axlen * v[3]]
        restyle = attr(
            x = [[c[1] for c in gi.corners], seg_x(gi.t), seg_x(gi.nh),
                 seg_x(gi.bh), seg_x(gi.spin), [gi.p[1]]],
            y = [[c[2] for c in gi.corners], seg_y(gi.t), seg_y(gi.nh),
                 seg_y(gi.bh), seg_y(gi.spin), [gi.p[2]]],
            z = [[c[3] for c in gi.corners], seg_z(gi.t), seg_z(gi.nh),
                 seg_z(gi.bh), seg_z(gi.spin), [gi.p[3]]],
        )
        relayout = attr(annotations = [attr(
            text = readout(i), x = 0.5, y = -0.02, xref = "paper",
            yref = "paper", showarrow = false,
            font = attr(color = "#ccc", size = 11))])
        # Label only a few steps so the slider tick row stays legible.
        lbl = (k - 1) % max(1, nsteps ÷ 8) == 0 ? @sprintf("%.2f", sv[i]) : " "
        push!(steps, attr(method = "update", label = lbl,
                          args = [restyle, relayout, cursor_idx]))
    end

    # Bounding box over everything the scene can show: the path, the cursor
    # plane/triad at any s, and floated labels. Axes stay fixed while scrubbing.
    ext = max(axlen, halfw * sqrt(2)) + 0.03 * diam
    ranges = ((minimum(xv) - ext, maximum(xv) + ext),
              (minimum(yv) - ext, maximum(yv) + ext),
              (minimum(zv) - ext, maximum(zv) + max(ext, 0.08 * diam)))
    layout = _dark3d(title; ranges = ranges)
    layout[:annotations] = [attr(text = readout(1), x = 0.5, y = -0.02,
        xref = "paper", yref = "paper", showarrow = false,
        font = attr(color = "#ccc", size = 11))]
    layout[:sliders] = [attr(steps = steps, active = 0, pad = attr(t = 28),
        currentvalue = attr(prefix = "s = ", suffix = " m",
                            font = attr(color = "#ccc")),
        font = attr(color = "#888"))]
    return plot(traces, layout)
end

# ---------------------------------------------------------------------
# V2 / V4 — side-by-side variant rows (3D and 2D)
# ---------------------------------------------------------------------

"""
    dh_variant_row(variants; spacing = 2.5, title = "", gray_terminal = false)

Side-by-side comparison of path variants in 3D. `variants` is a vector of
`(label, path, red)` where `red` lists 1-based placed-segment indices to draw
in the subject color (use `Int[]` for none). Variants are offset along +x by
`spacing`; non-red segments draw green; with `gray_terminal = true` each
subpath's terminal connector draws faint gray instead.
"""
function dh_variant_row(variants::Vector; spacing::Real = 2.5,
                     title::AbstractString = "", gray_terminal::Bool = false,
                     height::Int = 560)
    traces = GenericTrace[]
    sx = Float64[]; sy = Float64[]; sz = Float64[]
    lx = Float64[]; ly = Float64[]; lz = Float64[]; ltxt = String[]
    for (k, (label, path, red)) in enumerate(variants)
        dx = (k - 1) * Float64(spacing)
        zmax = -Inf
        for (i, sg) in enumerate(_placed(path))
            x, y, z = _segment_curve(path, sg.s0, sg.s1)
            color = i in red ? COLOR_SUBJECT :
                    (gray_terminal && sg.terminal ? COLOR_FAINT : COLOR_CONTEXT)
            push!(traces, scatter3d(x = x .+ dx, y = y, z = z, mode = "lines",
                line = attr(color = color, width = 6),
                name = "$(label) — seg $(i)", showlegend = false))
            zmax = max(zmax, maximum(z))
        end
        p0 = position(path, 0.0)
        push!(sx, _nom(p0[1]) + dx); push!(sy, _nom(p0[2])); push!(sz, _nom(p0[3]))
        push!(lx, dx); push!(ly, 0.0); push!(lz, zmax + 0.18)
        push!(ltxt, String(label))
    end
    push!(traces, scatter3d(x = sx, y = sy, z = sz, mode = "markers",
        marker = attr(size = 5, color = COLOR_CONTEXT), name = "start",
        showlegend = false))
    push!(traces, scatter3d(x = lx, y = ly, z = lz, mode = "text", text = ltxt,
        textfont = attr(color = "#ddd", size = 13), showlegend = false,
        hoverinfo = "skip"))
    # Start 25% more zoomed out than Plotly's default eye of (1.25, 1.25, 1.25).
    return plot(traces, _dark3d(title; height = height,
                                camera_eye = (1.5625, 1.5625, 1.5625)))
end

"""
    dh_variant_row_2d(variants; spacing = 2.5, title = "")

The 2D (x–z plane) version of [`dh_variant_row`](@ref) for paths lying in y = 0:
variant columns offset along x, red/green semantic coloring, white dots at
segment ends, per-variant labels at the top.
"""
function dh_variant_row_2d(variants::Vector; spacing::Real = 2.5,
                        title::AbstractString = "", height::Int = 470)
    traces = GenericTrace[]
    annos = []
    ex = Float64[]; ez = Float64[]
    for (k, (label, path, red)) in enumerate(variants)
        dx = (k - 1) * Float64(spacing)
        zmax = -Inf
        for (i, sg) in enumerate(_placed(path))
            x, _, z = _segment_curve(path, sg.s0, sg.s1)
            color = i in red ? COLOR_SUBJECT : COLOR_CONTEXT
            push!(traces, scatter(x = x .+ dx, y = z, mode = "lines",
                line = attr(color = color, width = 2.5), showlegend = false,
                name = "$(label) — seg $(i)"))
            push!(ex, x[end] + dx); push!(ez, z[end])
            zmax = max(zmax, maximum(z))
        end
        p0 = position(path, 0.0)
        push!(traces, scatter(x = [_nom(p0[1]) + dx], y = [_nom(p0[3])],
            mode = "markers", marker = attr(size = 9, color = COLOR_CONTEXT),
            showlegend = false))
        push!(annos, attr(x = dx, y = zmax, yshift = 16, text = String(label),
            showarrow = false, font = attr(color = "#ddd", size = 12)))
    end
    push!(traces, scatter(x = ex, y = ez, mode = "markers",
        marker = attr(size = 5, color = "#fff",
                      line = attr(color = "#000", width = 0.6)),
        showlegend = false))
    layout = _dark2d(title; height = height)
    layout[:annotations] = annos
    return plot(traces, layout)
end

# ---------------------------------------------------------------------
# V3 — baseline-vs-modified overlay (2D, x–z)
# ---------------------------------------------------------------------

"""
    dh_overlay_compare(baseline, modified; title = "", legend = :topleft)

Overlay two built paths in the x–z plane on a light background: baseline
black, modified red, open circles at every placed-segment boundary, and a
legend reporting both total lengths and the percent change.
"""
function dh_overlay_compare(baseline, modified; title::AbstractString = "",
                         legend::Symbol = :topleft, height::Int = 470)
    Lb = _nom(path_length(baseline))
    Lm = _nom(path_length(modified))
    traces = GenericTrace[]
    for (path, color, name) in (
        (baseline, COLOR_BASE, @sprintf("baseline:  L = %.4f m", Lb)),
        (modified, COLOR_SUBJECT,
         @sprintf("modified:  L = %.4f m  (Δ = %+.1f%%)", Lm, 100 * (Lm / Lb - 1))),
    )
        first_seg = true
        ex = Float64[]; ez = Float64[]
        for sg in _placed(path)
            x, _, z = _segment_curve(path, sg.s0, sg.s1)
            push!(traces, scatter(x = x, y = z, mode = "lines",
                line = attr(color = color, width = 2.5), name = name,
                legendgroup = name, showlegend = first_seg))
            append!(ex, (x[1], x[end])); append!(ez, (z[1], z[end]))
            first_seg = false
        end
        push!(traces, scatter(x = ex, y = ez, mode = "markers",
            marker = attr(size = 6, color = "#fafafa",
                          line = attr(color = color, width = 1.5)),
            legendgroup = name, showlegend = false))
    end
    layout = _light2d(title; height = height)
    layout[:legend] = legend === :bottomright ?
        attr(x = 0.98, y = 0.02, xanchor = "right", yanchor = "bottom",
             bgcolor = "rgba(255,255,255,0.8)") :
        attr(x = 0.02, y = 0.98, xanchor = "left", yanchor = "top",
             bgcolor = "rgba(255,255,255,0.8)")
    return plot(traces, layout)
end

# ---------------------------------------------------------------------
# V5 — failure-tolerant sweeps
# ---------------------------------------------------------------------

"""
    dh_try_build(f) -> (value, error)

Run the zero-argument builder `f`, returning `(result, nothing)` on success or
`(nothing, err)` if the build throws (e.g. a `min_bend_radius` infeasibility).
"""
function dh_try_build(f)
    try
        return (f(), nothing)
    catch err
        return (nothing, err)
    end
end

# ---------------------------------------------------------------------
# V6 — ensemble post-processing and scatter plots
# ---------------------------------------------------------------------

"""
    dh_jones_to_stokes(J; input = [1, 0]) -> NamedTuple

Apply a scalar `ComplexF64` Jones matrix to `input`, normalize, and return
`(s1, s2, s3, dlp, angle_deg)` where `dlp = hypot(s1, s2)` and `angle_deg` is
the linear polarization angle folded to [0°, 180°).
"""
function dh_jones_to_stokes(J::AbstractMatrix; input = ComplexF64[1.0, 0.0])
    ψ = Matrix{ComplexF64}(J) * ComplexF64.(input)
    ψ ./= sqrt(abs2(ψ[1]) + abs2(ψ[2]))
    s1 = real(abs2(ψ[1]) - abs2(ψ[2]))
    s2 = 2 * real(ψ[1] * conj(ψ[2]))
    s3 = -2 * imag(ψ[1] * conj(ψ[2]))
    return (s1 = s1, s2 = s2, s3 = s3, dlp = hypot(s1, s2),
            angle_deg = mod(0.5 * atan(s2, s1) * 180 / π, 180.0))
end

_particle(x::Real, k::Int) =
    hasfield(typeof(x), :particles) ? getfield(x, :particles)[k] : Float64(x)
_particle(z::Complex, k::Int) =
    complex(_particle(real(z), k), _particle(imag(z), k))

"""
    dh_stokes_ensemble(J_p; input = [1, 0]) -> NamedTuple of vectors

Slice a Jones matrix with MCM `Particles` entries into its per-particle scalar
matrices and convert each to Stokes observables. Returns
`(s1, s2, s3, dlp, angle_deg, N)`.
"""
function dh_stokes_ensemble(J_p::AbstractMatrix; input = ComplexF64[1.0, 0.0])
    r = real(J_p[1, 1])
    N = hasfield(typeof(r), :particles) ? length(getfield(r, :particles)) : 1
    s1 = zeros(N); s2 = zeros(N); s3 = zeros(N)
    dlp = zeros(N); ang = zeros(N)
    for k in 1:N
        Jk = ComplexF64[_particle(J_p[1, 1], k) _particle(J_p[1, 2], k);
                        _particle(J_p[2, 1], k) _particle(J_p[2, 2], k)]
        st = dh_jones_to_stokes(Jk; input = input)
        s1[k] = st.s1; s2[k] = st.s2; s3[k] = st.s3
        dlp[k] = st.dlp; ang[k] = st.angle_deg
    end
    return (s1 = s1, s2 = s2, s3 = s3, dlp = dlp, angle_deg = ang, N = N)
end

"""
    dh_ensemble_scatter(x, series; color = x, title, xlab, ylab, yrange = nothing)

Scatter one or more observables against an uncertain input, marker color
encoding `color` (Viridis). `series` is a vector of `(name, values)` pairs.
"""
function dh_ensemble_scatter(x::AbstractVector, series::Vector;
                          color = x, title::AbstractString = "",
                          xlab::AbstractString = "Temperature (°C)",
                          ylab::AbstractString = "", yrange = nothing,
                          height::Int = 430)
    traces = GenericTrace[]
    for (j, (name, y)) in enumerate(series)
        push!(traces, scatter(x = x, y = y, mode = "markers", name = String(name),
            marker = attr(size = 5, color = color, colorscale = "Viridis",
                          opacity = 0.85, showscale = j == 1,
                          colorbar = attr(title = "T (°C)",
                                          tickfont = attr(color = "#ccc"))),
            ))
    end
    layout = Layout(
        title = attr(text = title, font = attr(color = "#ddd", size = 14)),
        paper_bgcolor = "#1a1a1a", plot_bgcolor = "#1a1a1a",
        font = attr(color = "#ccc"),
        xaxis = merge(_DARK_AX, attr(title = xlab)),
        yaxis = merge(_DARK_AX, attr(title = ylab)),
        height = height, margin = attr(l = 60, r = 20, t = 45, b = 50),
    )
    yrange === nothing || (layout[:yaxis][:range] = yrange)
    return plot(traces, layout)
end

"""
    dh_poincare_equatorial(s1, s2; color, title = "") -> Plot

S1–S2 equatorial projection of the Poincaré sphere: per-particle markers
colored by the uncertain input, axes locked to [−1, 1] at unit aspect, with
the unit circle (fully polarized linear states) drawn for reference.
"""
function dh_poincare_equatorial(s1::AbstractVector, s2::AbstractVector;
                             color = nothing, title::AbstractString = "",
                             height::Int = 470)
    θ = range(0, 2π; length = 181)
    traces = GenericTrace[
        scatter(x = cos.(θ), y = sin.(θ), mode = "lines",
                line = attr(color = "#444", width = 1), showlegend = false,
                hoverinfo = "skip"),
        scatter(x = s1, y = s2, mode = "markers",
            marker = attr(size = 5, color = color, colorscale = "Viridis",
                          opacity = 0.8, showscale = color !== nothing,
                          colorbar = attr(title = "T (°C)",
                                          tickfont = attr(color = "#ccc"))),
            showlegend = false),
    ]
    layout = Layout(
        title = attr(text = title, font = attr(color = "#ddd", size = 14)),
        paper_bgcolor = "#1a1a1a", plot_bgcolor = "#1a1a1a",
        font = attr(color = "#ccc"),
        xaxis = merge(_DARK_AX, attr(title = "S1", range = [-1.05, 1.05])),
        yaxis = merge(_DARK_AX, attr(title = "S2", range = [-1.05, 1.05],
                                     scaleanchor = "x", scaleratio = 1)),
        height = height, margin = attr(l = 60, r = 20, t = 45, b = 50),
    )
    return plot(traces, layout)
end

# ---------------------------------------------------------------------
# V7 — adaptive step-doubling panels
# ---------------------------------------------------------------------

"""
    dh_adaptive_panels(records, K_norm, s0, s1; rtol, atol, components = [])

Three stacked panels over arc length for `collect_adaptive_steps` output:
(1) step size h (log y) — accepted green dots, rejected red ✕, scaled ‖K(s)‖
as a shaded band; (2) err/tol per trial (log y) with the acceptance threshold
dashed at 1; (3) the generator's labeled component coefficients.
`components` is a vector of `(label, fn, color)`.
"""
function dh_adaptive_panels(records::Vector, K_norm, s0::Real, s1::Real;
                         rtol::Float64 = 1e-6, atol::Float64 = 1e-9,
                         components::Vector = [], height::Int = 760,
                         title::AbstractString = "Adaptive step-doubling")
    acc = [r for r in records if r.accepted]
    rej = [r for r in records if !r.accepted]
    ss  = collect(range(Float64(s0), Float64(s1); length = 481))
    Kn  = [Float64(K_norm(s)) for s in ss]
    hmax = maximum(r.h for r in records)
    hmin = minimum(r.h for r in records)
    # Clamp the shaded band away from zero so the log axis keeps a sane range.
    band = max.(Kn ./ maximum(Kn) .* hmax, 0.5 * hmin)

    p = make_subplots(rows = 3, cols = 1, shared_xaxes = true,
                      vertical_spacing = 0.06,
                      row_heights = [0.42, 0.33, 0.25])
    add_trace!(p, scatter(x = ss, y = band, mode = "lines", fill = "tozeroy",
        line = attr(color = "rgba(120,150,220,0.5)", width = 1),
        fillcolor = "rgba(120,150,220,0.15)", name = "‖K(s)‖ (scaled)"),
        row = 1, col = 1)
    add_trace!(p, scatter(x = [r.s_start for r in acc], y = [r.h for r in acc],
        mode = "markers", marker = attr(size = 5, color = "#2ca02c"),
        name = "accepted step"), row = 1, col = 1)
    add_trace!(p, scatter(x = [r.s_start for r in rej], y = [r.h for r in rej],
        mode = "markers",
        marker = attr(size = 7, color = "#8b1a1a", symbol = "x"),
        name = "rejected step"), row = 1, col = 1)

    ratio(r) = r.err_abs / r.tol
    add_trace!(p, scatter(x = [r.s_start for r in acc], y = ratio.(acc),
        mode = "markers", marker = attr(size = 4, color = "#2ca02c"),
        name = "err/tol (accepted)", showlegend = false), row = 2, col = 1)
    add_trace!(p, scatter(x = [r.s_start for r in rej], y = ratio.(rej),
        mode = "markers", marker = attr(size = 5, color = "#d62728"),
        name = "err/tol (rejected)", showlegend = false), row = 2, col = 1)
    add_trace!(p, scatter(x = [Float64(s0), Float64(s1)], y = [1.0, 1.0],
        mode = "lines", line = attr(color = "#888", dash = "dash", width = 1),
        name = "threshold (err = tol)"), row = 2, col = 1)

    for (label, fn, color) in components
        add_trace!(p, scatter(x = ss, y = [Float64(fn(s)) for s in ss],
            mode = "lines", line = attr(color = color, width = 2),
            name = String(label)), row = 3, col = 1)
    end

    relayout!(p, height = height, title = attr(text = title),
        template = "plotly_white",
        legend = attr(orientation = "h", y = -0.08),
        yaxis  = attr(type = "log", title = "step size h"),
        yaxis2 = attr(type = "log", title = "err / tol"),
        yaxis3 = attr(title = "coefficient"),
        xaxis3 = attr(title = "s (arc length)"))
    return p
end

# ---------------------------------------------------------------------
# V8 — benchmark chart and table
# ---------------------------------------------------------------------

"""
    dh_benchmark_chart(results; title = "") -> Plot

Grouped log-scale bars of first-call (JIT + run) vs steady-state wall time.
`results` is a vector of NamedTuples `(label, n_particles, first_ms,
steady_ms)` as produced by the benchmark cells.
"""
function dh_benchmark_chart(results::Vector; title::AbstractString = "")
    labels = [String(r.label) for r in results]
    p = plot([
        bar(x = labels, y = [r.first_ms for r in results],
            name = "first-call (JIT + run)", marker_color = "#55aaff"),
        bar(x = labels, y = [r.steady_ms for r in results],
            name = "steady-state (post-JIT min)", marker_color = "#ff8800"),
    ], Layout(
        title = attr(text = title, font = attr(color = "#ddd", size = 14)),
        paper_bgcolor = "#1a1a1a", plot_bgcolor = "#1a1a1a",
        font = attr(color = "#ccc"), barmode = "group",
        xaxis = _DARK_AX,
        yaxis = merge(_DARK_AX, attr(title = "wall time (ms)", type = "log")),
        height = 430, margin = attr(l = 65, r = 20, t = 45, b = 80),
        legend = attr(font = attr(color = "#ccc")),
    ))
    return p
end

"""
    dh_benchmark_table(results) -> Markdown.MD

Markdown table of benchmark numbers with ratios against the first row.
"""
function dh_benchmark_table(results::Vector)
    io = IOBuffer()
    println(io, "| Variant | N particles | First-call (ms) | Steady-state (ms) ",
                "| First ratio | Steady ratio |")
    println(io, "| --- | --- | --- | --- | --- | --- |")
    for r in results
        @printf(io, "| %s | %s | %.1f | %.3f | %.1f× | %.1f× |\n",
                r.label, r.n_particles == 0 ? "—" : string(r.n_particles),
                r.first_ms, r.steady_ms,
                r.first_ms / results[1].first_ms,
                r.steady_ms / results[1].steady_ms)
    end
    return Markdown.parse(String(take!(io)))
end


# ---------------------------------------------------------------------
# V6 (legacy format) — MCM temperature rows, mirroring demo3mcm's layout
# ---------------------------------------------------------------------

"""
    dh_temperature_ptf_row(T_C, st; title = "") -> Plot

Two side-by-side panels in the legacy demo3mcm format: linear polarization
angle vs temperature (Viridis, shared colorbar) and S1/S2/S3/DLP vs
temperature (per-series Reds/Greens/Blues/Oranges colorscales, DLP as
diamonds). `st` is the named tuple from [`dh_stokes_ensemble`](@ref).
"""
function dh_temperature_ptf_row(T_C::AbstractVector, st; title::AbstractString = "")
    p = make_subplots(rows = 1, cols = 2, horizontal_spacing = 0.10,
        subplot_titles = ["Linear polarization angle vs T" "Stokes parameters vs T"])
    add_trace!(p, scatter(x = T_C, y = st.angle_deg, mode = "markers",
        marker = attr(size = 6, color = T_C, colorscale = "Viridis",
                      opacity = 0.85,
                      colorbar = attr(title = "T (°C)", x = 0.42, len = 0.85)),
        name = "pol. angle", showlegend = false), row = 1, col = 1)
    series = (("S1", st.s1, "Reds", "circle"), ("S2", st.s2, "Greens", "circle"),
              ("S3", st.s3, "Blues", "circle"), ("DLP", st.dlp, "Oranges", "diamond"))
    for (nm, y, cs, sym) in series
        add_trace!(p, scatter(x = T_C, y = y, mode = "markers", name = nm,
            marker = attr(size = 6, color = T_C, colorscale = cs, opacity = 0.85,
                          symbol = sym)), row = 1, col = 2)
    end
    relayout!(p, height = 470, paper_bgcolor = "#1a1a1a", plot_bgcolor = "#1a1a1a",
        font = attr(color = "#ccc"), title = attr(text = title),
        legend = attr(font = attr(color = "#ccc"), x = 1.0, y = 1.0),
        xaxis  = merge(_DARK_AX, attr(title = "Temperature (°C)")),
        xaxis2 = merge(_DARK_AX, attr(title = "Temperature (°C)")),
        yaxis  = merge(_DARK_AX, attr(title = "Angle (deg)", range = [0, 180])),
        yaxis2 = merge(_DARK_AX, attr(title = "Stokes parameter / DLP",
                                      range = [-1, 1])))
    return p
end

"""
    dh_temperature_scatter_row(T_C, st; title = "") -> Plot

Two side-by-side panels in the legacy demo3mcm scatter format: polarization
angle vs temperature, and the Poincaré equatorial projection (S1–S2, unit
circle for reference, axes locked to [−1, 1] at unit aspect), both colored by
temperature.
"""
function dh_temperature_scatter_row(T_C::AbstractVector, st;
                                    title::AbstractString = "")
    p = make_subplots(rows = 1, cols = 2, horizontal_spacing = 0.10,
        subplot_titles = ["Polarization angle vs T (MCM scatter)" "Poincaré equatorial projection (S1–S2)"])
    add_trace!(p, scatter(x = T_C, y = st.angle_deg, mode = "markers",
        marker = attr(size = 4, color = T_C, colorscale = "Viridis",
                      opacity = 0.7,
                      colorbar = attr(title = "T (°C)", x = 0.42, len = 0.85)),
        showlegend = false), row = 1, col = 1)
    θ = range(0, 2π; length = 181)
    add_trace!(p, scatter(x = cos.(θ), y = sin.(θ), mode = "lines",
        line = attr(color = "#444", width = 1), hoverinfo = "skip",
        showlegend = false), row = 1, col = 2)
    add_trace!(p, scatter(x = st.s1, y = st.s2, mode = "markers",
        marker = attr(size = 4, color = T_C, colorscale = "Viridis",
                      opacity = 0.7),
        showlegend = false), row = 1, col = 2)
    relayout!(p, height = 470, paper_bgcolor = "#1a1a1a", plot_bgcolor = "#1a1a1a",
        font = attr(color = "#ccc"), title = attr(text = title),
        xaxis  = merge(_DARK_AX, attr(title = "Temperature (°C)")),
        yaxis  = merge(_DARK_AX, attr(title = "Angle (deg)", range = [0, 180])),
        xaxis2 = merge(_DARK_AX, attr(title = "S1", range = [-1.05, 1.05])),
        yaxis2 = merge(_DARK_AX, attr(title = "S2", range = [-1.05, 1.05],
                                      scaleanchor = "x2", scaleratio = 1)))
    return p
end

end # module DemoHelper
