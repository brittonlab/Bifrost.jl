# Interactive Plotly HTML for a `SubpathBuilt`; see write_path_geometry_plot3d.
# This file depends only on the `PathGeometry` submodule and the Plotly CDN.

using LinearAlgebra

# Reuse the canonical Bifrost.PathGeometry submodule rather than re-including
# path-geometry.jl here. This file is loaded inside `Bifrost.Plots`, which can
# see its sibling submodules via `..`. The `PathGeometry.X` qualified references
# below resolve against that submodule.
using ..PathGeometry

# ---------------------------------------------------------------------------
# Sampling (uses `frame` from PathGeometry — analytic Frenet data on `Path`)
# ---------------------------------------------------------------------------

"""
    _expand(ps::PathGeometry.PathSample) → NamedTuple

Unpack a `PathSample` into the flat named arrays expected by the HTML template.
"""
function _expand(ps::PathGeometry.PathSample)
    n = ps.n
    s        = [smpl.s                    for smpl in ps.samples]
    x        = [smpl.position[1]          for smpl in ps.samples]
    y        = [smpl.position[2]          for smpl in ps.samples]
    z        = [smpl.position[3]          for smpl in ps.samples]
    tx       = [smpl.tangent[1]           for smpl in ps.samples]
    ty       = [smpl.tangent[2]           for smpl in ps.samples]
    tz       = [smpl.tangent[3]           for smpl in ps.samples]
    nx       = [smpl.normal[1]            for smpl in ps.samples]
    ny       = [smpl.normal[2]            for smpl in ps.samples]
    nz       = [smpl.normal[3]            for smpl in ps.samples]
    bx       = [smpl.binormal[1]          for smpl in ps.samples]
    by       = [smpl.binormal[2]          for smpl in ps.samples]
    bz       = [smpl.binormal[3]          for smpl in ps.samples]
    kappa    = [smpl.curvature            for smpl in ps.samples]
    tau_geom = [smpl.geometric_torsion    for smpl in ps.samples]
    tau_spin  = [smpl.spin_rate       for smpl in ps.samples]
    return (; s, x, y, z, tx, ty, tz, nx, ny, nz, bx, by, bz, kappa, tau_geom, tau_spin)
end

# ---------------------------------------------------------------------------
# Plotly HTML helpers (small JSON serializers for embedded numeric arrays)
# ---------------------------------------------------------------------------

function _plot_scalar(x::Real)
    if hasfield(typeof(x), :particles)
        particles = getfield(x, :particles)
        return Float64(sum(particles) / length(particles))
    end
    return Float64(x)
end

function _js_real(x::Real)
    xf = _plot_scalar(x)
    if isnan(xf)
        return "NaN"
    elseif xf == Inf
        return "Infinity"
    elseif xf == -Inf
        return "-Infinity"
    else
        return string(xf)
    end
end

_js_array(xs::AbstractVector{<:Real}) = "[" * join(_js_real.(xs), ", ") * "]"

function _js_string_array(xs::AbstractVector{<:AbstractString})
    escaped = replace.(xs, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
    return "[" * join(["\"" * x * "\"" for x in escaped], ", ") * "]"
end


# -----------------------------------------------------------------------
# Type-specific helpers: collect segment-boundary marker positions and
# segment nickname labels for a SubpathBuilt or PathBuilt within
# [s1f, s2f]. Returned as plain Float64/String vectors ready for HTML
# embedding.
# -----------------------------------------------------------------------

function _collect_segment_boundaries(path::PathGeometry.SubpathBuilt,
                                     s1f::Float64, s2f::Float64)
    seg_bx = Float64[]; seg_by = Float64[]; seg_bz = Float64[]
    seg_bound_hover = String[]
    placed = path.placed_segments
    if length(placed) >= 2
        for i in 2:length(placed)
            sj = Float64(PathGeometry._qc_nominalize(placed[i].s_offset_eff))
            (sj < s1f || sj > s2f) && continue
            p = PathGeometry.position(path, sj)
            push!(seg_bx, p[1]); push!(seg_by, p[2]); push!(seg_bz, p[3])
            push!(seg_bound_hover,
                  "Segment boundary<br>s = $(sj) m<br>x, y, z = $(p[1]), $(p[2]), $(p[3])")
        end
    end
    sj_t = Float64(PathGeometry._qc_nominalize(path.jumpto_placed.s_offset_eff))
    if sj_t >= s1f && sj_t <= s2f
        p = PathGeometry.position(path, sj_t)
        push!(seg_bx, p[1]); push!(seg_by, p[2]); push!(seg_bz, p[3])
        push!(seg_bound_hover,
              "Terminal connector start<br>s = $(sj_t) m<br>x, y, z = $(p[1]), $(p[2]), $(p[3])")
    end
    return (; seg_bx, seg_by, seg_bz, seg_bound_hover)
end

function _collect_segment_boundaries(path::PathGeometry.PathBuilt,
                                     s1f::Float64, s2f::Float64)
    seg_bx = Float64[]; seg_by = Float64[]; seg_bz = Float64[]
    seg_bound_hover = String[]
    offs = PathGeometry.s_offsets(path)
    n = length(path.subpaths)

    for (i, sp) in enumerate(path.subpaths)
        sp_off = offs[i]
        # Interior segment joins inside Subpath i.
        placed = sp.placed_segments
        if length(placed) >= 2
            for k in 2:length(placed)
                sj = sp_off + Float64(PathGeometry._qc_nominalize(placed[k].s_offset_eff))
                (sj < s1f || sj > s2f) && continue
                p = PathGeometry.position(path, sj)
                push!(seg_bx, p[1]); push!(seg_by, p[2]); push!(seg_bz, p[3])
                push!(seg_bound_hover,
                      "Subpath $i segment join<br>s = $(sj) m<br>x, y, z = $(p[1]), $(p[2]), $(p[3])")
            end
        end
        # Terminal connector start of Subpath i.
        sj_t = sp_off +
            Float64(PathGeometry._qc_nominalize(sp.jumpto_placed.s_offset_eff))
        if sj_t >= s1f && sj_t <= s2f
            p = PathGeometry.position(path, sj_t)
            push!(seg_bx, p[1]); push!(seg_by, p[2]); push!(seg_bz, p[3])
            push!(seg_bound_hover,
                  "Subpath $i terminal connector start<br>s = $(sj_t) m<br>x, y, z = $(p[1]), $(p[2]), $(p[3])")
        end
    end
    # Subpath-to-Subpath boundary markers (start of subpath i+1 = end of subpath i).
    if n >= 2
        for i in 2:n
            sj = offs[i]
            (sj < s1f || sj > s2f) && continue
            p = PathGeometry.position(path, sj)
            push!(seg_bx, p[1]); push!(seg_by, p[2]); push!(seg_bz, p[3])
            push!(seg_bound_hover,
                  "Subpath $(i-1) → $i boundary<br>s = $(sj) m<br>x, y, z = $(p[1]), $(p[2]), $(p[3])")
        end
    end
    return (; seg_bx, seg_by, seg_bz, seg_bound_hover)
end

function _collect_segment_labels(path::PathGeometry.SubpathBuilt,
                                 s1f::Float64, s2f::Float64, nudge::Float64)
    label_x = Float64[]; label_y = Float64[]; label_z = Float64[]
    label_strs = String[]
    placed = path.placed_segments
    for ps in vcat(placed, PathGeometry.PlacedSegment[path.jumpto_placed])
        nick = PathGeometry.segment_nickname(ps.segment)
        isnothing(nick) && continue
        s_lo = Float64(PathGeometry._qc_nominalize(ps.s_offset_eff))
        s_hi = s_lo + Float64(PathGeometry._qc_nominalize(
            PathGeometry.arc_length(ps.segment)))
        s_a = max(s_lo, s1f)
        s_b = min(s_hi, s2f)
        s_a >= s_b - 1e-15 && continue
        s_mid = (s_a + s_b) / 2
        fr = PathGeometry.frame(path, s_mid)
        r = collect(fr.position); N = collect(fr.normal)
        nn = norm(N)
        if nn >= 1e-12
            N ./= nn
            r .+= nudge .* N
        end
        push!(label_x, r[1]); push!(label_y, r[2]); push!(label_z, r[3])
        push!(label_strs, nick)
    end
    return (; label_x, label_y, label_z, label_strs)
end

function _collect_segment_labels(path::PathGeometry.PathBuilt,
                                 s1f::Float64, s2f::Float64, nudge::Float64)
    label_x = Float64[]; label_y = Float64[]; label_z = Float64[]
    label_strs = String[]
    offs = PathGeometry.s_offsets(path)
    for (i, sp) in enumerate(path.subpaths)
        sp_off = offs[i]
        for ps in vcat(sp.placed_segments,
                       PathGeometry.PlacedSegment[sp.jumpto_placed])
            nick = PathGeometry.segment_nickname(ps.segment)
            isnothing(nick) && continue
            s_lo = sp_off +
                Float64(PathGeometry._qc_nominalize(ps.s_offset_eff))
            s_hi = s_lo + Float64(PathGeometry._qc_nominalize(
                PathGeometry.arc_length(ps.segment)))
            s_a = max(s_lo, s1f)
            s_b = min(s_hi, s2f)
            s_a >= s_b - 1e-15 && continue
            s_mid = (s_a + s_b) / 2
            fr = PathGeometry.frame(path, s_mid)
            r = collect(fr.position); N = collect(fr.normal)
            nn = norm(N)
            if nn >= 1e-12
                N ./= nn
                r .+= nudge .* N
            end
            push!(label_x, r[1]); push!(label_y, r[2]); push!(label_z, r[3])
            push!(label_strs, nick)
        end
    end
    return (; label_x, label_y, label_z, label_strs)
end

"""
    write_path_geometry_plot3d(path, s1, s2; fidelity = 3.0, output, title, kwargs...)

Write a standalone interactive Plotly HTML file visualizing `path` (a
`PathGeometry.SubpathBuilt` or `PathGeometry.PathBuilt`) over the arc-length interval
`[s1, s2]`, and return the output path.

Horizontal mouse position (with no button pressed) scrubs arc length: left-to-right maps
linearly in `s` over the plotted `[s1, s2]` interval (not linearly in sample index). The
cursor snaps to the nearest polyline sample; the transverse normal–binormal plane and the
T/N/B frame axes update accordingly. The scene uses Plotly `aspectmode: "cube"` so one
meter maps to the same on-screen length on every axis, and the axis ranges are fixed at
generation time (from `PathGeometry.bounding_box`, padded) so scrubbing does not re-fit the
view.

Open-circle markers mark effective arc-length joins between authored segments that fall
within `[s1, s2]`. Each authored segment carrying a `Nickname` meta gets a 3D text label at
its midpoint arc length, nudged along the principal normal. A red arrow in the local N̂–B̂
plane points along `cos(Φ) N̂ + sin(Φ) B̂`, where `Φ` is `PathGeometry.total_spin` from
`s1` to the cursor arc length.

# Arguments

- `path`: a built `PathGeometry.SubpathBuilt` or `PathGeometry.PathBuilt`.
- `s1`, `s2`: arc-length bounds (m) of the plotted interval.

# Keywords

- `fidelity`: sampling-density multiplier passed to `PathGeometry.sample_path` (default 3.0).
- `output`: output HTML file path (default `"path_geometry_3d.html"`).
- `title`: plot title.
- `plane_extent_frac`, `axis_extent_frac`: half-width of the transverse square and the
  length of the T/N/B segments, as fractions of the bounding-box diagonal.
- `segment_label_nudge_frac`: normal offset of nickname labels, as a fraction of the
  bounding-box diagonal.
- `spin_n_quad`: quadrature point count passed to `total_spin` for the Φ overlay
  (default 128).

# Examples
```julia
using Bifrost.Plots

sb = PathGeometry.SubpathBuilder()
PathGeometry.start!(sb)
PathGeometry.straight!(sb; length = 0.2)
PathGeometry.bend!(sb; radius = 0.4, angle = π / 2)
PathGeometry.jumpto!(sb; point = (0.0, 0.0, 0.6))
b = PathGeometry.build(sb)
write_path_geometry_plot3d(b, 0.0, PathGeometry.s_end(b); title = "Demo",
                           fidelity = 1.0, output = "path.html")
```
"""
function write_path_geometry_plot3d(
    path::Union{PathGeometry.SubpathBuilt, PathGeometry.PathBuilt},
    s1::Real,
    s2::Real;
    fidelity::Float64 = 3.0,
    output::AbstractString = "path_geometry_3d.html",
    title::AbstractString = "Path geometry",
    plane_extent_frac::Float64 = 0.08,
    axis_extent_frac::Float64 = 0.06,
    segment_label_nudge_frac::Float64 = 0.035,
    spin_n_quad::Int = 128,
)
    path_sample = PathGeometry.sample_path(path, s1, s2; fidelity = fidelity)
    samples = _expand(path_sample)
    xs = samples.x
    ys = samples.y
    zs = samples.z

    bb = PathGeometry.bounding_box(path; n = max(path_sample.n, 512))
    lo = Vector{Float64}(bb.lo)
    hi = Vector{Float64}(bb.hi)
    lo = min.(lo, [minimum(xs), minimum(ys), minimum(zs)])
    hi = max.(hi, [maximum(xs), maximum(ys), maximum(zs)])
    diag = norm(hi - lo)
    diag = diag > 0 ? diag : 1.0
    plane_half = Float64(plane_extent_frac) * diag
    axis_len = max(Float64(axis_extent_frac) * diag, 1e-9)
    pad = sqrt(2) * plane_half + axis_len
    lo_plot = lo .- pad
    hi_plot = hi .+ pad
    x_range_js = _js_array([lo_plot[1], hi_plot[1]])
    y_range_js = _js_array([lo_plot[2], hi_plot[2]])
    z_range_js = _js_array([lo_plot[3], hi_plot[3]])

    s1f = Float64(s1)
    s2f = Float64(s2)
    nudge = Float64(segment_label_nudge_frac) * diag
    bnd = _collect_segment_boundaries(path, s1f, s2f)
    seg_bx = bnd.seg_bx; seg_by = bnd.seg_by; seg_bz = bnd.seg_bz
    seg_bound_hover = bnd.seg_bound_hover
    lbl = _collect_segment_labels(path, s1f, s2f, nudge)
    label_x = lbl.label_x; label_y = lbl.label_y; label_z = lbl.label_z
    label_strs = lbl.label_strs

    title_html = replace(replace(title, "&" => "&amp;"), "<" => "&lt;")

    s_samples = Vector{Float64}(samples.s)
    # Φ(s) = ∫_{s1}^{s} τ_spin(s') ds' for the red ∫τ_spin overlay arrow.
    integrated_tau_spin = [
        s <= s1f ? 0.0 :
            Float64(PathGeometry.total_spin(
                path; s_start = s1f, s_end = s, rtol = 1e-6))
        for s in s_samples
    ]

    html = """
    <!--
      Main 3D plot (legend names match Plotly traces):
      - path: centerline polyline; faint gray line; small markers along the curve colored by
        arc length s (Turbo colorscale).
      - segment joins: open circles at effective arc-length boundaries between authored
        segments (within the plotted s interval).
      - start / end: filled markers at the first and last sample points of the plotted
        interval.
      - cursor: filled marker at the scrub position along the path.
      - normal–binormal plane: semi-transparent square in the local plane spanned by N̂ and B̂
        at the cursor (moves with scrub).
      - T̂: orange segment, unit tangent at the cursor.
      - N̂: blue segment, principal normal at the cursor.
      - B̂: green segment, binormal T̂×N̂ at the cursor.
      - ∫τ_spin: red arrow in the N̂–B̂ plane at the cursor; Φ = total_spin(path; s_start
        = plot start, s_end = cursor) (same length scale as T̂/N̂/B̂ axes).
      - segment labels (optional): 3D text for each authored segment that has a nickname, when
        that segment overlaps the plotted s-interval.
    -->
    <!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="utf-8" />
      <meta name="viewport" content="width=device-width, initial-scale=1" />
      <title>$title_html</title>
      <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
      <style>
        html, body {
          width: 100%;
          height: 100%;
          margin: 0;
          font-family: sans-serif;
        }
        #viewer {
          position: relative;
          width: 100%;
          height: 100%;
          overflow: hidden;
        }
        #plot {
          width: 100%;
          height: 100%;
        }
        #status {
          position: absolute;
          top: 16px;
          right: 16px;
          width: min(28vw, 340px);
          min-width: 220px;
          padding: 10px 12px;
          box-sizing: border-box;
          border: 1px solid rgba(0, 0, 0, 0.15);
          background: rgba(255, 255, 255, 0.92);
          box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
          font-size: 13px;
          line-height: 1.4;
          white-space: pre-line;
        }
        #scrub-help {
          position: absolute;
          left: 16px;
          bottom: 18px;
          padding: 8px 10px;
          font-size: 13px;
          background: rgba(255, 255, 255, 0.86);
          border: 1px solid rgba(0, 0, 0, 0.12);
        }
      </style>
    </head>
    <body>
      <div id="viewer">
        <div id="plot"></div>
        <div id="status"></div>
        <div id="scrub-help">Move the mouse left-to-right over the main view to move the plane along the path (arc length).</div>
      </div>
      <script>
        const xs = $(_js_array(samples.x));
        const ys = $(_js_array(samples.y));
        const zs = $(_js_array(samples.z));
        const ss = $(_js_array(samples.s));
        const tx = $(_js_array(samples.tx));
        const ty = $(_js_array(samples.ty));
        const tz = $(_js_array(samples.tz));
        const nx = $(_js_array(samples.nx));
        const ny = $(_js_array(samples.ny));
        const nz = $(_js_array(samples.nz));
        const bx = $(_js_array(samples.bx));
        const by = $(_js_array(samples.by));
        const bz = $(_js_array(samples.bz));
        const kappa = $(_js_array(samples.kappa));
        const tauGeom = $(_js_array(samples.tau_geom));
        const tauSpin = $(_js_array(samples.tau_spin));
        const integratedTau = $(_js_array(integrated_tau_spin));
        const planeHalf = $(_js_real(plane_half));
        const axisLen = $(_js_real(axis_len));

        const xRange = $x_range_js;
        const yRange = $y_range_js;
        const zRange = $z_range_js;

        const segBoundX = $(_js_array(seg_bx));
        const segBoundY = $(_js_array(seg_by));
        const segBoundZ = $(_js_array(seg_bz));
        const segBoundHover = $(_js_string_array(seg_bound_hover));

        const labelX = $(_js_array(label_x));
        const labelY = $(_js_array(label_y));
        const labelZ = $(_js_array(label_z));
        const labelTexts = $(_js_string_array(label_strs));

        const sceneAxisLayout = {
          "scene.aspectmode": "cube",
          "scene.xaxis.autorange": false,
          "scene.xaxis.range": xRange,
          "scene.yaxis.autorange": false,
          "scene.yaxis.range": yRange,
          "scene.zaxis.autorange": false,
          "scene.zaxis.range": zRange
        };

        /** Map horizontal fraction t in [0,1] to nearest sample index so s = ss[0] + t*(ss[end]-ss[0]). */
        function scrubIndexFromT(t) {
          const n = ss.length;
          if (n <= 1) return 0;
          const s0 = ss[0];
          const sEnd = ss[n - 1];
          const span = sEnd - s0;
          if (!(span > 0)) return 0;
          const sTarget = s0 + t * span;
          if (sTarget <= s0) return 0;
          if (sTarget >= sEnd) return n - 1;
          let lo = 0;
          let hi = n - 1;
          while (lo < hi - 1) {
            const mid = (lo + hi) >> 1;
            if (ss[mid] <= sTarget) lo = mid;
            else hi = mid;
          }
          return (sTarget - ss[lo]) <= (ss[hi] - sTarget) ? lo : hi;
        }

        function planeMesh(i) {
          const r = [xs[i], ys[i], zs[i]];
          const N = [nx[i], ny[i], nz[i]];
          const B = [bx[i], by[i], bz[i]];
          const h = planeHalf;
          const corners = [
            [r[0] + h * ( N[0] + B[0]), r[1] + h * ( N[1] + B[1]), r[2] + h * ( N[2] + B[2])],
            [r[0] + h * (-N[0] + B[0]), r[1] + h * (-N[1] + B[1]), r[2] + h * (-N[2] + B[2])],
            [r[0] + h * (-N[0] - B[0]), r[1] + h * (-N[1] - B[1]), r[2] + h * (-N[2] - B[2])],
            [r[0] + h * ( N[0] - B[0]), r[1] + h * ( N[1] - B[1]), r[2] + h * ( N[2] - B[2])]
          ];
          const x = corners.map(c => c[0]);
          const y = corners.map(c => c[1]);
          const z = corners.map(c => c[2]);
          const iTri = [0, 0];
          const jTri = [1, 2];
          const kTri = [2, 3];
          return { x, y, z, i: iTri, j: jTri, k: kTri };
        }

        function axisSegments(i) {
          const r = [xs[i], ys[i], zs[i]];
          const T = [tx[i], ty[i], tz[i]];
          const N = [nx[i], ny[i], nz[i]];
          const B = [bx[i], by[i], bz[i]];
          const L = axisLen;
          const seg = (v) => [
            [r[0], r[0] + L * v[0]],
            [r[1], r[1] + L * v[1]],
            [r[2], r[2] + L * v[2]]
          ];
          const Ta = seg(T);
          const Na = seg(N);
          const Ba = seg(B);
          return {
            T: { x: Ta[0], y: Ta[1], z: Ta[2] },
            N: { x: Na[0], y: Na[1], z: Na[2] },
            B: { x: Ba[0], y: Ba[1], z: Ba[2] }
          };
        }

        function spinArrowInNBPlane(i) {
          const phi = integratedTau[i];
          const r = [xs[i], ys[i], zs[i]];
          const N = [nx[i], ny[i], nz[i]];
          const B = [bx[i], by[i], bz[i]];
          const ux = Math.cos(phi) * N[0] + Math.sin(phi) * B[0];
          const uy = Math.cos(phi) * N[1] + Math.sin(phi) * B[1];
          const uz = Math.cos(phi) * N[2] + Math.sin(phi) * B[2];
          const L = axisLen;
          return {
            x: [r[0], r[0] + L * ux],
            y: [r[1], r[1] + L * uy],
            z: [r[2], r[2] + L * uz]
          };
        }

        const pathTrace = {
          type: "scatter3d",
          mode: "lines+markers",
          x: xs,
          y: ys,
          z: zs,
          hoverinfo: "skip",
          line: {
            width: 4,
            color: "rgba(30, 30, 30, 0.25)"
          },
          marker: {
            size: 3.5,
            color: ss,
            colorscale: "Turbo",
            cmin: ss[0],
            cmax: ss[ss.length - 1],
            showscale: false
          },
          name: "path"
        };

        const segmentBoundaryTrace = {
          type: "scatter3d",
          mode: "markers",
          x: segBoundX,
          y: segBoundY,
          z: segBoundZ,
          text: segBoundHover,
          hoverinfo: "text",
          marker: {
            size: 8,
            symbol: "circle-open",
            color: "#7f7f7f",
            line: { width: 2, color: "#444444" }
          },
          name: "segment joins",
          visible: segBoundX.length > 0,
          showlegend: true
        };

        const startTrace = {
          type: "scatter3d",
          mode: "markers",
          x: [$(samples.x[1])],
          y: [$(samples.y[1])],
          z: [$(samples.z[1])],
          hovertemplate: "Start<br>x=%{x} m<br>y=%{y} m<br>z=%{z} m<extra></extra>",
          marker: { size: 7, color: "#2ca02c", symbol: "circle" },
          name: "start"
        };

        const endTrace = {
          type: "scatter3d",
          mode: "markers",
          x: [$(samples.x[end])],
          y: [$(samples.y[end])],
          z: [$(samples.z[end])],
          hovertemplate: "End<br>x=%{x} m<br>y=%{y} m<br>z=%{z} m<extra></extra>",
          marker: { size: 7, color: "#d62728", symbol: "circle" },
          name: "end"
        };

        const cursorTrace = {
          type: "scatter3d",
          mode: "markers",
          x: [xs[0]],
          y: [ys[0]],
          z: [zs[0]],
          hoverinfo: "skip",
          marker: { size: 6, color: "#111111", symbol: "circle" },
          name: "cursor"
        };

        const pm0 = planeMesh(0);
        const ax0 = axisSegments(0);

        const planeTrace = {
          type: "mesh3d",
          x: pm0.x,
          y: pm0.y,
          z: pm0.z,
          i: pm0.i,
          j: pm0.j,
          k: pm0.k,
          opacity: 0.38,
          color: "#6ba3ff",
          flatshading: true,
          hoverinfo: "skip",
          name: "normal–binormal plane"
        };

        const traceT = {
          type: "scatter3d",
          mode: "lines",
          x: ax0.T.x,
          y: ax0.T.y,
          z: ax0.T.z,
          line: { width: 5, color: "#ff7f0e" },
          hoverinfo: "skip",
          name: "T̂"
        };
        const traceN = {
          type: "scatter3d",
          mode: "lines",
          x: ax0.N.x,
          y: ax0.N.y,
          z: ax0.N.z,
          line: { width: 5, color: "#1f77b4" },
          hoverinfo: "skip",
          name: "N̂"
        };
        const traceB = {
          type: "scatter3d",
          mode: "lines",
          x: ax0.B.x,
          y: ax0.B.y,
          z: ax0.B.z,
          line: { width: 5, color: "#2ca02c" },
          hoverinfo: "skip",
          name: "B̂"
        };

        const tw0 = spinArrowInNBPlane(0);
        const traceTauInt = {
          type: "scatter3d",
          mode: "lines+markers",
          x: tw0.x,
          y: tw0.y,
          z: tw0.z,
          line: { width: 5, color: "#e31a1c" },
          marker: {
            size: [2, 9],
            color: "#e31a1c",
            symbol: ["circle", "diamond"],
            line: { width: 1, color: "#7f0000" }
          },
          hoverinfo: "skip",
          name: "∫τ_spin ds"
        };

        const layout = {
          title: $(repr(title)),
          uirevision: "path-geometry-plot",
          showlegend: true,
          legend: { x: 0.02, y: 0.98 },
          margin: { l: 0, r: 0, b: 0, t: 48 },
          scene: {
            xaxis: { title: "x (m)", range: xRange, autorange: false, showspikes: false },
            yaxis: { title: "y (m)", range: yRange, autorange: false, showspikes: false },
            zaxis: { title: "z (m)", range: zRange, autorange: false, showspikes: false },
            camera: { eye: { x: 1.6, y: 1.6, z: 0.9 } },
            aspectmode: "cube"
          }
        };

        const plotTraces = [pathTrace, segmentBoundaryTrace, startTrace, endTrace, cursorTrace, planeTrace, traceT, traceN, traceB, traceTauInt];
        if (labelX.length > 0) {
          plotTraces.push({
            type: "scatter3d",
            mode: "text",
            x: labelX,
            y: labelY,
            z: labelZ,
            text: labelTexts,
            textposition: "top center",
            textfont: { size: 14, color: "rgba(35, 35, 35, 0.92)" },
            hoverinfo: "text",
            showlegend: false
          });
        }
        Plotly.newPlot("plot", plotTraces, layout, {
          responsive: true,
          scrollZoom: true,
          displaylogo: false
        });

        const statusBox = document.getElementById("status");
        let activeIndex = 0;

        function formatStatus(index) {
          const deg = integratedTau[index] * (180 / Math.PI);
          return [
            "Arc length s = " + ss[index].toFixed(5) + " m",
            "x, y, z = " + xs[index].toFixed(4) + ", " + ys[index].toFixed(4) + ", " + zs[index].toFixed(4) + " m",
            "κ = " + kappa[index].toExponential(4) + " 1/m",
            "τ_geom = " + tauGeom[index].toExponential(4) + " rad/m",
            "τ_spin  = " + tauSpin[index].toExponential(4) + " rad/m",
            "∫τ_spin ds = " + integratedTau[index].toFixed(5) + " rad (" + deg.toFixed(2) + "°)"
          ].join("\\n");
        }

        function updateCursor(index) {
          activeIndex = Math.max(0, Math.min(xs.length - 1, index));
          const pm = planeMesh(activeIndex);
          const ax = axisSegments(activeIndex);
          statusBox.textContent = formatStatus(activeIndex);
          const chain = (p, fn) => (p && typeof p.then === "function" ? p.then(fn) : fn());
          let p = Plotly.restyle("plot", {
            x: [[xs[activeIndex]]],
            y: [[ys[activeIndex]]],
            z: [[zs[activeIndex]]]
          }, [4]);
          p = chain(p, () => Plotly.restyle("plot", {
            x: [pm.x],
            y: [pm.y],
            z: [pm.z],
            i: [pm.i],
            j: [pm.j],
            k: [pm.k]
          }, [5]));
          p = chain(p, () => Plotly.restyle("plot", {
            x: [ax.T.x],
            y: [ax.T.y],
            z: [ax.T.z]
          }, [6]));
          p = chain(p, () => Plotly.restyle("plot", {
            x: [ax.N.x],
            y: [ax.N.y],
            z: [ax.N.z]
          }, [7]));
          p = chain(p, () => Plotly.restyle("plot", {
            x: [ax.B.x],
            y: [ax.B.y],
            z: [ax.B.z]
          }, [8]));
          p = chain(p, () => {
            const tw = spinArrowInNBPlane(activeIndex);
            return Plotly.restyle("plot", {
              x: [tw.x],
              y: [tw.y],
              z: [tw.z]
            }, [9]);
          });
          chain(p, () => Plotly.relayout("plot", sceneAxisLayout));
        }

        const viewer = document.getElementById("viewer");
        viewer.addEventListener("mousemove", event => {
          if (event.buttons !== 0) return;
          const rect = viewer.getBoundingClientRect();
          const t = Math.max(0, Math.min(1, (event.clientX - rect.left) / rect.width));
          const idx = scrubIndexFromT(t);
          if (idx !== activeIndex) updateCursor(idx);
        });

        viewer.addEventListener("touchmove", event => {
          if (event.touches.length === 0) return;
          const rect = viewer.getBoundingClientRect();
          const t = Math.max(0, Math.min(1, (event.touches[0].clientX - rect.left) / rect.width));
          const idx = scrubIndexFromT(t);
          if (idx !== activeIndex) updateCursor(idx);
        }, { passive: true });

        updateCursor(0);
      </script>
    </body>
    </html>
    """

    open(output, "w") do io
        write(io, html)
    end
    return output
end
