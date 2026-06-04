# =====================================================================
# demo-path-geometry.jl
#
# Geometry-only demos used as a human-in-the-loop visual debug tool for
# `path-geometry*.jl`. Loads only the geometry layer (no fiber, no
# cross-section, no propagation). Each demo function authors a Subpath
# (or PathBuilt) using the new SubpathBuilder API, builds it, and writes
# an interactive Plotly HTML to `julia-port/output/`.
#
# Run all demos:
#
#     include("julia-port/demo-path-geometry.jl")
#     demo_path_geometry_all()
# =====================================================================

using Bifrost
using Bifrost.Plots
using Bifrost.Plots.PlotRuntime
using Bifrost.PathGeometry: _qc_nominalize
include("demo-index-helpers.jl")

const _OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "output")

_path_html(name::AbstractString) = joinpath(_OUTPUT_DIR, name)

# Convenience wrapper: SubpathBuilder built object → plot full domain.
function _plot_full(b::PathGeometry.SubpathBuilt; output, title, fidelity)
    s_hi = Float64(PathGeometry._qc_nominalize(PathGeometry.s_end(b)))
    return write_path_geometry_plot3d(b, 0.0, s_hi;
        output = output, title = title, fidelity = fidelity)
end

function _plot_full(p::PathGeometry.PathBuilt; output, title, fidelity)
    s_hi = PathGeometry.s_end(p)
    return write_path_geometry_plot3d(p, 0.0, s_hi;
        output = output, title = title, fidelity = fidelity)
end

# =====================================================================
# 2.1 Simple multi-segment Subpath (mirrors old demo_path_geometry).
# =====================================================================

function demo_path_geometry_simple(;
    output::AbstractString = _path_html("path-geometry.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Fiber path geometry: bends, catenary",
)
    PG = PathGeometry
    sb = PG.SubpathBuilder()
    PG.start!(sb)
    PG.straight!(sb; length = 0.10, meta = [PG.Nickname("Straight")])
    PG.bend!(sb; radius = 0.05, angle = π / 2, meta = [PG.Nickname("Bend")])
    PG.straight!(sb; length = 0.12, meta = [PG.Nickname("Straight")])
    PG.catenary!(sb; a = 0.03, length = 0.10, axis_angle = 0.0,
                 meta = [PG.Nickname("Catenary")])
    PG.bend!(sb; radius = 0.06, angle = π / 3, meta = [PG.Nickname("Bend")])
    PG.straight!(sb; length = 0.08, meta = [PG.Nickname("Straight")])
    PG.seal!(sb)
    sub = PG.Subpath(sb)
    b = PG.build(sub)
    println("Arc length: ", PG.path_length(b))
    println("Writhe:     ", PG.writhe(b; n = 128))
    plot_path = _plot_full(b; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = b, plot_path)
end

# =====================================================================
# 2.2 Segment labels (mirrors old demo_path_geometry_segment_labels).
# =====================================================================

function demo_path_geometry_segment_labels(;
    output::AbstractString = _path_html("path-geometry-segment-labels.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Path geometry: segment nicknames",
)
    PG = PathGeometry
    sb = PG.SubpathBuilder()
    PG.start!(sb)
    PG.straight!(sb; length = 0.08, meta = [PG.Nickname("lead-in")])
    PG.bend!(sb; radius = 0.06, angle = π / 2, meta = [PG.Nickname("90° bend")])
    PG.straight!(sb; length = 0.06, meta = [PG.Nickname("spacer")])
    PG.catenary!(sb; a = 0.04, length = 0.08, axis_angle = 0.0,
                 meta = [PG.Nickname("sag")])
    PG.helix!(sb; radius = 0.025, pitch = 0.015, turns = 1.2,
              axis_angle = 0.0, meta = [PG.Nickname("spin section")])
    PG.straight!(sb; length = 0.06, meta = [PG.Nickname("lead-out")])
    b = PG.build(PG.Subpath(PG.seal!(sb)))
    println("Arc length: ", PG.path_length(b), " m")
    plot_path = _plot_full(b; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = b, plot_path)
end

# =====================================================================
# 2.3 Helix demos for axis_angle = 0, π/3, 2π/3.
# =====================================================================

function _demo_helix(axis_angle::Float64; output, fidelity, title)
    PG = PathGeometry
    sb = PG.SubpathBuilder()
    PG.start!(sb)
    PG.straight!(sb; length = 0.05, meta = [PG.Nickname("Straight")])
    PG.helix!(sb; radius = 0.03, pitch = 0.02, turns = 2.0,
              axis_angle = axis_angle, meta = [PG.Nickname("Helix")])
    PG.straight!(sb; length = 0.05, meta = [PG.Nickname("Straight")])
    b = PG.build(PG.Subpath(PG.seal!(sb)))
    println("Helix axis_angle=", axis_angle, ": arc_length=",
            round(PG.path_length(b); digits = 4), " m")
    plot_path = _plot_full(b; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = b, plot_path)
end

demo_path_geometry_helix_0(;
    output::AbstractString = _path_html("path-geometry-helix-0.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = 0",
) = _demo_helix(0.0; output, fidelity, title)

demo_path_geometry_helix_pi_3(;
    output::AbstractString = _path_html("path-geometry-helix-pi-3.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = π/3",
) = _demo_helix(π / 3; output, fidelity, title)

demo_path_geometry_helix_2pi_3(;
    output::AbstractString = _path_html("path-geometry-helix-2pi-3.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = 2π/3",
) = _demo_helix(2π / 3; output, fidelity, title)

# =====================================================================
# 2.4 jumps_min_radius — paddle pattern realized as a PathBuilt.
#
# Old demo: alternating straight / jumpto / straight / jumpto / ...
# In the new architecture each `jumpto!` seals a Subpath, so the
# original 4-segment 4-jump pattern becomes 5 Subpaths.
# =====================================================================

function demo_path_geometry_jumps_min_radius(;
    output::AbstractString = _path_html("path-geometry-jumps-min-radius.html"),
    fidelity::Float64 = 4.0,
    title::AbstractString = "JumpBy/JumpTo paddle: PathBuilt of 5 Subpaths",
)
    PG = PathGeometry
    # Subpath 1: straight up to (0,0,1), seal with jumpto landing at (1,0,1)
    # with incoming tangent (0,0,-1) (heading down at landing). This gives
    # a transverse chord; min_bend_radius=0.4 keeps the connector smooth.
    sb1 = PG.SubpathBuilder()
    PG.start!(sb1)
    PG.straight!(sb1; length = 1.0, meta = [PG.Nickname("Sub1 straight")])
    PG.jumpto!(sb1; point = (1.0, 0.0, 1.0),
               incoming_tangent = (0.0, 0.0, -1.0),
               min_bend_radius = 0.4)

    # Subpath 2: starts at (1,0,1) heading -z, straight to (1,0,0), seals
    # to (2,0,0) with incoming tangent (0,0,1).
    sb2 = PG.SubpathBuilder()
    PG.start!(sb2; point = (1.0, 0.0, 1.0),
                  outgoing_tangent = (0.0, 0.0, -1.0))
    PG.straight!(sb2; length = 1.0, meta = [PG.Nickname("Sub2 straight")])
    PG.jumpto!(sb2; point = (2.0, 0.0, 0.0),
               incoming_tangent = (0.0, 0.0, 1.0),
               min_bend_radius = 0.1)

    # Subpath 3: starts at (2,0,0) heading +z, straight to (2,0,1),
    # seals to (3,0,1) with incoming tangent (0,0,-1).
    sb3 = PG.SubpathBuilder()
    PG.start!(sb3; point = (2.0, 0.0, 0.0),
                  outgoing_tangent = (0.0, 0.0, 1.0))
    PG.straight!(sb3; length = 1.0, meta = [PG.Nickname("Sub3 straight")])
    PG.jumpto!(sb3; point = (3.0, 0.0, 1.0),
               incoming_tangent = (0.0, 0.0, -1.0),
               min_bend_radius = 0.05)

    # Subpath 4: starts at (3,0,1) heading -z, straight + interior JumpBy.
    # JumpBy delta is in the local frame; after the straight, local +z is
    # global -z, so delta=(-1,0,0) (local) lands at (2,0,0) heading +z.
    # Seal with an explicit jumpto! at that landing point so the endpoint is
    # declared up front — Subpath 5 starts there directly, no probe build.
    sb4 = PG.SubpathBuilder()
    PG.start!(sb4; point = (3.0, 0.0, 1.0),
                  outgoing_tangent = (0.0, 0.0, -1.0))
    PG.straight!(sb4; length = 1.0, meta = [PG.Nickname("Sub4 straight")])
    PG.jumpby!(sb4; delta = (-1.0, 0.0, 0.0),
               tangent = (0.0, 0.0, -1.0),
               min_bend_radius = 0.1,
               meta = [PG.Nickname("Sub4 JumpBy")])
    PG.jumpto!(sb4; point = (2.0, 0.0, 0.0),
               incoming_tangent = (0.0, 0.0, 1.0))

    # Subpath 5: continues straight from sb4's landing point (2,0,0), +z.
    sb5 = PG.SubpathBuilder()
    PG.start!(sb5; point = (2.0, 0.0, 0.0),
                  outgoing_tangent = (0.0, 0.0, 1.0))
    PG.straight!(sb5; length = 1.0, meta = [PG.Nickname("Sub5 straight")])
    sb5 = PG.seal!(sb5)

    p = PG.build([PG.Subpath(sb1), PG.Subpath(sb2), PG.Subpath(sb3),
                  PG.Subpath(sb4), PG.Subpath(sb5)])
    println("PathBuilt arc length: ", PG.path_length(p), " m")
    plot_path = _plot_full(p; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = p, plot_path)
end

# =====================================================================
# 2.5 NEW: Multi-Subpath demo showing PathBuilt assembly.
# =====================================================================

function demo_path_geometry_pathbuilt(;
    output::AbstractString = _path_html("path-geometry-pathbuilt.html"),
    fidelity::Float64 = 2.0,
    title::AbstractString = "PathBuilt: three Subpaths (straight, bend, helix)",
)
    PG = PathGeometry

    # Subpath 1: straight, sealed at (0,0,0.2) with tangent +z.
    sb1 = PG.SubpathBuilder(meta = [PG.Nickname("Subpath 1: straight")])
    PG.start!(sb1)
    PG.straight!(sb1; length = 0.2, meta = [PG.Nickname("Straight")])
    PG.jumpto!(sb1; point = (0.0, 0.0, 0.2),
               incoming_tangent = (0.0, 0.0, 1.0))

    # Subpath 2: starts at (0,0,0.2) tangent +z, quarter bend (axis_angle=0),
    # so end position (R, 0, 0.2 + R) and end tangent +x.
    R2 = 0.05
    sb2 = PG.SubpathBuilder(meta = [PG.Nickname("Subpath 2: bend")])
    PG.start!(sb2; point = (0.0, 0.0, 0.2),
                  outgoing_tangent = (0.0, 0.0, 1.0))
    PG.bend!(sb2; radius = R2, angle = π / 2,
             meta = [PG.Nickname("90° bend")])
    PG.jumpto!(sb2; point = (R2, 0.0, 0.2 + R2),
               incoming_tangent = (1.0, 0.0, 0.0))

    # Subpath 3: starts at (R2, 0, 0.2 + R2) tangent +x, helix in transverse
    # plane (axis_angle=0). Seal at the helix's natural exit.
    sb3 = PG.SubpathBuilder(meta = [PG.Nickname("Subpath 3: helix")])
    PG.start!(sb3; point = (R2, 0.0, 0.2 + R2),
                  outgoing_tangent = (1.0, 0.0, 0.0))
    PG.helix!(sb3; radius = 0.025, pitch = 0.02, turns = 1.5,
              axis_angle = 0.0, meta = [PG.Nickname("Helix")])
    sb3 = PG.seal!(sb3)

    p = PG.build([PG.Subpath(sb1), PG.Subpath(sb2), PG.Subpath(sb3)])
    println("PathBuilt: ", length(p.subpaths), " subpaths")
    println("PathBuilt arc length: ", PG.path_length(p), " m")
    plot_path = _plot_full(p; output, title, fidelity)
    println("Wrote ", plot_path)
    return (; path = p, plot_path)
end

# =====================================================================
# 2.6 Aggregator + index HTML
# =====================================================================

# Registry of demos. Each entry: (group, fn, desc). Used by
# `demo_path_geometry_all` to run every demo and emit an index HTML that
# links to each output with a short description.
const _DEMO_PATH_GEOMETRY_INDEX = [
    (group = "Subpath",
     fn    = demo_path_geometry_simple,
     desc  = "Single Subpath: straight + bend + straight + catenary + bend + straight, sealed at the natural exit."),
    (group = "Subpath",
     fn    = demo_path_geometry_segment_labels,
     desc  = "Same shape as 'simple' but every segment carries a Nickname so labels render in the plot."),
    (group = "Subpath — helix",
     fn    = demo_path_geometry_helix_0,
     desc  = "HelixSegment with axis_angle = 0."),
    (group = "Subpath — helix",
     fn    = demo_path_geometry_helix_pi_3,
     desc  = "HelixSegment with axis_angle = π/3."),
    (group = "Subpath — helix",
     fn    = demo_path_geometry_helix_2pi_3,
     desc  = "HelixSegment with axis_angle = 2π/3."),
    (group = "PathBuilt",
     fn    = demo_path_geometry_jumps_min_radius,
     desc  = "Paddle pattern: 5 Subpaths joined at jumpto endpoints, each with its own min_bend_radius. Includes one interior JumpBy segment."),
    (group = "PathBuilt",
     fn    = demo_path_geometry_pathbuilt,
     desc  = "Three Subpaths (straight, quarter-bend, helix) stitched into a PathBuilt. Demonstrates the conformity check at each Subpath boundary."),
]

function demo_path_geometry_all(;
    index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT,
)
    isdir(_OUTPUT_DIR) || mkpath(_OUTPUT_DIR)

    return _write_demo_index(
        [(title = "Path-geometry-only demos",
          entries = demo_path_geometry_entries(),
          group_titles = Dict{String, String}())];
        index_output,
    )
end

function demo_path_geometry_entries()
    entries = Tuple{String, String, String, String}[]
    for d in _DEMO_PATH_GEOMETRY_INDEX
        println("[ demo ] $(nameof(d.fn))")
        result = d.fn()
        push!(entries, (d.group, basename(result.plot_path),
                        result.plot_path, d.desc))
    end
    return entries
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_path_geometry_all()
end
