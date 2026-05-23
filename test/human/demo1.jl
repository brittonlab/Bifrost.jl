using LinearAlgebra
using MonteCarloMeasurements
using Bifrost
using Bifrost.Plots
using Bifrost.Plots.PlotRuntime


const DEMO_FIBER_CROSS_SECTION = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like demo"
)
const DEMO_λ_M = 1550e-9
const DEMO_T_K = 297.15

# =====================================================================
# Path-geometry demos
# =====================================================================

function demo_path_geometry(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "path-geometry.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Fiber path geometry: bends, catenary, twist",
)
    PG = PathGeometry
    spec = PG.PathSpecBuilder()
    L_lead = 0.10
    R1 = 0.05
    θ1 = π / 2
    L_spacer = 0.12
    PG.straight!(spec; length = L_lead, meta = [PG.Nickname("Straight")])
    PG.bend!(spec; radius = R1, angle = θ1, meta = [PG.Nickname("Bend")])
    PG.straight!(spec; length = L_spacer, meta = [PG.Nickname("Straight")])
    PG.catenary!(spec; a = 0.03, length = 0.10, axis_angle = 0.0, meta = [PG.Nickname("Catenary")])
    PG.bend!(spec; radius = 0.06, angle = π / 3, meta = [PG.Nickname("Bend")])
    PG.straight!(spec; length = 0.08, meta = [PG.Nickname("Straight")])
    # TODO: twist refactor — PG.twist!(spec; s_start = s_twist, length = L_spacer, rate = 35.0)
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path))
    println("Writhe: ", PG.writhe(path))
    plot_path = write_path_geometry_plot3d(
        path,
        path.spec.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    println("Wrote path geometry plot to: ", plot_path)
    return (; path, plot_path)
end

function demo_path_geometry_segment_labels(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "path-geometry-segment-labels.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Path geometry: segment nicknames",
)
    PG = PathGeometry
    spec = PG.PathSpecBuilder()
    PG.straight!(spec; length = 0.08, meta = [PG.Nickname("lead-in")])
    PG.bend!(spec; radius = 0.06, angle = π / 2, meta = [PG.Nickname("90° bend")])
    PG.straight!(spec; length = 0.06, meta = [PG.Nickname("spacer")])
    PG.catenary!(spec; a = 0.04, length = 0.08, axis_angle = 0.0, meta = [PG.Nickname("sag")])
    PG.helix!(spec; radius = 0.025, pitch = 0.015, turns = 1.2, axis_angle = 0.0, meta = [PG.Nickname("twist section")])
    PG.straight!(spec; length = 0.06, meta = [PG.Nickname("lead-out")])
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path), " m")
    plot_path = write_path_geometry_plot3d(
        path,
        path.spec.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    println("Wrote segment-label demo to: ", plot_path)
    return (; path, plot_path)
end

function demo_path_geometry_helix_0(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "path-geometry-helix-0.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = 0",
)
    PG = PathGeometry
    spec = PG.PathSpecBuilder()
    PG.straight!(spec; length = 0.05, meta = [PG.Nickname("Straight")])
    PG.helix!(spec; radius = 0.03, pitch = 0.02, turns = 2.0, axis_angle = 0.0, meta = [PG.Nickname("Helix")])
    PG.straight!(spec; length = 0.05, meta = [PG.Nickname("Straight")])
    path = PG.build(spec)
    println("Helix axis_angle=0: arc_length=$(round(PG.path_length(path), digits=4)) m")
    plot_path = write_path_geometry_plot3d(
        path, path.spec.s_start, path.s_end;
        fidelity = fidelity, output = output, title = title,
    )
    println("Wrote helix demo to: ", plot_path)
    return (; path, plot_path)
end

function demo_path_geometry_helix_pi_3(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "path-geometry-helix-pi-3.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = π/3",
)
    PG = PathGeometry
    spec = PG.PathSpecBuilder()
    PG.straight!(spec; length = 0.05, meta = [PG.Nickname("Straight")])
    PG.helix!(spec; radius = 0.03, pitch = 0.02, turns = 2.0, axis_angle = π/3, meta = [PG.Nickname("Helix")])
    PG.straight!(spec; length = 0.05, meta = [PG.Nickname("Straight")])
    path = PG.build(spec)
    println("Helix axis_angle=π/3: arc_length=$(round(PG.path_length(path), digits=4)) m")
    plot_path = write_path_geometry_plot3d(
        path, path.spec.s_start, path.s_end;
        fidelity = fidelity, output = output, title = title,
    )
    println("Wrote helix demo to: ", plot_path)
    return (; path, plot_path)
end

function demo_path_geometry_helix_2pi_3(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "path-geometry-helix-2pi-3.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: axis_angle = 2π/3",
)
    PG = PathGeometry
    spec = PG.PathSpecBuilder()
    PG.straight!(spec; length = 0.05, meta = [PG.Nickname("Straight")])
    PG.helix!(spec; radius = 0.03, pitch = 0.02, turns = 2.0, axis_angle = 2π/3, meta = [PG.Nickname("Helix")])
    PG.straight!(spec; length = 0.05, meta = [PG.Nickname("Straight")])
    path = PG.build(spec)
    println("Helix axis_angle=2π/3: arc_length=$(round(PG.path_length(path), digits=4)) m")
    plot_path = write_path_geometry_plot3d(
        path, path.spec.s_start, path.s_end;
        fidelity = fidelity, output = output, title = title,
    )
    println("Wrote helix demo to: ", plot_path)
    return (; path, plot_path)
end

function demo_path_geometry_jumps_min_radius(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "path-geometry-jumps-min-radius.html"),
    fidelity::Float64 = 4.0,
    title::AbstractString = "JumpBy and JumpTo: Hermite connectors, min_bend_radius",
)
    PG = PathGeometry
    spec = PG.PathSpecBuilder()

    PG.straight!(spec; length = 1, meta = [PG.Nickname("Straight")])
    PG.jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.4, meta = [PG.Nickname("JumpTo")]) # T1: fails only if min_bend_radius >0.5
    PG.straight!(spec; length = 1, meta = [PG.Nickname("Straight")])
    PG.jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
            min_bend_radius = 0.1, meta = [PG.Nickname("JumpTo")]) # T2: fails only if min_bend_radius >0.5
    PG.straight!(spec; length = 1, meta = [PG.Nickname("Straight")])
    PG.jumpto!(spec; destination = (3, 0.0, 1), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.05, meta = [PG.Nickname("JumpTo")]) # T3: fails only if min_bend_radius >0.50
    PG.straight!(spec; length = 1, meta = [PG.Nickname("Straight")])
    PG.jumpby!(spec; delta = (-1, 0.0, 0), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.1, meta = [PG.Nickname("JumpBy")])
    PG.straight!(spec; length = 1, meta = [PG.Nickname("Straight")])

    path = PG.build(spec)
    plot_path = write_path_geometry_plot3d(
        path, path.spec.s_start, path.s_end;
        fidelity = fidelity, output = output, title = title)
    return (; path, plot_path)
end

# =====================================================================
# Modify demos — helpers shared by all 12 rows
# =====================================================================

# Sample a single placed segment's centerline in the global frame.
function _sample_segment_xyz(path::PathSpecCached, seg_index::Int; n::Int = 128)
    ps = path.placed_segments[seg_index]
    s0 = ps.s_offset_eff
    s1 = s0 + arc_length(ps.segment)
    ss = range(s0, s1; length = n)
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for s in ss
        p = position(path, s)
        push!(xs, Float64(p[1]))
        push!(ys, Float64(p[2]))
        push!(zs, Float64(p[3]))
    end
    return (x = xs, y = ys, z = zs)
end

# Build the baseline 3-segment inverted-U and attach `target_meta` to the
# segment indexed by `target_idx` (1 = first straight, 2 = bend, 3 = second
# straight). Returns the modified path and the target segment index.
function _build_modify_variant(L::Float64, R::Float64,
                               target_idx::Int, target_meta::AbstractVector{<:AbstractMeta})
    spec = PathSpecBuilder()
    straight!(spec; length = L, meta = target_idx == 1 ? target_meta : AbstractMeta[])
    bend!(spec; radius = R, angle = π, axis_angle = 0.0,
          meta = target_idx == 2 ? target_meta : AbstractMeta[])
    straight!(spec; length = L, meta = target_idx == 3 ? target_meta : AbstractMeta[])
    path  = build(spec)
    fiber = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = DEMO_T_K)
    return modify(fiber)
end

# Baseline helix parameters for the 4-segment variant (inverted-U plus a helix
# between the bend and the second straight). `axis_angle = 0` tilts the helix
# axis into the x–z plane so the coil projects recognisably onto the viewer's
# plane.
const _MODIFY_HELIX_RADIUS = 0.15
const _MODIFY_HELIX_PITCH  = 0.25
const _MODIFY_HELIX_TURNS  = 1.5

# Build a 4-segment baseline (straight · bend π · helix · straight) and attach
# `target_meta` to the segment indexed by `target_idx`:
#   1 = first straight, 2 = bend, 3 = helix, 4 = second straight.
function _build_modify_variant_helix(L::Float64, R::Float64,
                                     target_idx::Int, target_meta::AbstractVector{<:AbstractMeta})
    spec = PathSpecBuilder()
    straight!(spec; length = L, meta = target_idx == 1 ? target_meta : AbstractMeta[])
    bend!(spec; radius = R, angle = π, axis_angle = 0.0,
          meta = target_idx == 2 ? target_meta : AbstractMeta[])
    helix!(spec; radius = _MODIFY_HELIX_RADIUS,
                 pitch  = _MODIFY_HELIX_PITCH,
                 turns  = _MODIFY_HELIX_TURNS,
                 axis_angle = 0.0,
                 meta = target_idx == 3 ? target_meta : AbstractMeta[])
    straight!(spec; length = L, meta = target_idx == 4 ? target_meta : AbstractMeta[])
    path  = build(spec)
    fiber = Fiber(path; cross_section = DEMO_FIBER_CROSS_SECTION, T_ref_K = DEMO_T_K)
    return modify(fiber)
end

# Render one experiment row: each variant's path laid out side-by-side along x,
# with the target segment drawn in red and the others in green.
function _modify_row_html(output::AbstractString, title::AbstractString;
                          L::Float64, R::Float64,
                          variants::Vector,
                          variant_spacing::Float64 = 2.5,
                          builder = _build_modify_variant)
    js_num(x)   = isnan(x) ? "NaN" : string(Float64(x))
    js_arr(xs)  = "[" * join(js_num.(xs), ",") * "]"
    js_strarr(xs) = "[" * join(("\"" * replace(string(x), "\"" => "\\\"") * "\"" for x in xs), ",") * "]"

    trace_strs  = String[]
    label_xs    = Float64[]; label_ys = Float64[]; label_zs = Float64[]
    label_texts = String[]
    start_xs    = Float64[]; start_ys = Float64[]; start_zs = Float64[]

    for (k, (label, target_idx, target_meta)) in enumerate(variants)
        path = builder(L, R, target_idx, target_meta)
        dx = (k - 1) * variant_spacing
        n_segs = length(path.placed_segments)
        for i in 1:n_segs
            s     = _sample_segment_xyz(path, i)
            xs    = s.x .+ dx
            color = (i == target_idx && !isempty(target_meta)) ? "#e55" : "#6c6"
            name  = "$(label) — seg $(i)"
            push!(trace_strs, """
{
  type: 'scatter3d',
  mode: 'lines',
  x: $(js_arr(xs)),
  y: $(js_arr(s.y)),
  z: $(js_arr(s.z)),
  line: {color: '$(color)', width: 6},
  name: $(repr(name)),
  showlegend: $(i == 1)
}""")
        end
        p0 = position(path, path.spec.s_start)
        push!(start_xs, Float64(p0[1]) + dx)
        push!(start_ys, Float64(p0[2]))
        push!(start_zs, Float64(p0[3]))

        push!(label_xs, dx + R)
        push!(label_ys, 0.0)
        push!(label_zs, L + R + 0.15)
        push!(label_texts, label)
    end

    start_trace = """
{
  type: 'scatter3d',
  mode: 'markers',
  x: $(js_arr(start_xs)),
  y: $(js_arr(start_ys)),
  z: $(js_arr(start_zs)),
  marker: {color: '#6c6', size: 6, symbol: 'circle'},
  name: 'start',
  showlegend: false
}"""

    labels_trace = """
{
  type: 'scatter3d',
  mode: 'text',
  x: $(js_arr(label_xs)),
  y: $(js_arr(label_ys)),
  z: $(js_arr(label_zs)),
  text: $(js_strarr(label_texts)),
  textfont: {color: '#ddd', size: 13},
  showlegend: false
}"""

    html = """<!DOCTYPE html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><title>$(title)</title>
<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>
<style>html,body{margin:0;padding:0;width:100%;height:100%;background:#111;color:#eee;font-family:sans-serif;}#plot{width:100%;height:100%;}</style>
</head><body><div id=\"plot\"></div>
<script>
const traces = [$(join(trace_strs, ",\n")),
$(start_trace)$(isempty(label_texts) ? "" : ",\n")$(labels_trace)];
const layout = {
  title: {text: $(repr(title)), font: {color: '#eee'}},
  paper_bgcolor: '#111',
  plot_bgcolor:  '#111',
  scene: {
    bgcolor: '#111',
    xaxis: {gridcolor: '#333', color: '#aaa'},
    yaxis: {gridcolor: '#333', color: '#aaa'},
    zaxis: {gridcolor: '#333', color: '#aaa'},
    aspectmode: 'data',
    camera: {eye: {x: 0.0, y: -2.5, z: 0.5}}
  },
  margin: {l:0, r:0, t:40, b:0},
  legend: {font: {color: '#eee'}}
};
Plotly.newPlot('plot', traces, layout);
</script>
</body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end

# =====================================================================
# Modify demos — one function per HTML file
# =====================================================================

function demo_modify_straight_length(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-straight-length.html"),
        "modify(:length) on the first straight segment";
        L = 1.0, R = 0.5,
        variants = [
            ("baseline",              1, AbstractMeta[]),
            ("MCMadd(:length, -0.4)", 1, [MCMadd(:length, -0.4)]),
        ],
    )
end

function demo_modify_bend_radius(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-bend-radius.html"),
        "modify(:radius) on the bend segment";
        L = 1.0, R = 0.5,
        variants = [
            ("baseline",               2, AbstractMeta[]),
            ("MCMadd(:radius, -0.25)", 2, [MCMadd(:radius, -0.25)]),
            ("MCMadd(:radius, +0.50)", 2, [MCMadd(:radius,  0.50)]),
        ],
    )
end

function demo_modify_bend_angle(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-bend-angle.html"),
        "modify(:angle) on the bend segment";
        L = 1.0, R = 0.5,
        variants = [
            ("baseline",            2, AbstractMeta[]),
            ("MCMadd(:angle, -π/2)", 2, [MCMadd(:angle, -π/2)]),
            ("MCMadd(:angle, +π/4)", 2, [MCMadd(:angle,  π/4)]),
            ("MCMadd(:angle, +π)",   2, [MCMadd(:angle,  π)]),
        ],
    )
end

function demo_modify_straight_length_mul(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-straight-length-mul.html"),
        "modify(:length) on the first straight segment — MCMmul";
        L = 1.0, R = 0.5,
        variants = [
            ("baseline",              1, AbstractMeta[]),
            ("MCMmul(:length, -0.4)", 1, [MCMmul(:length, -0.4)]),
            ("MCMmul(:length, +0.5)", 1, [MCMmul(:length,  0.5)]),
        ],
    )
end

function demo_modify_bend_radius_mul(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-bend-radius-mul.html"),
        "modify(:radius) on the bend segment — MCMmul";
        L = 1.0, R = 0.5,
        variants = [
            ("baseline",             2, AbstractMeta[]),
            ("MCMmul(:radius, 0.5)", 2, [MCMmul(:radius, 0.5)]),
            # neg radius not supported, throws @assertion error:
            # ("MCMmul(:radius, -0.5)", 2, [MCMmul(:radius, -0.5)]),
            ("MCMmul(:radius, 2.0)", 2, [MCMmul(:radius, 2.0)]),
        ],
    )
end

function demo_modify_bend_angle_mul(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-bend-angle-mul.html"),
        "modify(:angle) on the bend segment — MCMmul";
        L = 1.0, R = 0.5,
        variants = [
            ("baseline",            2, AbstractMeta[]),
            ("MCMmul(:angle, 0.5)",  2, [MCMmul(:angle, 0.5)]),
            ("MCMmul(:angle, 1.25)", 2, [MCMmul(:angle, 1.25)]),
        ],
    )
end

function demo_modify_helix_radius(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-helix-radius.html"),
        "modify(:radius) on the helix segment — MCMadd";
        L = 1.0, R = 0.5, variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",               3, AbstractMeta[]),
            ("MCMadd(:radius, -0.05)", 3, [MCMadd(:radius, -0.05)]),
            ("MCMadd(:radius, +0.10)", 3, [MCMadd(:radius,  0.10)]),
        ],
    )
end

function demo_modify_helix_pitch(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-helix-pitch.html"),
        "modify(:pitch) on the helix segment — MCMadd";
        L = 1.0, R = 0.5, variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",              3, AbstractMeta[]),
            ("MCMadd(:pitch, -0.10)", 3, [MCMadd(:pitch, -0.10)]),
            ("MCMadd(:pitch, +0.20)", 3, [MCMadd(:pitch,  0.20)]),
        ],
    )
end

function demo_modify_helix_turns(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-helix-turns.html"),
        "modify(:turns) on the helix segment — MCMadd";
        L = 1.0, R = 0.5, variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",             3, AbstractMeta[]),
            ("MCMadd(:turns, -0.5)", 3, [MCMadd(:turns, -0.5)]),
            ("MCMadd(:turns, +0.5)", 3, [MCMadd(:turns,  0.5)]),
        ],
    )
end

function demo_modify_helix_radius_mul(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-helix-radius-mul.html"),
        "modify(:radius) on the helix segment — MCMmul";
        L = 1.0, R = 0.5, variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",             3, AbstractMeta[]),
            ("MCMmul(:radius, 0.5)", 3, [MCMmul(:radius, 0.5)]),
            ("MCMmul(:radius, 2.0)", 3, [MCMmul(:radius, 2.0)]),
        ],
    )
end

function demo_modify_helix_pitch_mul(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-helix-pitch-mul.html"),
        "modify(:pitch) on the helix segment — MCMmul";
        L = 1.0, R = 0.5, variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",            3, AbstractMeta[]),
            ("MCMmul(:pitch, 0.5)", 3, [MCMmul(:pitch, 0.5)]),
            ("MCMmul(:pitch, 2.0)", 3, [MCMmul(:pitch, 2.0)]),
        ],
    )
end

function demo_modify_helix_turns_mul(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _modify_row_html(
        joinpath(output_dir, "modify-helix-turns-mul.html"),
        "modify(:turns) on the helix segment — MCMmul";
        L = 1.0, R = 0.5, variant_spacing = 3.5,
        builder = _build_modify_variant_helix,
        variants = [
            ("baseline",             3, AbstractMeta[]),
            ("MCMmul(:turns, 0.67)", 3, [MCMmul(:turns, 0.67)]),
            ("MCMmul(:turns, 1.5)",  3, [MCMmul(:turns, 1.5)]),
        ],
    )
end

# =====================================================================
# MCM temperature and twist demos
# =====================================================================

# =====================================================================
# Adaptive step-doubling demo
# =====================================================================

"""
    demo_adaptive_step_doubling(; output, rtol, atol, title)

Illustrate the adaptive step-doubling integrator on a smooth, noncommuting generator

    K(s) = α·i·σ_x·cos(π·s) + β·i·σ_z·sin(2π·s),   s ∈ [0, 2]

The two Pauli components oscillate at different frequencies, so ‖K(s)‖ varies along the path
and the integrator must work harder near the fast oscillation peaks. The output plot shows:

- **Top panel**: accepted (green) and rejected (red) step sizes vs position, with the
  generator norm ‖K(s)‖ as a shaded overlay — small steps should cluster where K varies fastest.
- **Bottom panel**: err/tol ratio for every trial, with the acceptance threshold at 1.

Calls `collect_adaptive_steps` and `write_adaptive_steps_plot` from `fiber-path-plot.jl`;
the production solver in `path-integral.jl` is not modified.
"""
function demo_adaptive_step_doubling(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "adaptive-step-doubling.html"),
    rtol::Float64 = 1e-6,
    atol::Float64 = 1e-9,
    title::AbstractString = "Adaptive step-doubling: noncommuting K(s)"
)
    SX = ComplexF64[0 1; 1 0]
    SZ = ComplexF64[1 0; 0 -1]
    α = 1.2
    β = 0.9
    s0, s1 = 0.0, 2.0

    K = s -> α * im * cos(π * s) .* SX + β * im * sin(2π * s) .* SZ
    K_norm = s -> opnorm(K(s))
    components = [
        ("α·cos(πs)  [σx coeff]", s -> α * cos(π * s),  "#1f77b4"),
        ("β·sin(2πs) [σz coeff]", s -> β * sin(2π * s), "#ff7f0e"),
    ]

    J0 = Matrix{ComplexF64}(I, 2, 2)
    J_final, records = collect_adaptive_steps(K, s0, s1, J0; rtol = rtol, atol = atol)

    n_acc = count(r.accepted for r in records)
    n_rej = count(!r.accepted for r in records)
    println("Accepted steps: $n_acc,  rejected: $n_rej")
    println("Final Jones matrix:\n", J_final)

    plot_path = write_adaptive_steps_plot(
        records, K_norm, s0, s1;
        output = output,
        title  = title,
        rtol   = rtol,
        atol   = atol,
        components = components,
    )
    println("Wrote adaptive step-doubling plot to: ", plot_path)
    return output
end

# =====================================================================
# Poincaré rendering helper
# =====================================================================

function _write_poincare_cloud_plot(
    rep_mean::NamedTuple,
    stokes_samples::Vector,
    sample_colors::Vector{Float64};
    output::AbstractString,
    title::AbstractString,
)
    sphere = render_poincare_sphere(rep_mean)
    xs = [s[1] for s in stokes_samples]
    ys = [s[2] for s in stokes_samples]
    zs = [s[3] for s in stokes_samples]

    js_num(x) = isnan(x) ? "NaN" : string(Float64(x))
    js_arr(xs) = "[" * join(js_num.(xs), ",") * "]"
    js_surf(xss) = "[" * join([js_arr(r) for r in xss], ",") * "]"
    js_strarr(xs) = "[" * join(["\"$(x)\"" for x in xs], ",") * "]"

    html = """
<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>$title</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;padding:0;width:100%;height:100%;font-family:sans-serif;}#plot{width:100%;height:100%;}</style>
</head><body><div id="plot"></div>
<script>
const sphU = $(js_surf(sphere.surface.x));
const sphV = $(js_surf(sphere.surface.y));
const sphW = $(js_surf(sphere.surface.z));
const eqXYx = $(js_arr(sphere.circles.xy.x));
const eqXYy = $(js_arr(sphere.circles.xy.y));
const eqXYz = $(js_arr(sphere.circles.xy.z));
const eqXZx = $(js_arr(sphere.circles.xz.x));
const eqXZy = $(js_arr(sphere.circles.xz.y));
const eqXZz = $(js_arr(sphere.circles.xz.z));
const eqYZx = $(js_arr(sphere.circles.yz.x));
const eqYZy = $(js_arr(sphere.circles.yz.y));
const eqYZz = $(js_arr(sphere.circles.yz.z));
const vecX = $(js_arr(sphere.vector.x));
const vecY = $(js_arr(sphere.vector.y));
const vecZ = $(js_arr(sphere.vector.z));
const labX = $(js_arr(sphere.labels.x));
const labY = $(js_arr(sphere.labels.y));
const labZ = $(js_arr(sphere.labels.z));
const labT = $(js_strarr(sphere.labels.text));
const cloudX = $(js_arr(xs));
const cloudY = $(js_arr(ys));
const cloudZ = $(js_arr(zs));
const cloudC = $(js_arr(sample_colors));
const traces = [
  {type:"surface", x:sphU, y:sphV, z:sphW, opacity:0.12, showscale:false,
   colorscale:[[0,"#dddddd"],[1,"#dddddd"]], hoverinfo:"skip", contours:{x:{highlight:false},y:{highlight:false},z:{highlight:false}}},
  {type:"scatter3d", mode:"lines", x:eqXYx, y:eqXYy, z:eqXYz, line:{width:2,color:"#888"}, hoverinfo:"skip", showlegend:false},
  {type:"scatter3d", mode:"lines", x:eqXZx, y:eqXZy, z:eqXZz, line:{width:2,color:"#888"}, hoverinfo:"skip", showlegend:false},
  {type:"scatter3d", mode:"lines", x:eqYZx, y:eqYZy, z:eqYZz, line:{width:2,color:"#888"}, hoverinfo:"skip", showlegend:false},
  {type:"scatter3d", mode:"markers", x:cloudX, y:cloudY, z:cloudZ,
   marker:{size:3, color:cloudC, colorscale:"Viridis", showscale:true, colorbar:{title:"T (K)"}, opacity:0.85},
   name:"samples"},
  {type:"scatter3d", mode:"lines", x:vecX, y:vecY, z:vecZ, line:{width:6,color:"#c00"}, name:"mean state"},
  {type:"scatter3d", mode:"markers", x:[vecX[1]], y:[vecY[1]], z:[vecZ[1]],
   marker:{size:6, color:"#c00"}, showlegend:false},
  {type:"scatter3d", mode:"text", x:labX, y:labY, z:labZ, text:labT, textfont:{size:14}, hoverinfo:"skip", showlegend:false}
];
const layout = {
  title: "$title",
  scene: {aspectmode:"cube",
          xaxis:{title:"S1", range:[-1.3,1.3]},
          yaxis:{title:"S2", range:[-1.3,1.3]},
          zaxis:{title:"S3", range:[-1.3,1.3]}},
  margin:{l:0,r:0,t:40,b:0}
};
Plotly.newPlot("plot", traces, layout, {responsive:true});
</script></body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end

# =====================================================================
# Index
# =====================================================================

const DEMO_INDEX = [
    (group = "path-geometry", fn = demo_path_geometry, kwargs = NamedTuple(),
     desc = "Mixed-segment path: straight leads, circular bends, a catenary sag, and a " *
            "material twist overlay — illustrates the segment assembly API and " *
            "the Frenet–Serret sliding frame."),
    (group = "path-geometry", fn = demo_path_geometry_segment_labels, kwargs = NamedTuple(),
     desc = "Same style of path as the mixed-segment demo, but each segment has a `nickname` " *
            "string; the HTML plot shows those names as 3D labels offset in the osculating plane."),
    (group = "path-geometry", fn = demo_path_geometry_helix_0, kwargs = NamedTuple(),
     desc = "HelixSegment with axis_angle = 0 — winding plane aligned to fiber entry direction."),
    (group = "path-geometry", fn = demo_path_geometry_helix_pi_3, kwargs = NamedTuple(),
     desc = "HelixSegment with axis_angle = π/3 — winding plane rotated by 60°."),
    (group = "path-geometry", fn = demo_path_geometry_helix_2pi_3, kwargs = NamedTuple(),
     desc = "HelixSegment with axis_angle = 2π/3 — winding plane rotated by 120°."),
    (group = "path-geometry", fn = demo_path_geometry_jumps_min_radius, kwargs = NamedTuple(),
     desc = "Demonstrate `jumpby!` and `jumpto!` with focus on the min_bend_radius parameter."),
    (group = "modify", fn = demo_modify_straight_length, kwargs = NamedTuple(),
     desc = "MCMadd(:length) on the first straight of a 3-segment inverted-U baseline."),
    (group = "modify", fn = demo_modify_bend_radius, kwargs = NamedTuple(),
     desc = "MCMadd(:radius) on the bend of a 3-segment inverted-U baseline."),
    (group = "modify", fn = demo_modify_bend_angle, kwargs = NamedTuple(),
     desc = "MCMadd(:angle) on the bend of a 3-segment inverted-U baseline."),
    (group = "modify", fn = demo_modify_straight_length_mul, kwargs = NamedTuple(),
     desc = "MCMmul(:length) on the first straight of a 3-segment inverted-U baseline."),
    (group = "modify", fn = demo_modify_bend_radius_mul, kwargs = NamedTuple(),
     desc = "MCMmul(:radius) on the bend of a 3-segment inverted-U baseline."),
    (group = "modify", fn = demo_modify_bend_angle_mul, kwargs = NamedTuple(),
     desc = "MCMmul(:angle) on the bend of a 3-segment inverted-U baseline."),
    (group = "modify", fn = demo_modify_helix_radius, kwargs = NamedTuple(),
     desc = "MCMadd(:radius) on the helix of a 4-segment (straight · bend · helix · straight) baseline."),
    (group = "modify", fn = demo_modify_helix_pitch, kwargs = NamedTuple(),
     desc = "MCMadd(:pitch) on the helix of a 4-segment baseline."),
    (group = "modify", fn = demo_modify_helix_turns, kwargs = NamedTuple(),
     desc = "MCMadd(:turns) on the helix of a 4-segment baseline."),
    (group = "modify", fn = demo_modify_helix_radius_mul, kwargs = NamedTuple(),
     desc = "MCMmul(:radius) on the helix of a 4-segment baseline."),
    (group = "modify", fn = demo_modify_helix_pitch_mul, kwargs = NamedTuple(),
     desc = "MCMmul(:pitch) on the helix of a 4-segment baseline."),
    (group = "modify", fn = demo_modify_helix_turns_mul, kwargs = NamedTuple(),
     desc = "MCMmul(:turns) on the helix of a 4-segment baseline."),
    (group = "adaptive-step", fn = demo_adaptive_step_doubling, kwargs = NamedTuple(),
     desc = "Adaptive step-doubling diagnostic on a smooth noncommuting generator " *
            "K(s) = α·i·σx·cos(πs) + β·i·σz·sin(2πs). Top panel: accepted/rejected " *
            "step sizes with ‖K(s)‖ overlay; bottom panel: err/tol ratio vs threshold."),
]

"""
    demo_all(; index_output)

Run every demo in `DEMO_INDEX` and write an `demo1.html` that links to each output file
with a short description of what it illustrates.
"""
const _DEMO1_GROUP_TITLES = Dict(
    "path-geometry" => "Path geometry",
    "modify"        => "Modify (MCM parameter perturbations)",
    "adaptive-step" => "Adaptive step-doubling",
)

function demo_all(; index_output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "demo1.html"))
    # entries: (group, title, path, desc)
    entries = Tuple{String, String, String, String}[]

    for d in DEMO_INDEX
        println("[ demo ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        paths = result isa NamedTuple ? values(result) : (result,)
        for v in paths
            if v isa AbstractString && endswith(v, ".html")
                push!(entries, (d.group, basename(v), v, d.desc))
            elseif v isa AbstractVector
                for item in v
                    item isa AbstractString && endswith(item, ".html") &&
                        push!(entries, (d.group, basename(item), item, d.desc))
                end
            end
        end
    end

    # Collect groups in insertion order
    seen_groups = String[]
    for (g, _, _, _) in entries
        g in seen_groups || push!(seen_groups, g)
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST path-geometry demos</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; background: #111; color: #ddd; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #444; padding-bottom: 0.3em; }
    h2   { font-size: 1.15em; margin-top: 1.8em; color: #4db87a; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #4db87a; }
    p.desc { margin: 0.3em 0 0 0; color: #999; font-size: 0.95em; }
    nav.index-nav { font-size: 0.85em; margin-bottom: 1em; color: #666; }
    nav.index-nav a { font-weight: normal; color: #4db87a; margin-right: 0.8em; }
  </style>
</head>
<body>
  <nav class="index-nav">
    <a href="demo1.html">demo1</a>
    <a href="demo2.html">demo2</a>
    <a href="demo3mcm.html">demo3mcm</a>
    <a href="demo3benchmark.html">demo3benchmark</a>
  </nav>
  <h1>BIFROST path-geometry demos</h1>""")
        for g in seen_groups
            heading = get(_DEMO1_GROUP_TITLES, g, g)
            println(io, "  <h2>$(heading)</h2>")
            println(io, "  <ul>")
            for (eg, title, path, desc) in entries
                eg == g || continue
                println(io, "    <li>")
                println(io, "      <a href=\"$(path)\">$(title)</a>")
                println(io, "      <p class=\"desc\">$(desc)</p>")
                println(io, "    </li>")
            end
            println(io, "  </ul>")
        end
        println(io, """</body>
</html>""")
    end

    println("Wrote demo index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_all()
end
