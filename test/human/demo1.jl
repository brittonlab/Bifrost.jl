using LinearAlgebra
using MonteCarloMeasurements

using Bifrost
using Bifrost.Plots
using Bifrost.Plots.PlotRuntime
using Bifrost.PathGeometry: _qc_nominalize
include("demo-index-helpers.jl")


const DEMO_FIBER_CROSS_SECTION = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like demo"
)
const DEMO_λ_M = 1550e-9
const DEMO_T_K = 297.15

# =====================================================================
# Helpers
# =====================================================================

# Sample a single placed segment's centerline in the global frame.
# `seg_index` indexes into the list returned by `_all_placed(path)` —
# interior segments first, then the terminal connector.
function _all_placed(path::Union{SubpathBuilt, PathBuilt})
    if path isa SubpathBuilt
        return vcat(path.placed_segments,
                    PlacedSegment[path.jumpto_placed])
    else
        # PathBuilt: flatten across subpaths.
        result = PlacedSegment[]
        offs = s_offsets(path)
        for (i, sp) in enumerate(path.subpaths)
            for ps in sp.placed_segments
                push!(result, PlacedSegment(ps.segment,
                    offs[i] + ps.s_offset_eff, ps.origin, ps.frame,
                    ps.bishop_e1))
            end
            push!(result, PlacedSegment(sp.jumpto_placed.segment,
                offs[i] + sp.jumpto_placed.s_offset_eff,
                sp.jumpto_placed.origin, sp.jumpto_placed.frame,
                sp.jumpto_placed.bishop_e1))
        end
        return result
    end
end

function _sample_segment_xyz(path::Union{SubpathBuilt, PathBuilt},
                             seg_index::Int; n::Int = 128)
    placed = _all_placed(path)
    ps = placed[seg_index]
    s0 = Float64(_qc_nominalize(ps.s_offset_eff))
    s1 = s0 + Float64(_qc_nominalize(arc_length(ps.segment)))
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

# =====================================================================
# Modify demos — helpers shared by all 12 rows
# =====================================================================

# Build the baseline 3-segment inverted-U and attach `target_meta` to the
# segment indexed by `target_idx` (1 = first straight, 2 = bend, 3 = second
# straight). Returns the modified path (a SubpathBuilt).
function _build_modify_variant(L::Float64, R::Float64,
                               target_idx::Int, target_meta::AbstractVector{<:AbstractMeta})
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = target_idx == 1 ? target_meta : AbstractMeta[])
    bend!(sb; radius = R, angle = π, axis_angle = 0.0,
          meta = target_idx == 2 ? target_meta : AbstractMeta[])
    straight!(sb; length = L, meta = target_idx == 3 ? target_meta : AbstractMeta[])
    seal!(sb)
    # Fiber(builder) applies the segment meta (field-level MCM here) during the
    # single build; no separate modify step.
    return Fiber(sb; cross_section = DEMO_FIBER_CROSS_SECTION,
                 T_ref_K = DEMO_T_K).path
end

# Baseline helix parameters for the 4-segment variant (inverted-U plus a helix
# between the bend and the second straight).
const _MODIFY_HELIX_RADIUS = 0.15
const _MODIFY_HELIX_PITCH  = 0.25
const _MODIFY_HELIX_TURNS  = 1.5

# Build a 4-segment baseline (straight · bend π · helix · straight) and attach
# `target_meta` to the segment indexed by `target_idx`:
#   1 = first straight, 2 = bend, 3 = helix, 4 = second straight.
function _build_modify_variant_helix(L::Float64, R::Float64,
                                     target_idx::Int, target_meta::AbstractVector{<:AbstractMeta})
    sb = SubpathBuilder(); start!(sb)
    straight!(sb; length = L, meta = target_idx == 1 ? target_meta : AbstractMeta[])
    bend!(sb; radius = R, angle = π, axis_angle = 0.0,
          meta = target_idx == 2 ? target_meta : AbstractMeta[])
    helix!(sb; radius = _MODIFY_HELIX_RADIUS,
                 pitch  = _MODIFY_HELIX_PITCH,
                 turns  = _MODIFY_HELIX_TURNS,
                 axis_angle = 0.0,
                 meta = target_idx == 3 ? target_meta : AbstractMeta[])
    straight!(sb; length = L, meta = target_idx == 4 ? target_meta : AbstractMeta[])
    seal!(sb)
    # Fiber(builder) applies the segment meta (field-level MCM here) during the
    # single build; no separate modify step.
    return Fiber(sb; cross_section = DEMO_FIBER_CROSS_SECTION,
                 T_ref_K = DEMO_T_K).path
end

# Render one experiment row: each variant's path laid out side-by-side along x,
# with the target segment drawn in red and the others in green. The terminal
# connector is rendered in faint gray.
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
        all_segs = _all_placed(path)
        n_interior = length(path.placed_segments)
        for i in eachindex(all_segs)
            s     = _sample_segment_xyz(path, i)
            xs    = s.x .+ dx
            color = if i > n_interior
                "#444"   # terminal connector — faint
            elseif i == target_idx && !isempty(target_meta)
                "#e55"   # red highlight on perturbed segment
            else
                "#6c6"   # baseline green
            end
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
        p0 = position(path, 0.0)
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
# Helix with material spin (ported from the old demo.jl)
# =====================================================================

function demo_helix_mcm_spin(;
    output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "helix-mcm-spin.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Helix with MCM spin",
)
    # Geometry-only demo. Uses the PathGeometry module wrapper from
    # path-geometry-plot.jl so the dispatch on `PathGeometry.SubpathBuilt`
    # in `write_path_geometry_plot3d` resolves correctly.
    PG = PathGeometry
    sb = PG.SubpathBuilder(); PG.start!(sb; spin_rate = 2π)
    PG.straight!(sb; length = 1.0, meta = [PG.Nickname("lead-in")])
    PG.helix!(sb; radius = 0.5, pitch = 0.05, turns = 4.0, axis_angle = 0.0,
              meta = [PG.Nickname("helix")])
    PG.straight!(sb; length = 1.0, meta = [PG.Nickname("lead-out")])
    # Seal at the natural exit of the trailing straight — no terminal connector
    # bending. After the axis_angle=0 helix on +z the exit tangent has x and z
    # components; the trailing straight advances along it and seal! ends there.
    PG.seal!(sb)
    b = PG.build(sb)
    return write_path_geometry_plot3d(
        b, 0.0, Float64(PG._qc_nominalize(PG.s_end(b)));
        fidelity = fidelity,
        output = output,
        title = title,
    )
end

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
# Index
# =====================================================================

const DEMO_INDEX = [
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
    (group = "spin", fn = demo_helix_mcm_spin, kwargs = NamedTuple(),
     desc = "Helix with a constant material spin rate set via start!(; spin_rate). " *
            "Demonstrates that material spin propagates through the geometry " *
            "layer's Frenet frame independent of segment torsion."),
    (group = "adaptive-step", fn = demo_adaptive_step_doubling, kwargs = NamedTuple(),
     desc = "Adaptive step-doubling diagnostic on a smooth noncommuting generator " *
            "K(s) = α·i·σx·cos(πs) + β·i·σz·sin(2πs). Top panel: accepted/rejected " *
            "step sizes with ‖K(s)‖ overlay; bottom panel: err/tol ratio vs threshold."),
]

"""
    demo_all(; index_output)

Run every demo in `DEMO_INDEX` and write `demo-index.html` linking to each
output file with a short description.

Geometry-only demos live in `demo-path-geometry.jl`; this file covers
modify-pipeline demos plus the adaptive step-doubling and helix-mcm-spin
diagnostics.
"""
const _DEMO1_GROUP_TITLES = Dict(
    "modify"        => "Modify (MCM parameter perturbations)",
    "spin"         => "Spin",
    "adaptive-step" => "Adaptive step-doubling",
)

function demo_entries()
    entries = Tuple{String, String, String, String}[]

    for d in DEMO_INDEX
        println("[ demo ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        desc = _demo_result_desc(result, d)
        for path in _demo_html_paths(result)
            push!(entries, (d.group, basename(path), path, desc))
        end
    end
    return entries
end

function demo_all(; index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT)
    return _write_demo_index(
        [(title = "Modify + diagnostics demos",
          entries = demo_entries(),
          group_titles = _DEMO1_GROUP_TITLES)];
        index_output,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_all()
end
