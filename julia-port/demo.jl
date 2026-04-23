include("material-properties.jl")
include("fiber-cross-section.jl")
include("path-geometry.jl")
include("path-geometry-plot.jl")



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

# """
#     demofiber1()
# The shape of the fiber in this file is that of polarization paddles.
# The first loop is L/4, the second is L/2, and the third is L/4, where 
# L is the length of fiber in each loop. The angles of the loops are 0, π/2, 
# and 0, respectively. The input state is vertical polarization.
# """
# function demofiber1()
#     lead_in = 4.0
#     spacer = 8.5
#     lead_out = 4.0
#     loop_turns = [1.0, 1.0, 1.0]
#     target_retardance = [π / 2, π, π / 2]
#     loop_radii = 2π .* loop_turns ./ target_retardance
#     loop_angles = [0.0, π / 2, 0.0]

#     total_length = lead_in + sum(2π .* loop_radii .* loop_turns) + 2 * spacer + lead_out
#     spec = FiberSpec(0.0, total_length; cross_section = DEMO_FIBER_CROSS_SECTION, λ_m = DEMO_λ_M)

#     s = 0.0
#     bend!(spec, s, s + lead_in; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += lead_in

#     loop1_angle = 2π * loop_turns[1]
#     loop1_length = loop_radii[1] * loop1_angle
#     bend!(spec, s, s + loop1_length; angle = loop1_angle, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += loop1_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[2] - loop_angles[1]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop2_angle = 2π * loop_turns[2]
#     loop2_length = loop_radii[2] * loop2_angle
#     bend!(spec, s, s + loop2_length; angle = loop2_angle, axis = loop_angles[2], T_K = DEMO_T_K)
#     s += loop2_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[2], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[3] - loop_angles[2]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop3_angle = 2π * loop_turns[3]
#     loop3_length = loop_radii[3] * loop3_angle
#     bend!(spec, s, s + loop3_length; angle = loop3_angle, axis = loop_angles[3], T_K = DEMO_T_K)
#     s += loop3_length

#     bend!(spec, s, s + lead_out; angle = 0.0, axis = loop_angles[3], T_K = DEMO_T_K)

#     return (
#         fiber = build(spec),
#         title = "Fiber polarization control paddles",
#         input_state = ComplexF64[0.0 + 0.0im, 1.0 + 0.0im],
#         output = "fiber-path-plot.html",
#         n = 4001
#     )
# end


# """
#     demofiber2()
# This fiber is a more complex version of the first demo, with four loops instead of three, 
# and with the loop angles changing by π/3 between each loop. The retardances of the loops 
# are π/2, π, π/2, and π/3. The input state is diagonal polarization.  
# """
# function demofiber2()
#     lead_in = 3.0
#     spacer = 6.0
#     lead_out = 3.0
#     loop_turns = [0.75, 1.0, 0.75, 0.5]
#     target_retardance = [π / 2, π, π / 2, π / 3]
#     loop_radii = 2π .* loop_turns ./ target_retardance
#     loop_angles = [0.0, π / 3, 2π / 3, π / 6]

#     total_length = lead_in + sum(2π .* loop_radii .* loop_turns) + 3 * spacer + lead_out
#     spec = FiberSpec(0.0, total_length; cross_section = DEMO_FIBER_CROSS_SECTION, λ_m = DEMO_λ_M)

#     s = 0.0
#     bend!(spec, s, s + lead_in; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += lead_in

#     loop1_angle = 2π * loop_turns[1]
#     loop1_length = loop_radii[1] * loop1_angle
#     bend!(spec, s, s + loop1_length; angle = loop1_angle, axis = loop_angles[1], T_K = DEMO_T_K)
#     s += loop1_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[1], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[2] - loop_angles[1]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop2_angle = 2π * loop_turns[2]
#     loop2_length = loop_radii[2] * loop2_angle
#     bend!(spec, s, s + loop2_length; angle = loop2_angle, axis = loop_angles[2], T_K = DEMO_T_K)
#     s += loop2_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[2], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[3] - loop_angles[2]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop3_angle = 2π * loop_turns[3]
#     loop3_length = loop_radii[3] * loop3_angle
#     bend!(spec, s, s + loop3_length; angle = loop3_angle, axis = loop_angles[3], T_K = DEMO_T_K)
#     s += loop3_length

#     bend!(spec, s, s + spacer; angle = 0.0, axis = loop_angles[3], T_K = DEMO_T_K)
#     twist!(spec, s, s + spacer; rate = (loop_angles[4] - loop_angles[3]) / spacer, T_K = DEMO_T_K)
#     s += spacer

#     loop4_angle = 2π * loop_turns[4]
#     loop4_length = loop_radii[4] * loop4_angle
#     bend!(spec, s, s + loop4_length; angle = loop4_angle, axis = loop_angles[4], T_K = DEMO_T_K)
#     s += loop4_length

#     bend!(spec, s, s + lead_out; angle = 0.0, axis = loop_angles[4], T_K = DEMO_T_K)

#     return (
#         fiber = build(spec),
#         title = "Four-loop twisted fiber controller",
#         input_state = ComplexF64[inv(sqrt(2)) + 0.0im, 0.0 + inv(sqrt(2)) * im],
#         output = "fiber-path-plot.html",
#         n = 4001
#     )
# end

# function run_demo(demo ; output::Union{Nothing,AbstractString} = nothing)
#     fiber = demo.fiber
#     title = demo.title
#     input_state = demo.input_state
#     plot_output = isnothing(output) ? demo.output : output
#     n = demo.n

#     J_final, stats = propagate_fiber(
#         fiber;
#         rtol = 1e-10,
#         atol = 1e-12,
#         h_init = 1e-2
#     )

#     println("Final Jones matrix:")
#     println(J_final)

#     println("\nInterval stats:")
#     for (i, st) in enumerate(stats)
#         println("interval $i: accepted=$(st.accepted_steps), rejected=$(st.rejected_steps)")
#     end

#     plot_path = write_fiber_input_plot3d(
#         fiber,
#         fiber.s_start,
#         fiber.s_end;
#         n = n,
#         output = plot_output,
#         title = title,
#         input_state = input_state
#     )

#     println("\nWrote interactive 3D plot to:")
#     println(plot_path)

#     return (; demo..., J_final, stats, plot_path)
# end

"""
    demo_fiber_path(; output, n, title)

Demonstrate the analytic 3D path layer in [`path-geometry.jl`](path-geometry.jl) (via the
`PathGeometry` module in [`path-geometry-plot.jl`](path-geometry-plot.jl)): a short centerline
with circular bends, a catenary segment, and a material twist overlay on a straight spacer.
Writes an interactive Plotly HTML file.
"""
function demo_fiber_path(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Fiber path geometry: bends, catenary, twist",
)
    PG = PathGeometry
    spec = PG.PathSpec()
    L_lead = 0.10
    R1 = 0.05
    θ1 = π / 2
    L_spacer = 0.12
    PG.straight!(spec; length = L_lead)
    PG.bend!(spec; radius = R1, angle = θ1)
    PG.straight!(spec; length = L_spacer)
    PG.catenary!(spec; a = 0.03, length = 0.10, axis_angle = 0.0)
    PG.bend!(spec; radius = 0.06, angle = π / 3)
    PG.straight!(spec; length = 0.08)
    s_twist = L_lead + R1 * abs(θ1)
    PG.twist!(spec; s_start = s_twist, length = L_spacer, rate = 35.0)
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path))
    println("Writhe: ", PG.writhe(path))
    plot_path = write_path_geometry_plot3d(
        path,
        path.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    println("Wrote path geometry plot to:")
    println(plot_path)
    return (; path, plot_path)
end

"""
    demo_fiber_path_segment_labels(; output, fidelity, title)

Illustrate optional **segment nicknames**: `straight!`, `bend!`, `catenary!`, and `helix!` are
called with `nickname` strings. [`write_path_geometry_plot3d`](path-geometry-plot.jl) renders each
name as 3D text near the segment midpoint (offset along the principal normal, in the osculating
plane). Uses only straight / bend / catenary / helix segments so the demo does not depend on jump
connectors or `min_bend_radius`.
"""
function demo_fiber_path_segment_labels(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-segment-labels-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Path geometry: segment nicknames",
)
    PG = PathGeometry
    spec = PG.PathSpec()
    PG.straight!(spec; length = 0.08, nickname = "lead-in")
    PG.bend!(spec; radius = 0.06, angle = π / 2, nickname = "90° bend")
    PG.straight!(spec; length = 0.06, nickname = "spacer")
    PG.catenary!(spec; a = 0.04, length = 0.08, axis_angle = 0.0, nickname = "sag")
    PG.helix!(spec; radius = 0.025, pitch = 0.015, turns = 1.2, axis_angle = 0.0, nickname = "twist section")
    PG.straight!(spec; length = 0.06, nickname = "lead-out")
    path = PG.build(spec)
    println("Arc length (effective): ", PG.path_length(path), " m")
    plot_path = write_path_geometry_plot3d(
        path,
        path.s_start,
        path.s_end;
        fidelity = fidelity,
        output = output,
        title = title,
    )
    println("Wrote segment-label demo to:")
    println(plot_path)
    return (; path, plot_path)
end

"""
    demo_fiber_path_helix(; output, n, title)

Demonstrate `HelixSegment` with three separate single-segment paths, each a helix
with a different `axis_angle` (0, π/3, 2π/3).  Each helix has the same radius and
pitch, so the shape is the same but the winding plane rotates around the fiber axis.
Writes an interactive Plotly HTML file via `write_path_geometry_plot3d`.
"""
function demo_fiber_path_helix(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-helix-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "HelixSegment: three axis angles",
)
    PG = PathGeometry
    paths = map([0.0, π/3, 2π/3]) do axis_angle
        spec = PG.PathSpec()
        PG.straight!(spec; length = 0.05)
        PG.helix!(spec; radius = 0.03, pitch = 0.02, turns = 2.0, axis_angle)
        PG.straight!(spec; length = 0.05)
        PG.build(spec)
    end
    for (i, path) in enumerate(paths)
        println("Helix $(i) (axis_angle=$(round((i-1)*π/3, digits=3)) rad): " *
                "arc_length=$(round(PG.path_length(path), digits=4)) m")
    end
    base, ext = splitext(output)
    plot_paths = map(enumerate(paths)) do (i, path)
        angle_label = ["0", "pi_3", "2pi_3"][i]
        out = "$(base)_$(angle_label)$(ext)"
        write_path_geometry_plot3d(
            path, path.s_start, path.s_end;
            fidelity = fidelity,
            output = out,
            title = "$(title) — axis_angle=$(["0", "π/3", "2π/3"][i])",
        )
    end
    for p in plot_paths
        println("Wrote helix demo plot to: ", p)
    end
    return (; paths, plot_paths)
end

"""
    demo_fiber_path_jumps_min_radius(; output, fidelity, title)

Demonstrate `jumpby!` and `jumpto!` with focus on the min_bend_radius parameter.
"""
function demo_fiber_path_jumps_min_radius(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-jumps-demo.html"),
    fidelity::Float64 = 4.0,
    title::AbstractString = "JumpBy and JumpTo: Hermite connectors, min_bend_radius",
)
    PG = PathGeometry
    spec = PG.PathSpec()

    PG.straight!(spec; length = 1)
    PG.jumpto!(spec; destination = (1, 0.0, 1), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.4) # T1: fails only if min_bend_radius >0.5
    PG.straight!(spec; length = 1)
    PG.jumpto!(spec; destination = (2, 0.0, 0), tangent = (0.0, 0.0, 1.0),
            min_bend_radius = 0.1) # T2: fails only if min_bend_radius >0.5
    PG.straight!(spec; length = 1)
    PG.jumpto!(spec; destination = (3, 0.0, 1), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.2) # T3: fails only if min_bend_radius >0.50
    PG.straight!(spec; length = 1)
    PG.jumpby!(spec; delta = (-1, 0.0, 0), tangent = (0.0, 0.0, -1.0),
            min_bend_radius = 0.1)
    PG.straight!(spec; length = 1)

    path = PG.build(spec)

    plot_path = write_path_geometry_plot3d(
        path, path.s_start, path.s_end;
        fidelity = fidelity, output   = output, title    = title )
    return (; path, plot_path)
end

"""
    demo_fiber_path_jumps2(; output, fidelity, title)

Demonstrate a more complex routing scenario with multiple `jumpby!` and `jumpto!`
calls interspersed with straights, a helix, and a circular bend.

Path layout:
  1. Short straight lead-in.
  2. **JumpBy** with an explicit `tangent_out` — the Hermite connector prescribes
     the exit direction (pointing diagonally upward), not the chord direction.
  3. Single-turn helix section.
  4. Short straight spacer.
  5. **JumpBy** without `tangent_out` — exit tangent follows the chord.
  6. Quarter-circle bend.
  7. **JumpTo** to a fixed lab-frame waypoint with no prescribed tangent.
  8. Short straight.
  9. **JumpTo** to a second waypoint with a prescribed outgoing tangent (global +z).
  10. Straight lead-out.
"""
function demo_fiber_path_jumps2(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "path-geometry-jumps2-demo.html"),
    fidelity::Float64 = 1.0,
    title::AbstractString = "Multi-jump routing: JumpBy (with/without tangent_out) + JumpTo",
)
    PG = PathGeometry
    spec = PG.PathSpec()

    PG.straight!(spec; length = 0.05)

    # Prescribed tangent_out + min_bend_radius: handle length grows to keep κ ≤ 1/R.
    PG.jumpby!(spec; delta = (0.03, 0.0, 0.04), tangent = (0.0, 0.0, 1.0),
               min_bend_radius = 0.02)

    PG.helix!(spec; radius = 0.02, pitch = 0.01, turns = 1.0, axis_angle = 0.0)
    PG.straight!(spec; length = 0.04)

    # No tangent_out → exit tangent follows chord direction.
    PG.jumpby!(spec; delta = (0.05, 0.0, 0.03))

    PG.bend!(spec; radius = 0.04, angle = π / 2, axis_angle = 0.0)
    PG.jumpto!(spec; destination = (0.22, 0.0, 0.20))
    PG.straight!(spec; length = 0.03)
    PG.jumpto!(spec; destination = (0.22, 0.0, 0.28), tangent = (0.0, 0.0, 1.0),
               min_bend_radius = 0.1)
    PG.straight!(spec; length = 0.05)

    path = PG.build(spec)
    println("Arc length : ", round(PG.path_length(path) * 100; digits = 2), " cm")
    println("Start      : ", round.(PG.start_point(path) .* 100; digits = 1), " cm")
    println("End        : ", round.(PG.end_point(path) .* 100; digits = 1), " cm")
    println("Writhe     : ", round(PG.writhe(path); digits = 6))

    plot_path = write_path_geometry_plot3d(
        path, path.s_start, path.s_end;
        fidelity = fidelity,
        output   = output,
        title    = title,
    )
    println("Wrote jumps2 demo to: ", plot_path)
    return (; path, plot_path)
end

const DEMO_INDEX = [
    (
        fn   = demo_fiber_path,
        kwargs = NamedTuple(),
        desc = "Mixed-segment path: straight leads, circular bends, a catenary sag, and a " *
               "material twist overlay — illustrates the segment assembly API and " *
               "the Frenet–Serret sliding frame.",
    ),
    (
        fn   = demo_fiber_path_segment_labels,
        kwargs = NamedTuple(),
        desc = "Same style of path as the mixed-segment demo, but each segment has a `nickname` " *
               "string; the HTML plot shows those names as 3D labels offset in the osculating plane.",
    ),
    (
        fn   = demo_fiber_path_helix,
        kwargs = NamedTuple(),
        desc = "Three HelixSegment paths sharing the same radius and pitch but with " *
               "axis_angle ∈ {0, π/3, 2π/3} — shows how the winding plane rotates " *
               "while the entry tangent stays aligned with the incoming fiber direction.",
    ),
    (
        fn   = demo_fiber_path_jumps_min_radius,
        kwargs = NamedTuple(),
        desc = "Demonstrate `jumpby!` and `jumpto!` with focus on the min_bend_radius parameter." 
    ),
    (
        fn   = demo_fiber_path_jumps2,
        kwargs = NamedTuple(),
        desc = "Multi-jump routing: two JumpBy calls (one with, one without tangent_out) " *
               "and two JumpTo calls (one with, one without tangent_out) interleaved with " *
               "a helix, a circular bend, and straight spacers — exercises every " *
               "combination of jump type × tangent prescription.",
    ),
]

"""
    demo_all(; index_output)

Run every demo in `DEMO_INDEX` and write an `index.html` that links to each output file
with a short description of what it illustrates.
"""
function demo_all(; index_output::AbstractString = joinpath(@__DIR__, "..", "output", "index.html"))
    entries = Tuple{String, String, String}[]   # (title, rel_path, desc)

    for d in DEMO_INDEX
        result = d.fn(; d.kwargs...)
        # Collect all HTML output paths from the result (handles single or vector).
        paths = result isa NamedTuple ? values(result) : (result,)
        html_paths = String[]
        for v in paths
            if v isa AbstractString && endswith(v, ".html")
                push!(html_paths, v)
            elseif v isa AbstractVector
                for item in v
                    item isa AbstractString && endswith(item, ".html") && push!(html_paths, item)
                end
            end
        end
        for p in html_paths
            push!(entries, (basename(p), p, d.desc))
        end
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST path-geometry demos</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; color: #222; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #ccc; padding-bottom: 0.3em; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #1a6; }
    p.desc { margin: 0.3em 0 0 0; color: #555; font-size: 0.95em; }
  </style>
</head>
<body>
  <h1>BIFROST path-geometry demos</h1>
  <ul>""")
        for (title, path, desc) in entries
            println(io, "    <li>")
            println(io, "      <a href=\"$(path)\">$(title)</a>")
            println(io, "      <p class=\"desc\">$(desc)</p>")
            println(io, "    </li>")
        end
        println(io, """  </ul>
</body>
</html>""")
    end

    println("Wrote demo index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_all()
end
