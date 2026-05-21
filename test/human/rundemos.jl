using LinearAlgebra
using MonteCarloMeasurements
using Bifrost
using Bifrost.Plots
using Bifrost.Plots.PlotRuntime

# Include every demo file. demo3benchmark transitively includes demo4mcm,
# and demo2 transitively includes demo1, but listing them explicitly keeps
# the load order deterministic and the file dependencies obvious here.
include(joinpath(@__DIR__, "demo1.jl"))
include(joinpath(@__DIR__, "demo2.jl"))
include(joinpath(@__DIR__, "demo4mcm.jl"))
include(joinpath(@__DIR__, "demo3benchmark.jl"))

# =====================================================================
# Aggregate index — one section per source demo file. Each section's
# entries come from that file's own DEMO*_INDEX constant so the
# descriptions stay next to the implementations.
# =====================================================================

# Source files in display order. Each tuple:
#   (source_basename, section_heading, index_constant, group_titles_or_nothing)
const _RUNDEMOS_SOURCES = [
    ("demo1.jl",           "demo1 — path geometry, modify, adaptive-step", DEMO_INDEX,            _DEMO1_GROUP_TITLES),
    ("demo2.jl",           "demo2 — JumpBy / JumpTo",                      DEMO2_INDEX,           nothing),
    ("demo4mcm.jl",        "demo4mcm — MCM temperature PTF",               DEMO4MCM_INDEX,        nothing),
    ("demo3benchmark.jl",  "demo3benchmark — MCM propagation timing",      DEMO3BENCHMARK_INDEX,  nothing),
]

"""
    demo_all(; index_output)

Run every demo registered in each demo file's `DEMO*_INDEX` constant and
write a single `index.html` linking to every output file.
"""
function demo_all(; index_output::AbstractString = joinpath(@__DIR__, "..", "..", "output", "index.html"))
    # entries_by_source: source_basename => Vector{(group, title, path, desc)}
    entries_by_source = Dict{String, Vector{Tuple{String, String, String, String}}}()

    for (source, _heading, index, _groups) in _RUNDEMOS_SOURCES
        entries = Tuple{String, String, String, String}[]
        for d in index
            println("[ $(source) ] $(d.fn)")
            result = d.fn(; d.kwargs...)
            desc_inline = (result isa NamedTuple && haskey(result, :desc)) ?
                          String(result.desc) : nothing
            desc_entry  = hasproperty(d, :desc) ? d.desc : ""
            desc        = isnothing(desc_inline) ? desc_entry : desc_inline
            group       = hasproperty(d, :group) ? d.group : ""

            paths = result isa NamedTuple ? values(result) : (result,)
            for v in paths
                if v isa AbstractString && endswith(v, ".html")
                    push!(entries, (group, basename(v), v, desc))
                elseif v isa AbstractVector
                    for item in v
                        item isa AbstractString && endswith(item, ".html") &&
                            push!(entries, (group, basename(item), item, desc))
                    end
                end
            end
        end
        entries_by_source[source] = entries
    end

    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>BIFROST demos</title>
  <style>
    body { font-family: sans-serif; max-width: 800px; margin: 2em auto; background: #111; color: #ddd; }
    h1   { font-size: 1.5em; border-bottom: 1px solid #444; padding-bottom: 0.3em; }
    h2   { font-size: 1.2em;  margin-top: 1.8em; color: #4db87a; }
    h3   { font-size: 1.0em;  margin-top: 1.2em; color: #9bb; font-weight: normal; }
    ul   { padding-left: 1.2em; }
    li   { margin: 0.8em 0; }
    a    { font-weight: bold; color: #4db87a; }
    p.desc { margin: 0.3em 0 0 0; color: #999; font-size: 0.95em; }
  </style>
</head>
<body>
  <h1>BIFROST demos</h1>""")
        for (source, heading, _index, group_titles) in _RUNDEMOS_SOURCES
            entries = entries_by_source[source]
            isempty(entries) && continue
            println(io, "  <h2>$(heading)</h2>")
            # Collect this source's groups in insertion order.
            seen_groups = String[]
            for (g, _, _, _) in entries
                g in seen_groups || push!(seen_groups, g)
            end
            for g in seen_groups
                if !isempty(g)
                    heading2 = group_titles isa Dict ? get(group_titles, g, g) : g
                    println(io, "  <h3>$(heading2)</h3>")
                end
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
        end
        println(io, "</body>\n</html>")
    end

    println("Wrote demo index to: ", index_output)
    return index_output
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_all()
end
