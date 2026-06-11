# =====================================================================
# demo-index.jl
#
# Main entry point for running every visual demo registry and writing one
# monolithic index:
#
#     output/demo-index.html
#
# Individual demo functions still produce their own HTML artifacts. This file
# replaces the former partial list pages such as demo1.html, demo2.html,
# demo3mcm.html, demo3benchmark.html, and demo-path-geometry-index.html.
# =====================================================================

include("demo-index-helpers.jl")

const DEMO_INDEX_DEMO_FILES = [
    "demo1.jl",
    "demo2.jl",
#    "demo3mcm.jl",
#    "demo3benchmark.jl",
    "demo-path-geometry.jl",
    "demo-frame-birefringence.jl",
]

const DEMO_INDEX_SECTION_REGISTRY = Dict(
    "demo-path-geometry.jl" => (
        title = "Path-geometry-only demos",
        entry_fn_name = :demo_path_geometry_entries,
        group_titles_name = nothing,
        group_titles = Dict{String, String}(),
    ),
    "demo1.jl" => (
        title = "Modify + diagnostics demos",
        entry_fn_name = :demo_entries,
        group_titles_name = :_DEMO1_GROUP_TITLES,
        group_titles = Dict{String, String}(),
    ),
    "demo2.jl" => (
        title = "JumpBy / JumpTo demos",
        entry_fn_name = :demo2_entries,
        group_titles_name = :_DEMO2_GROUP_TITLES,
        group_titles = Dict{String, String}(),
    ),
    "demo3mcm.jl" => (
        title = "MCM temperature PTF demos",
        entry_fn_name = :demo3mcm_entries,
        group_titles_name = nothing,
        group_titles = Dict("MCM" => "MCM / temperature-dependent PTF"),
    ),
    "demo3benchmark.jl" => (
        title = "MCM propagation benchmarks",
        entry_fn_name = :demo3benchmark_entries,
        group_titles_name = nothing,
        group_titles = Dict("benchmarks" => "Speed benchmarks"),
    ),
    "demo-frame-birefringence.jl" => (
        title = "Transported (Bishop) frame and birefringence",
        entry_fn_name = :demo_frame_birefringence_entries,
        group_titles_name = nothing,
        group_titles = Dict(
            "frame" => "Frame gauge: primer, sources, pathologies, anchor"),
    ),
)

function _include_demo_file_once(file::AbstractString)
    include(joinpath(@__DIR__, file))
    return nothing
end

function include_demo_files()
    for file in DEMO_INDEX_DEMO_FILES
        _include_demo_file_once(file)
    end
    return nothing
end

include_demo_files()

function _demo_index_resolve_section(file::AbstractString)
    haskey(DEMO_INDEX_SECTION_REGISTRY, file) ||
        error("No demo-index section registered for $(file)")

    spec = DEMO_INDEX_SECTION_REGISTRY[file]
    isdefined(Main, spec.entry_fn_name) ||
        error("$(file) did not define $(spec.entry_fn_name)")

    group_titles = isnothing(spec.group_titles_name) ? spec.group_titles :
                   getfield(Main, spec.group_titles_name)

    return (
        title = spec.title,
        source_file = file,
        entry_fn = getfield(Main, spec.entry_fn_name),
        group_titles = group_titles,
    )
end

function demo_index_sections()
    section_specs = [_demo_index_resolve_section(file) for file in DEMO_INDEX_DEMO_FILES]
    sections = []
    for section in section_specs
        println("[ demo-index ] running $(section.source_file)")
        push!(sections, (
            title = section.title,
            source_file = section.source_file,
            entries = section.entry_fn(),
            group_titles = section.group_titles,
        ))
    end
    return sections
end

function demo_index_all(;
    index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT,
)
    mkpath(dirname(index_output))
    return _write_demo_index(demo_index_sections(); index_output)
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_index_all()
end
