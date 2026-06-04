if !isdefined(@__MODULE__, :DEMO_MONOLITHIC_INDEX_OUTPUT)
    const DEMO_MONOLITHIC_INDEX_OUTPUT = joinpath(@__DIR__, "..", "..", "output", "demo-index.html")
end

function _demo_html_paths(result)
    values_to_scan = result isa NamedTuple ? values(result) : (result,)
    paths = String[]
    for v in values_to_scan
        if v isa AbstractString && endswith(v, ".html")
            push!(paths, String(v))
        elseif v isa AbstractVector
            for item in v
                item isa AbstractString && endswith(item, ".html") && push!(paths, String(item))
            end
        end
    end
    return paths
end

function _demo_result_desc(result, entry)
    desc_inline = (result isa NamedTuple && haskey(result, :desc)) ? String(result.desc) : ""
    desc_entry = hasproperty(entry, :desc) ? String(entry.desc) : ""
    return isempty(desc_inline) ? desc_entry : desc_inline
end

function _demo_seen_groups(entries)
    seen_groups = String[]
    for (g, _, _, _) in entries
        g in seen_groups || push!(seen_groups, g)
    end
    return seen_groups
end

function _write_demo_index(
    sections;
    index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT,
    title::AbstractString = "BIFROST demo index",
)
    isdir(dirname(index_output)) || mkpath(dirname(index_output))
    open(index_output, "w") do io
        println(io, """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>$(title)</title>
  <style>
    body { font-family: sans-serif; max-width: 900px; margin: 2em auto; background: #111; color: #ddd; }
    h1   { font-size: 1.6em; border-bottom: 1px solid #444; padding-bottom: 0.3em; }
    h2   { font-size: 1.25em; margin-top: 2.0em; color: #4db87a; }
    h3   { font-size: 1.05em; margin-top: 1.5em; color: #ddd; }
    ul   { padding-left: 1.2em; }
    li   { margin: 1em 0; }
    a    { font-weight: bold; color: #4db87a; }
    p.desc { margin: 0.3em 0 0 0; color: #999; font-size: 0.95em; }
    p.source { margin: -0.6em 0 1.2em 0; color: #888; font-size: 0.9em; }
    code { color: #bbb; }
  </style>
</head>
<body>
  <h1>$(title)</h1>""")
        for section in sections
            entries = section.entries
            isempty(entries) && continue
            println(io, "  <h2>$(section.title)</h2>")
            if haskey(section, :source_file)
                println(io, "  <p class=\"source\">Source: <code>$(section.source_file)</code></p>")
            end
            for group in _demo_seen_groups(entries)
                heading = get(section.group_titles, group, group)
                println(io, "  <h3>$(heading)</h3>")
                println(io, "  <ul>")
                for (eg, link_title, path, desc) in entries
                    eg == group || continue
                    println(io, "    <li>")
                    println(io, "      <a href=\"$(basename(path))\">$(link_title)</a>")
                    isempty(desc) || println(io, "      <p class=\"desc\">$(desc)</p>")
                    println(io, "    </li>")
                end
                println(io, "  </ul>")
            end
        end
        println(io, """</body>
</html>""")
    end
    println("Wrote demo index to: ", index_output)
    return index_output
end
