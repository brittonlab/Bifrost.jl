using Documenter
using Bifrost

# Hide private (underscore-prefixed) symbols from the rendered API. The `Filter`
# in each `@autodocs` block receives the documented object; `nameof` recovers a
# function/type/module name, while values that lack one (e.g. constants) are
# public by convention and kept.
function is_public_api(obj)
    try
        return !startswith(string(nameof(obj)), "_")
    catch
        return true
    end
end

makedocs(
    sitename = "Bifrost.jl",
    authors  = "Britton",
    modules  = [
        Bifrost,
        Bifrost.MaterialProperties,
        Bifrost.PathGeometry,
        Bifrost.FiberCS,
        Bifrost.FiberPath,
        Bifrost.PathIntegral,
        Bifrost.Plots,
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = ["assets/expand-docstrings.js"],
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Usage" => "usage.md",
        "Theory" => [
            "Deriving the generators" => "generators.md",
            "Path geometry" => "path-geometry.md",
            "Monte Carlo Measurements" => "mcm.md",
        ],
        "API" => [
            "MaterialProperties" => "api/material-properties.md",
            "PathGeometry" => "api/path-geometry.md",
            "FiberCS" => "api/fiber-cross-section.md",
            "FiberPath" => "api/fiber-path.md",
            "PathIntegral" => "api/path-integral.md",
            "Plots" => "api/plots.md",
        ],
    ],
    # Report exported-but-undocumented symbols without failing the initial build.
    checkdocs = :exports,
    warnonly  = true,
    doctest   = false,
)

deploydocs(
    repo = "github.com/brittonlab/BIFROST.git",
    devbranch = "documentor-support",  # TEMP: demo deploy from this branch; restore "main" before merge
    push_preview = true,
)
