using Documenter
using Bifrost

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
        assets = ["assets/expand-modules.js"],
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
