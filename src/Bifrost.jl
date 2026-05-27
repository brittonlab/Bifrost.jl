module Bifrost

# Helper: export every public (non-underscore) binding defined directly in
# the calling submodule. Skips imported names from other modules and the
# submodule's own name. Call at the bottom of each submodule.
function _export_public!(mod::Module)
    self_name = nameof(mod)
    for n in names(mod; all = true)
        s = String(n)
        startswith(s, "#") && continue          # gensyms, anonymous funcs
        startswith(s, "_") && continue          # private convention
        n === self_name && continue
        n === :eval && continue
        n === :include && continue
        # Only export bindings defined in `mod` itself, not re-imports.
        isdefined(mod, n) || continue
        try
            owner = parentmodule(getfield(mod, n))
            owner === mod || continue
        catch
            # parentmodule may fail for some bindings (e.g. constants); export anyway.
        end
        Core.eval(mod, :(export $n))
    end
end

module MaterialProperties
    using LinearAlgebra
    using Printf
    include("material-properties.jl")
    
    include("material/silica.jl")
    include("material/germania.jl")
    include("material/silica-germania.jl")
    include("material/silica-fluorinated.jl")
    
    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

module PathGeometry
    using LinearAlgebra
    # Extend Base.position rather than shadow it, so callers that do
    # `using Bifrost.PathGeometry` keep `position(::PathSpecCached, ...)`
    # working alongside other Base methods.
    import Base: position
    # path-geometry.jl already includes path-geometry-connector.jl internally.
    include("geometry/path-geometry.jl")
    # fiber-path-meta.jl only defines concrete AbstractMeta subtypes
    # (Nickname, MCMadd, MCMmul) plus segment_nickname — it's path-level
    # metadata, not fiber-specific, despite the legacy filename.
    include("fiber/fiber-path-meta.jl")
    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

module FiberCS
    using LinearAlgebra
    using ..MaterialProperties
    include("fiber/fiber-cross-section.jl")

    include("fiber-cross-section/step-index.jl")
    include("fiber-cross-section/graded-index.jl")

    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

module FiberPath
    using LinearAlgebra
    using ..MaterialProperties
    using ..PathGeometry
    # Internal cross-module references used by modify(): these live in
    # PathGeometry and are not exported (underscore-prefixed).
    using ..PathGeometry: _resolve_at_placement, _resolve_twists,
                          _build_quintic_connector, _safe_normalize,
                          _qc_nominalize
    using ..FiberCS
    include("fiber/fiber-path.jl")
    include("fiber/fiber-path-modify.jl")
    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

module PathIntegral
    using LinearAlgebra
    using Printf
    using ..PathGeometry
    using ..FiberCS
    using ..FiberPath
    include("path-integral.jl")
    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

# Umbrella re-export: each core submodule's exported names are surfaced at
# the top level so `using Bifrost` is enough for typical use.
using .MaterialProperties
using .FiberCS
using .PathGeometry
using .FiberPath
using .PathIntegral

for m in (MaterialProperties, FiberCS, PathGeometry, FiberPath, PathIntegral)
    for n in names(m)
        n === nameof(m) && continue
        @eval export $n
    end
end

# Also export the submodule names themselves so callers can write
# `PG = PathGeometry` (or qualified `PathGeometry.X`) after `using Bifrost`.
export MaterialProperties, FiberCS, PathGeometry, FiberPath, PathIntegral

# Plotting lives in a separate `Plots` submodule, opt-in via
# `using Bifrost.Plots`, so plain `using Bifrost` does not pull plotting
# symbols into scope.
module Plots
    using LinearAlgebra
    using ..PathGeometry
    using ..FiberCS
    using ..FiberPath
    using ..PathIntegral
    # Internal helper used by the adaptive-step diagnostic plot.
    using ..PathIntegral: _frobenius_norm
    include("geometry/path-geometry-plot.jl")
    include("fiber/fiber-path-plot.jl")
    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

end # module Bifrost
