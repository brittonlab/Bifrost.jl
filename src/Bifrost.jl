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
    # path-geometry.jl includes path-geometry-connector.jl AND
    # path-geometry-meta.jl internally (the Subpath constructors reference
    # MCMadd/MCMmul for validation, so meta must load as part of the geometry
    # layer). We therefore do NOT include path-geometry-meta.jl separately here —
    # that would double-define Nickname/MCMadd/MCMmul.
    include("geometry/path-geometry.jl")
    # Material-agnostic perturbation mechanism used by build(...; perturb=true)
    # and by consuming layers (the fiber) for isotropic scaling. Included after
    # path-geometry.jl so its segment types and meta vocabulary are in scope.
    include("geometry/path-geometry-perturb.jl")
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
    # Internal cross-module references. The geometry-layer perturbation mechanism
    # (used by the fiber's thermal :T_K interpretation) and a couple of helpers
    # are underscore-prefixed and not exported.
    using ..PathGeometry: _scale_length_fields, _meta_without, _length_fields,
                          _qc_nominalize, _resolve_inherited_start
    using ..FiberCS
    include("fiber/fiber-path.jl")
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
