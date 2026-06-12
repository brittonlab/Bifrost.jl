"""
Birefringence simulation for optical fiber: build three-dimensional fiber paths, bind
them to cross-section optics, and propagate Jones matrices with adaptive error control.

`using Bifrost` re-exports the public names of the core submodules
(`MaterialProperties`, `PathGeometry`, `FiberCS`, `FiberPath`, `PathIntegral`).
Plotting helpers are opt-in via `using Bifrost.Plots`.

Typical use: author geometry on a `SubpathBuilder` (`start!` → segment calls →
`jumpto!`/`seal!`), bind it with `Fiber(spec; cross_section, T_ref_K)`, then call
`propagate_fiber(fiber; λ_m)` or `propagate_fiber_sensitivity(fiber; λ_m)`.
"""
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

"""
Intrinsic optical and mechanical properties of fiber glasses.

Each concrete glass subtypes `AbstractMaterial` and implements the property interface
documented there (`refractive_index`, `cte`, `softening_temperature`, `poisson_ratio`,
`photoelastic_constants`, `youngs_modulus`, …). Concrete materials: `SiO2`, `GeO2`,
`SilicaGermaniaGlass`, and `SilicaFluorinatedGlass`.

Units are SI throughout: wavelengths in metres, temperatures in kelvin, moduli in
pascals; refractive indices, Poisson ratios, and photoelastic constants are
dimensionless.
"""
module MaterialProperties
    using LinearAlgebra
    using Printf
    include("material/material-properties.jl")

    include("material/silica.jl")
    include("material/germania.jl")
    include("material/silica-germania.jl")
    include("material/silica-fluorinated.jl")

    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

"""
Construct and query three-dimensional smooth space curves (fiber paths).

Authoring uses a small bang-DSL on a mutable `SubpathBuilder`: `start!` → segment
calls (`straight!`, `bend!`, `helix!`, `catenary!`, `jumpby!`) → seal (`jumpto!`
toward a global target, or `seal!` to end at the natural exit). `Subpath` freezes
the authored data; `build` compiles one Subpath to an immutable `SubpathBuilt`, and
`build(::Vector{Subpath})` concatenates independent Subpaths into a `PathBuilt`
under a shared global arc length, checking endpoint conformity between neighbors.

Built paths answer differential-geometry queries along arc length:

    arc_length(seg_or_path)
    arc_length(path, s1, s2)
    curvature(seg_or_path, s)
    geometric_torsion(seg_or_path, s)
    spin_rate(path, s)
    position(path, s)
    tangent(path, s)
    normal(path, s)
    binormal(path, s)
    frame(path, s)
    breakpoints(path)
    sample(path, s_values)
    sample_uniform(path; n)

Material spin is one spec per Subpath, set at `start!` via the `spin_rate` keyword
(constant rate, function of Subpath-local arc length, or `:inherit`); the
accumulated spin phase is continuous across Subpath boundaries.

Per-segment annotations use the `AbstractMeta` vocabulary (`Nickname`, `MCMadd`,
`MCMmul`). The geometry layer interprets only meta naming a segment's own fields
(applied by `build(...; perturb = true)`) and carries any foreign annotation
through untouched — interpretation of foreign meta (such as the fiber layer's
`:T_K` and `:tension`) is a consuming layer's job.
"""
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

"""
Transverse fiber cross-section optics.

Converts material properties into guided-mode quantities (normalized frequency,
propagation constant, dispersion, effective area) and local birefringence response
magnitudes (bending, twist, core ellipticity, asymmetric thermal stress, axial
tension) for a single transverse slice of fiber. Concrete types subtype
`FiberCrossSection`; see `StepIndexCrossSection`.
"""
module FiberCS
    using LinearAlgebra
    using ..MaterialProperties
    include("fiber-cross-section/cross-section.jl")

    include("fiber-cross-section/step-index.jl")
    include("fiber-cross-section/graded-index.jl")

    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

"""
Path-backed fiber assembly.

`Fiber` binds authored path geometry to a `FiberCrossSection` and assembles the
local Jones generators `K(s)` (via `generator_K`) and `Kω(s)` (via `generator_Kω`)
consumed by the propagation layer. This module is the sole interpreter of the
foreign segment meta `:T_K` (temperature excursion, K) and `:tension` (axial
tension, N), which the geometry layer carries inertly.
"""
module FiberPath
    using LinearAlgebra
    using ..MaterialProperties
    using ..PathGeometry
    # Internal cross-module references. The geometry-layer perturbation mechanism
    # (used by the fiber's thermal :T_K interpretation) and a couple of helpers
    # are underscore-prefixed and not exported.
    using ..PathGeometry: _scale_length_fields, _scale_inverse_twist_rate, _meta_without,
                          _length_fields, _qc_nominalize, _resolve_inherited_start
    using ..FiberCS
    include("fiber/fiber-path.jl")
    import ..Bifrost: _export_public!
    _export_public!(@__MODULE__)
end

"""
Lossless Jones-matrix propagation and DGD sensitivity integration.

Advances `dJ/ds = K(s)·J` for a callable local generator `K` with adaptive
step-doubling exponential-midpoint integration, respecting caller-supplied
breakpoints and optional lumped jumps, and integrates the coupled sensitivity
system for `G = ∂ωJ`. Fiber-level entry points are `propagate_fiber` and
`propagate_fiber_sensitivity`; DGD extraction uses `output_dgd` and the
MCM-friendly `output_dgd_2x2`.

The implementation assumes lossless SU(2) Jones dynamics. Error control is
phase-insensitive, and MCM-compatible code paths avoid scalar coercions,
particle-dependent branching, and generic matrix exponentials that do not lift
through `MonteCarloMeasurements.Particles`.
"""
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
"""
Visual diagnostics for built path geometry and fiber propagation.

Generates interactive Plotly HTML renderings: 3D centerlines with frame cursors,
polarization evolution along a fiber, Poincaré-sphere views, and adaptive-step
solver diagnostics. Opt-in via `using Bifrost.Plots`.
"""
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
