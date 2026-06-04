# =====================================================================
# demo2.jl — JumpBy / JumpTo visual experiments
# =====================================================================
#
# Companion file to demo1.jl, dedicated to illustrating `JumpBy` and
# `JumpTo`. Two flavours:
#
#   * 2D scenes (inline SVG, no external library):
#     demo_jumpby_2d_tangent_out / demo_jumpto_2d_* / demo_jumpto_2d_min_radius
#     — the path centerline projected onto the x–z plane. Useful as a
#     quick-look that needs no JavaScript.
#
#   * 3D scenes (Plotly): demo_jumpby_* / demo_jumpto_* — the same
#     scenes, with the resolved QuinticConnector drawn in red and the
#     surrounding fixed segments drawn in green.
#
# `demo2_all()` runs every demo in `DEMO2_INDEX` and writes
# `output/demo-index.html` with the groups under separate headings.
#
# This file expects to be `include`d after demo1.jl is in scope (it
# reuses `_sample_segment_xyz` and the path-builder API).
#
# Each demo function builds its path inline using PathSpecBuilder — no
# shared builder helpers. Code duplication across variants is intentional
# per AGENTS.md §5 (Visual Tests): clarity over abstraction.

if !isdefined(Main, :_sample_segment_xyz)
    include(joinpath(@__DIR__, "demo1.jl"))
end

# The compatibility proxies below add methods to the SubpathBuilder authoring
# verbs for proxy receiver types. Import the functions from Bifrost so these
# definitions EXTEND the real generics rather than shadow them with new
# `Main`-local functions (which would break `straight!(::SubpathBuilder)`).
import Bifrost: straight!, bend!, helix!, catenary!, jumpby!, jumpto!, build, Fiber

# ---------------------------------------------------------------------
# Compatibility proxy: lets demo closures keep using the old
# `PathSpecBuilder()` + `straight!/bend!/jumpby!/jumpto!` + `build(spec)`
# idiom while producing a SubpathBuilt or PathBuilt under the new API.
#
# Each `jumpto!` call seals the current Subpath and starts a new one at
# the just-sealed endpoint; the proxy trial-builds to copy the actual end
# tangent into the next Subpath's `start_outgoing_tangent` so the
# PathBuilt conformity check passes by construction.
#
# Old kwargs are translated:
#   destination= → point=
#   tangent=     → incoming_tangent=
#   curvature_out= → incoming_curvature=
#
# `meta` on the terminal jumpto is discarded — the new architecture's
# terminal connector is structural and carries no meta.
# ---------------------------------------------------------------------

mutable struct _PathProxy
    subpaths::Vector{SubpathBuilder}
    current::Union{SubpathBuilder, Nothing}
end

function PathSpecBuilder(; spin_rate = nothing)
    sb = SubpathBuilder(); start!(sb; spin_rate = spin_rate)
    return _PathProxy(SubpathBuilder[], sb)
end

# Forward interior-segment builders to the current Subpath.
function straight!(p::_PathProxy; kwargs...)
    @assert !isnothing(p.current) "PathSpecBuilder proxy already finalized"
    straight!(p.current; kwargs...)
    return p
end
function bend!(p::_PathProxy; kwargs...)
    @assert !isnothing(p.current) "PathSpecBuilder proxy already finalized"
    bend!(p.current; kwargs...)
    return p
end
function helix!(p::_PathProxy; kwargs...)
    @assert !isnothing(p.current) "PathSpecBuilder proxy already finalized"
    helix!(p.current; kwargs...)
    return p
end
function catenary!(p::_PathProxy; kwargs...)
    @assert !isnothing(p.current) "PathSpecBuilder proxy already finalized"
    catenary!(p.current; kwargs...)
    return p
end
function jumpby!(p::_PathProxy; kwargs...)
    @assert !isnothing(p.current) "PathSpecBuilder proxy already finalized"
    jumpby!(p.current; kwargs...)
    return p
end

# jumpto! seals the current Subpath and starts a new one.
function jumpto!(p::_PathProxy;
                 destination = nothing,
                 point = nothing,
                 tangent = nothing,
                 incoming_tangent = nothing,
                 curvature_out = nothing,
                 incoming_curvature = nothing,
                 min_bend_radius = nothing,
                 meta = AbstractMeta[])
    @assert !isnothing(p.current) "PathSpecBuilder proxy already finalized"
    pt = isnothing(destination) ? point : destination
    @assert !isnothing(pt) "jumpto!: must provide destination or point"
    inc_tan_user  = isnothing(tangent)       ? incoming_tangent   : tangent
    inc_curv_user = isnothing(curvature_out) ? incoming_curvature : curvature_out

    jumpto!(p.current;
        point                = pt,
        incoming_tangent     = inc_tan_user,
        incoming_curvature   = inc_curv_user,
        min_bend_radius      = min_bend_radius,
        meta                 = meta,
    )

    # Use the user-supplied jumpto_point exactly for the next Subpath's
    # start_point so the conformity check passes despite the connector
    # solver's finite tolerance (~1e-6 rel).
    pt_tup = (Float64(pt[1]), Float64(pt[2]), Float64(pt[3]))

    # For the start tangent: if the user supplied an explicit tangent, use
    # that (it matches the jumpto_incoming_tangent exactly). Otherwise
    # trial-build to read the actual connector exit tangent, which equals
    # the chord direction.
    sealed = Subpath(p.current)
    sealed_built = build(sealed)
    out_tan = if isnothing(inc_tan_user)
        end_tan_global = collect(end_tangent(sealed_built))
        (Float64(end_tan_global[1]),
         Float64(end_tan_global[2]),
         Float64(end_tan_global[3]))
    else
        (Float64(inc_tan_user[1]),
         Float64(inc_tan_user[2]),
         Float64(inc_tan_user[3]))
    end

    push!(p.subpaths, p.current)

    next = SubpathBuilder()
    out_curv = isnothing(inc_curv_user) ?
        (0.0, 0.0, 0.0) :
        (Float64(inc_curv_user[1]),
         Float64(inc_curv_user[2]),
         Float64(inc_curv_user[3]))
    start!(next;
        point              = pt_tup,
        outgoing_tangent   = out_tan,
        outgoing_curvature = out_curv,
    )
    p.current = next
    return p
end

# build(::_PathProxy): seal the trailing Subpath at its natural exit (if
# it has any segments), then assemble. Returns SubpathBuilt for one Subpath
# or PathBuilt for many.
function build(p::_PathProxy)
    @assert !isnothing(p.current) "build: proxy already finalized"
    if !isempty(p.current.segments)
        seal!(p.current)
        push!(p.subpaths, p.current)
    end
    p.current = nothing  # mark finalized

    if isempty(p.subpaths)
        error("build(::_PathProxy): no segments authored")
    end

    subs = [Subpath(sb) for sb in p.subpaths]
    return length(subs) == 1 ? build(subs[1]) : build(subs)
end

# Fiber(::_PathProxy; ...): seal the trailing Subpath (mirroring
# build(::_PathProxy)), then hand the resulting Subpath / Vector{Subpath} to the
# library Fiber constructor so it applies thermal `:T_K` and field-level MCM.
function Fiber(p::_PathProxy; cross_section, T_ref_K)
    @assert !isnothing(p.current) "Fiber: proxy already finalized"
    if !isempty(p.current.segments)
        seal!(p.current)
        push!(p.subpaths, p.current)
    end
    p.current = nothing  # mark finalized

    if isempty(p.subpaths)
        error("Fiber(::_PathProxy): no segments authored")
    end

    subs = [Subpath(sb) for sb in p.subpaths]
    spec = length(subs) == 1 ? subs[1] : subs
    return Fiber(spec; cross_section = cross_section, T_ref_K = T_ref_K)
end

# =====================================================================
# 2D GROUP — inline-SVG scenes
# =====================================================================
#
# Renderer + scenes for the 2D demos (x–z projection). No JavaScript,
# no external dependencies. Scene functions live below the renderer.
# ---------------------------------------------------------------------
# 2D inline-SVG renderer
# ---------------------------------------------------------------------
#
# All paths in this demo set lie in the y = 0 plane (we only ever use
# axis_angle = 0 and never set y-components in delta/destination), so
# the natural 2D projection is (x, z). The SVG is emitted directly —
# no JS, no external dependency.

# Sample every placed segment of `path`. Returns a Vector of
# (xs, zs, is_red) tuples — one per segment.
function _sample_path_2d(path, red_indices::Vector{Int})
    rows = Tuple{Vector{Float64}, Vector{Float64}, Bool}[]
    n = length(_all_placed(path))
    for i in 1:n
        s = _sample_segment_xyz(path, i)
        push!(rows, (s.x, s.z, i in red_indices))
    end
    return rows
end

# Render a single 2D scene to `output` as a standalone SVG-in-HTML file.
# `variants` is a list of (label, build_fn) where build_fn returns
# `(path, red_indices)`. Variants are offset along x by `variant_spacing`.
# Setting variant_spacing=0 overlays them.
function _jump_row_svg(output::AbstractString, title::AbstractString;
                       variants::Vector,
                       variant_spacing::Float64 = 2.5,
                       width::Int = 1100, height::Int = 520,
                       margin::Int = 60)
    # Collect all sampled curves across all variants, applying x-offset.
    all_curves = Vector{Tuple{Vector{Float64}, Vector{Float64}, Bool, Int}}()
    starts     = Vector{Tuple{Float64, Float64}}()
    labels     = Vector{Tuple{Float64, Float64, String}}()

    xmin =  Inf; xmax = -Inf
    zmin =  Inf; zmax = -Inf

    for (k, (label, build_fn)) in enumerate(variants)
        path, red_spec = build_fn()
        red_indices = red_spec isa Integer ? Int[red_spec] :
                      Int.(collect(red_spec))
        dx   = (k - 1) * variant_spacing
        rows = _sample_path_2d(path, red_indices)
        for (xs, zs, is_red) in rows
            xs2 = xs .+ dx
            push!(all_curves, (xs2, zs, is_red, k))
            xmin = min(xmin, minimum(xs2)); xmax = max(xmax, maximum(xs2))
            zmin = min(zmin, minimum(zs));  zmax = max(zmax, maximum(zs))
        end
        p0 = position(path, 0.0)
        push!(starts, (Float64(p0[1]) + dx, Float64(p0[3])))
        # Label sits just above each variant.
        x_label = dx
        push!(labels, (x_label, zmax, label))
    end

    # Padding so the curves don't touch the axes.
    pad_x = 0.05 * max(xmax - xmin, 1e-6)
    pad_z = 0.10 * max(zmax - zmin, 1e-6)
    xmin -= pad_x; xmax += pad_x
    zmin -= pad_z; zmax += pad_z + 0.15  # extra room for labels

    # Map data (x, z) to SVG (px, py). SVG y grows downward, so flip z.
    plot_w = width  - 2*margin
    plot_h = height - 2*margin
    sx = plot_w / (xmax - xmin)
    sz = plot_h / (zmax - zmin)
    s  = min(sx, sz)  # equal aspect
    # Center the data in the plot region.
    ox = margin + 0.5 * (plot_w - s * (xmax - xmin)) - s * xmin
    oz = margin + 0.5 * (plot_h - s * (zmax - zmin)) + s * zmax
    px(x) = ox + s * x
    pz(z) = oz - s * z

    # Build the SVG body.
    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(width)" height="$(height)" viewBox="0 0 $(width) $(height)">""")
    println(io, """<rect width="100%" height="100%" fill="#111"/>""")

    # Axis box.
    bx0 = px(xmin); bz0 = pz(zmax); bw = s*(xmax-xmin); bh = s*(zmax-zmin)
    println(io, """<rect x="$(round(bx0;digits=2))" y="$(round(bz0;digits=2))" width="$(round(bw;digits=2))" height="$(round(bh;digits=2))" fill="none" stroke="#333" stroke-width="1"/>""")

    # Light grid at integer coordinates (data units of metres).
    for ix in ceil(Int, xmin):floor(Int, xmax)
        x = px(Float64(ix))
        println(io, """<line x1="$(round(x;digits=2))" y1="$(round(bz0;digits=2))" x2="$(round(x;digits=2))" y2="$(round(bz0+bh;digits=2))" stroke="#222" stroke-width="0.5"/>""")
    end
    for iz in ceil(Int, zmin):floor(Int, zmax)
        z = pz(Float64(iz))
        println(io, """<line x1="$(round(bx0;digits=2))" y1="$(round(z;digits=2))" x2="$(round(bx0+bw;digits=2))" y2="$(round(z;digits=2))" stroke="#222" stroke-width="0.5"/>""")
    end

    # Axis labels.
    println(io, """<text x="$(round(bx0+bw/2;digits=2))" y="$(height-15)" fill="#aaa" font-family="sans-serif" font-size="13" text-anchor="middle">x (m)</text>""")
    println(io, """<text x="18" y="$(round(bz0+bh/2;digits=2))" fill="#aaa" font-family="sans-serif" font-size="13" text-anchor="middle" transform="rotate(-90 18,$(round(bz0+bh/2;digits=2)))">z (m)</text>""")
    println(io, """<text x="$(round(width/2;digits=2))" y="22" fill="#eee" font-family="sans-serif" font-size="15" text-anchor="middle">$(title)</text>""")

    # Curves.
    for (xs, zs, is_red, _k) in all_curves
        color = is_red ? "#e55" : "#6c6"
        # Build a polyline points string.
        parts = String[]
        for i in eachindex(xs)
            push!(parts, string(round(px(xs[i]); digits=2), ",", round(pz(zs[i]); digits=2)))
        end
        println(io, """<polyline points="$(join(parts, ' '))" fill="none" stroke="$(color)" stroke-width="2.5" stroke-linejoin="round"/>""")
    end

    # Start markers (green dots).
    for (sxv, szv) in starts
        println(io, """<circle cx="$(round(px(sxv);digits=2))" cy="$(round(pz(szv);digits=2))" r="4" fill="#6c6" stroke="#9d9" stroke-width="1"/>""")
    end

    # White circles at the end of every atomic segment (segment boundaries).
    for (xs, zs, _is_red, _k) in all_curves
        ex = xs[end]; ez = zs[end]
        println(io, """<circle cx="$(round(px(ex);digits=2))" cy="$(round(pz(ez);digits=2))" r="3" fill="#fff" stroke="#000" stroke-width="0.6"/>""")
    end

    # Variant labels (top of each column).
    for (lx, lz, txt) in labels
        println(io, """<text x="$(round(px(lx);digits=2))" y="$(round(pz(lz)-6;digits=2))" fill="#ddd" font-family="sans-serif" font-size="12" text-anchor="middle">$(txt)</text>""")
    end

    println(io, "</svg>")
    svg_body = String(take!(io))

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>$(title)</title>
<style>html,body{margin:0;padding:0;background:#111;color:#eee;font-family:sans-serif;}
.wrap{display:flex;align-items:center;justify-content:center;min-height:100vh;}
svg{max-width:100%;height:auto;}</style>
</head><body><div class="wrap">$(svg_body)</div></body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end

# Helper: build a path from an explicit list of builder calls (one per
# atomic segment) and return it in the (path, red_indices) shape that
# `_jump_row_svg` expects. Use it like:
#
#     _path() do spec
#         straight!(spec; length = 1.0)
#         jumpby!(spec;   delta  = (0.2, 0.0, 0.5))
#         straight!(spec; length = 1.0)
#     end
#
# `red` is optional and lists which atomic segments to colour red; with
# the default `Int[]` everything draws in the single non-red colour.
function _path(f::Function; red::Vector{Int} = Int[])
    spec = PathSpecBuilder()
    f(spec)
    return (build(spec), red)
end

# Build a path one builder call at a time, trapping the first call whose
# `build()` throws. Returns `(partial_path_or_nothing, n_built, error_or_nothing)`.
# Used by `min_bend_radius` demos to display a partial path with the
# infeasible step missing.
#
# Each `f(probe)` invocation should call atomic builders. We replay the call
# log step-by-step, building a fresh proxy each iteration and stopping at
# the first failure.
function _build_with_failure(f::Function)
    # Record builder calls via a recorder proxy.
    recorded = Vector{Tuple{Symbol, Tuple, Dict{Symbol, Any}}}()
    rec = _RecorderProxy(recorded)
    f(rec)

    n = length(recorded)
    last_good = nothing
    last_n = 0
    for i in 1:n
        # Replay the first i calls.
        spec = PathSpecBuilder()
        ok = true
        local err
        for j in 1:i
            (sym, args, kw) = recorded[j]
            try
                _replay_call(spec, sym, args, kw)
            catch e
                ok = false
                err = e
                break
            end
        end
        if !ok
            return (last_good, last_n, err)
        end
        try
            last_good = build(spec)
            last_n = i
        catch e
            return (last_good, last_n, e)
        end
    end
    return (last_good, n, nothing)
end

# Recorder proxy: captures (symbol, args, kwargs) for each builder call.
mutable struct _RecorderProxy
    log::Vector{Tuple{Symbol, Tuple, Dict{Symbol, Any}}}
end
straight!(p::_RecorderProxy; kwargs...) = (push!(p.log, (:straight, (), Dict{Symbol,Any}(kwargs))); p)
bend!(p::_RecorderProxy; kwargs...)     = (push!(p.log, (:bend,     (), Dict{Symbol,Any}(kwargs))); p)
helix!(p::_RecorderProxy; kwargs...)    = (push!(p.log, (:helix,    (), Dict{Symbol,Any}(kwargs))); p)
catenary!(p::_RecorderProxy; kwargs...) = (push!(p.log, (:catenary, (), Dict{Symbol,Any}(kwargs))); p)
jumpby!(p::_RecorderProxy; kwargs...)   = (push!(p.log, (:jumpby,   (), Dict{Symbol,Any}(kwargs))); p)
jumpto!(p::_RecorderProxy; kwargs...)   = (push!(p.log, (:jumpto,   (), Dict{Symbol,Any}(kwargs))); p)

function _replay_call(spec::_PathProxy, sym::Symbol, args::Tuple, kw::Dict{Symbol, Any})
    f = sym === :straight ? straight! :
        sym === :bend     ? bend! :
        sym === :helix    ? helix! :
        sym === :catenary ? catenary! :
        sym === :jumpby   ? jumpby! :
        sym === :jumpto   ? jumpto! :
        error("unknown builder symbol $sym")
    f(spec; kw...)
    return spec
end

# =====================================================================
# 2D demos — individual functions, one per HTML file
# =====================================================================

function demo_jumpby_2d_tangent_out(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_svg(
        joinpath(output_dir, "jumpby-2d-tangent-out.html"),
        "2D — JumpBy: outgoing tangent (local frame)";
        variants = [
            ("t = (+1,0,1)/√2", () -> _path() do spec
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length  = 1.0)
            end),
            ("t = (0,0,1)", () -> _path() do spec
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (0.0, 0.0, 1.0))
                straight!(spec; length  = 1.0)
            end),
            ("t = (-1,0,1)/√2", () -> _path() do spec
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (-1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length  = 1.0)
            end),
        ],
        variant_spacing = 2.0,
    )
end

function demo_jumpto_2d_destination(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-destination.html"),
        "2D — JumpTo(destination=(x, 0, 1.5)): transverse sweep";
        variants = [
            ("x = 0.0", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.0, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
            ("x = 0.3", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.3, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
            ("x = 0.6", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.6, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
            ("x = 1.0", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (1.0, 0.0, 1.5))
                straight!(spec; length      = 1.0)
            end),
        ],
        variant_spacing = 2.5,
    )
end

function demo_jumpto_2d_tangent_global(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-tangent-global.html"),
        "2D — JumpTo: outgoing tangent (GLOBAL frame)";
        variants = [
            ("t = (+1,0,0)", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (1.0, 0.0, 0.0))
                straight!(spec; length      = 1.0)
            end),
            ("t = (0,0,+1)", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (0.0, 0.0, 1.0))
                straight!(spec; length      = 1.0)
            end),
            ("t = (-1,0,1)/√2", () -> _path() do spec
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (-1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length      = 1.0)
            end),
        ],
        variant_spacing = 2.5,
    )
end

function demo_jumpto_2d_routing(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-routing.html"),
        "2D — JumpTo routing: T1/T2/T3 composite (anti-parallel tangents)";
        variants = [
            ("composite", () -> _path(red = [2, 4, 6]) do spec
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (1.0, 0.0, 1.0),
                                tangent         = (0.0, 0.0, -1.0),
                                min_bend_radius = 0.1)
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (2.0, 0.0, 0.0),
                                tangent         = (0.0, 0.0, 1.0),
                                min_bend_radius = 0.1)
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (3.0, 0.0, 1.0),
                                tangent         = (0.0, 0.0, -1.0),
                                min_bend_radius = 0.1)
            end),
        ],
        variant_spacing = 0.0,
        height = 420,
    )
end

"""
    demo_jumpto_2d_min_radius(; output_dir = …)

Sweep `min_bend_radius` on a `JumpTo` whose chord is transverse and whose
outgoing tangent is anti-parallel to the incoming one — the canonical
infeasibility case from the unit tests (threshold ≈ 0.5 m).

Small values succeed; values past the threshold intentionally fail. The
`ArgumentError` is trapped per-variant and the SVG renders the partial
path with the failed jump (and any subsequent segments) omitted; the
variant label shows `(infeasible)`.
"""
function demo_jumpto_2d_min_radius(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    variants = []
    for mbr in (0.10, 0.30, 0.49, 0.51, 0.70)
        partial, n_built, err = _build_with_failure() do spec
            straight!(spec; length = 1.0)
            jumpto!(spec;   destination     = (1.0, 0.0, 1.0),
                            tangent         = (0.0, 0.0, -1.0),
                            min_bend_radius = mbr)
            straight!(spec; length = 1.0)
        end
        partial === nothing && continue
        red   = length(_all_placed(partial)) >= 2 ? [2] : Int[]
        label = err === nothing ?
                "mbr = $(mbr)" :
                "mbr = $(mbr) (infeasible — $(n_built)/3 built)"
        push!(variants, (label, let p = partial, r = red; () -> (p, r); end))
    end

    return _jump_row_svg(
        joinpath(output_dir, "jumpto-2d-min-radius.html"),
        "2D — JumpTo min_bend_radius sweep (transverse chord, anti-parallel tangents; threshold ≈ 0.5 m)";
        variants = variants,
        variant_spacing = 2.0,
        height = 460,
    )
end

# =====================================================================
# meta and JumpTo interplay (2D)
# =====================================================================
#
# Three demos illustrating how `JumpTo` confines the propagation of
# meta-induced perturbations:
#
#   1. S-curve + JumpBy: no anchor — a `:radius` perturbation on the
#      bends shifts every downstream position; the endpoint drifts.
#   2. S-curve + JumpTo: same s-curve geometry as (1), but the JumpTo
#      destination is a lab-frame invariant. The s-curve interior swings
#      visibly while the connector absorbs the slack and the final
#      endpoint stays pinned.
#   3. Helix + JumpTo + `:T_K`: thermal expansion on the helix on top of
#      a radius perturbation. The connector arc length tracks
#      `τ·baseline` while the destination stays pinned (TD001
#      length-constrained resolve).

const _DEMO2_MODIFY_XS = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),     # pure silica cladding → α_lin = SILICA_CTE
    8.2e-6,
    125e-6,
)
const _DEMO2_T_REF = 297.15

# Sample (x, z) per placed segment. Returns Vector{NamedTuple}, one entry
# per segment with `.x` and `.z` arrays.
function _demo2_path_segments_xz(path)
    segs = _all_placed(path)
    return [_sample_segment_xyz(path, i) for i in eachindex(segs)]
end

# Render baseline (black) and modified (red) overlay on a light background.
# Open-circle markers are drawn at the start and end of every placed
# segment (so consecutive segments share a marker at the boundary). Total
# path lengths are shown in a legend whose position is selected by the
# `legend_position` kwarg (`:top_left` or `:bottom_right`).
function _modify_overlay_svg(output::AbstractString, title::AbstractString;
                              baseline, modified,
                              legend_position::Symbol = :top_left,
                              width::Int = 900, height::Int = 540,
                              margin::Int = 70)
    segs_b = _demo2_path_segments_xz(baseline)
    segs_m = _demo2_path_segments_xz(modified)

    xs_b = reduce(vcat, [s.x for s in segs_b])
    zs_b = reduce(vcat, [s.z for s in segs_b])
    xs_m = reduce(vcat, [s.x for s in segs_m])
    zs_m = reduce(vcat, [s.z for s in segs_m])

    xmin = min(minimum(xs_b), minimum(xs_m))
    xmax = max(maximum(xs_b), maximum(xs_m))
    zmin = min(minimum(zs_b), minimum(zs_m))
    zmax = max(maximum(zs_b), maximum(zs_m))

    pad_x = 0.10 * max(xmax - xmin, 1e-6)
    pad_z = 0.10 * max(zmax - zmin, 1e-6)
    xmin -= pad_x; xmax += pad_x
    zmin -= pad_z; zmax += pad_z

    plot_w = width - 2*margin
    plot_h = height - 2*margin
    sx = plot_w / (xmax - xmin)
    sz = plot_h / (zmax - zmin)
    s  = min(sx, sz)
    ox = margin + 0.5*(plot_w - s*(xmax - xmin)) - s*xmin
    oz = margin + 0.5*(plot_h - s*(zmax - zmin)) + s*zmax
    px(x) = ox + s*x
    pz(z) = oz - s*z

    L_b = path_length(baseline)
    L_m = path_length(modified)

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(width)" height="$(height)" viewBox="0 0 $(width) $(height)">""")
    println(io, """<rect width="100%" height="100%" fill="#fafafa"/>""")

    bx0 = px(xmin); bz0 = pz(zmax); bw = s*(xmax-xmin); bh = s*(zmax-zmin)
    println(io, """<rect x="$(round(bx0;digits=2))" y="$(round(bz0;digits=2))" width="$(round(bw;digits=2))" height="$(round(bh;digits=2))" fill="none" stroke="#999" stroke-width="1"/>""")

    for ix in floor(Int, xmin):ceil(Int, xmax)
        x = px(Float64(ix))
        println(io, """<line x1="$(round(x;digits=2))" y1="$(round(bz0;digits=2))" x2="$(round(x;digits=2))" y2="$(round(bz0+bh;digits=2))" stroke="#e5e5e5" stroke-width="0.5"/>""")
    end
    for iz in floor(Int, zmin):ceil(Int, zmax)
        z = pz(Float64(iz))
        println(io, """<line x1="$(round(bx0;digits=2))" y1="$(round(z;digits=2))" x2="$(round(bx0+bw;digits=2))" y2="$(round(z;digits=2))" stroke="#e5e5e5" stroke-width="0.5"/>""")
    end

    println(io, """<text x="$(round(bx0+bw/2;digits=2))" y="$(height-15)" fill="#444" font-family="sans-serif" font-size="13" text-anchor="middle">x (m)</text>""")
    println(io, """<text x="22" y="$(round(bz0+bh/2;digits=2))" fill="#444" font-family="sans-serif" font-size="13" text-anchor="middle" transform="rotate(-90 22,$(round(bz0+bh/2;digits=2)))">z (m)</text>""")
    println(io, """<text x="$(round(width/2;digits=2))" y="22" fill="#222" font-family="sans-serif" font-size="15" text-anchor="middle">$(title)</text>""")

    poly_str(xs, zs) = join(
        (string(round(px(xs[i]); digits=2), ",", round(pz(zs[i]); digits=2))
         for i in eachindex(xs)),
        ' ')
    println(io, """<polyline points="$(poly_str(xs_b, zs_b))" fill="none" stroke="#000" stroke-width="2.5" stroke-linejoin="round"/>""")
    println(io, """<polyline points="$(poly_str(xs_m, zs_m))" fill="none" stroke="#d22" stroke-width="2.5" stroke-linejoin="round"/>""")

    # Open-circle markers at the start and end of every placed segment.
    # Consecutive segments share a boundary, so they overlay one marker
    # per join; first and last segments contribute the path's true
    # endpoints.
    for (segs, color) in ((segs_b, "#000"), (segs_m, "#d22"))
        for s in segs
            for (x, z) in ((s.x[1], s.z[1]), (s.x[end], s.z[end]))
                cx = round(px(x); digits=2)
                cy = round(pz(z); digits=2)
                println(io, """<circle cx="$(cx)" cy="$(cy)" r="3" fill="#fafafa" stroke="$(color)" stroke-width="1.5"/>""")
            end
        end
    end

    legend_w, legend_h = 220, 50
    if legend_position === :bottom_right
        lx = bx0 + bw - 12 - legend_w + 6
        ly = bz0 + bh - 12 - legend_h + 15
    else
        lx = bx0 + 12
        ly = bz0 + 22
    end
    println(io, """<rect x="$(round(lx-6;digits=2))" y="$(round(ly-15;digits=2))" width="$(legend_w)" height="$(legend_h)" fill="#ffffffcc" stroke="#bbb" stroke-width="0.5"/>""")
    println(io, """<line x1="$(round(lx;digits=2))" y1="$(round(ly;digits=2))" x2="$(round(lx+24;digits=2))" y2="$(round(ly;digits=2))" stroke="#000" stroke-width="2.5"/>""")
    println(io, """<text x="$(round(lx+30;digits=2))" y="$(round(ly+4;digits=2))" fill="#222" font-family="sans-serif" font-size="12">baseline:  L = $(round(L_b;digits=4)) m</text>""")
    println(io, """<line x1="$(round(lx;digits=2))" y1="$(round(ly+18;digits=2))" x2="$(round(lx+24;digits=2))" y2="$(round(ly+18;digits=2))" stroke="#d22" stroke-width="2.5"/>""")
    println(io, """<text x="$(round(lx+30;digits=2))" y="$(round(ly+22;digits=2))" fill="#222" font-family="sans-serif" font-size="12">modified:  L = $(round(L_m;digits=4)) m  (Δ = $(round(100*(L_m/L_b - 1);digits=1))%)</text>""")

    println(io, "</svg>")
    svg_body = String(take!(io))

    html = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>$(title)</title>
<style>html,body{margin:0;padding:0;background:#fafafa;color:#222;font-family:sans-serif;}
.wrap{display:flex;align-items:center;justify-content:center;min-height:100vh;}
svg{max-width:100%;height:auto;}</style>
</head><body><div class="wrap">$(svg_body)</div></body></html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end

# ---------------------------------------------------------------------
# (1) S-curve + JumpBy → no anchor → endpoint drifts
# ---------------------------------------------------------------------

function demo_modify_jumpby_drift_2d(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc_short = "Meta on S-curve + JumpBy: no anchor — downstream drifts"
    desc_long  = "S-curve followed by a JumpBy (no anchor). A `:radius` " *
                 "perturbation on the bends inflates total path length by " *
                 "~10%, but since the JumpBy delta is fiber-relative, the " *
                 "entire downstream trajectory drifts and the endpoint " *
                 "visibly separates from the baseline."

    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = 0.0)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = π)
    straight!(spec; length = 0.3)
    jumpby!(spec; delta = (0.0, 0.0, 0.8))
    straight!(spec; length = 1)
    baseline = build(spec)

    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = 0.0,
          meta = [MCMmul(:radius, 1.25)])
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = π,
          meta = [MCMmul(:radius, 1.25)])
    straight!(spec; length = 0.3)
    jumpby!(spec; delta = (0.0, 0.0, 0.8))
    straight!(spec; length = 1)
    modified = Fiber(spec;
        cross_section = _DEMO2_MODIFY_XS, T_ref_K = _DEMO2_T_REF).path

    out = joinpath(output_dir, "modify-jumpby-drift-2d.html")
    path = _modify_overlay_svg(out, desc_short;
        baseline = baseline, modified = modified,
        legend_position = :bottom_right)
    return (path = path, desc = desc_long)
end

# ---------------------------------------------------------------------
# (2) S-curve + JumpTo → anchor → endpoint pinned
# ---------------------------------------------------------------------

function demo_modify_jumpto_anchor_2d(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc_short = "Meta on S-curve + JumpTo: anchor pins endpoint, slack absorbed by connector"
    desc_long  = "Same s-curve geometry as (1), but the trailing JumpBy " *
                 "is replaced with a JumpTo to a lab-frame waypoint. The " *
                 "same scale of `:radius` perturbation on the bends " *
                 "causes the s-curve interior to swing wide — the " *
                 "connector chord changes — but the endpoint stays " *
                 "pinned at the JumpTo destination."

    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = 0.0)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = π)
    straight!(spec; length = 0.3)
    pre_baseline = build(spec)
    pre_pos      = end_point(pre_baseline)
    dz           = 4.000 - path_length(pre_baseline)
    destination  = (pre_pos[1], pre_pos[2], pre_pos[3] + dz)

    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = 0.0)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = π)
    straight!(spec; length = 0.3)
    jumpto!(spec; destination = destination)
    baseline = build(spec)

    spec = PathSpecBuilder()
    straight!(spec; length = 0.3)
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = 0.0,
          meta = [MCMmul(:radius, 1.5)])
    bend!(spec; radius = 0.5, angle = π/2, axis_angle = π,
          meta = [MCMmul(:radius, 1.5)])
    straight!(spec; length = 0.3)
    jumpto!(spec; destination = destination)
    modified = Fiber(spec;
        cross_section = _DEMO2_MODIFY_XS, T_ref_K = _DEMO2_T_REF).path

    out = joinpath(output_dir, "modify-jumpto-anchor-2d.html")
    path = _modify_overlay_svg(out, desc_short;
        baseline = baseline, modified = modified)
    return (path = path, desc = desc_long)
end

# ---------------------------------------------------------------------
# (3) Helix + JumpTo with :T_K on the helix → length tracks τ·baseline,
#     endpoint pinned
# ---------------------------------------------------------------------

function demo_modify_jumpto_anchor_thermal_2d(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    desc_short = "Meta + :T_K + JumpTo carrying :T_K: connector thermally expands (issue #33)"
    desc_long  = "Lead-in straight with `:T_K`, sealed by a transverse " *
                 "JumpTo that also carries `:T_K`. Per issue #33 the terminal " *
                 "connector thermally expands (its arc length scales by τ) " *
                 "while still landing at the fixed destination, so the whole " *
                 "Subpath grows by τ and the extra length shows up as " *
                 "connector curvature."

    α_lin = cte(_DEMO2_MODIFY_XS.cladding_material, _DEMO2_T_REF)
    # 5% length perturbation. The transverse-chord geometry has plenty of
    # connector slack so this stays feasible.
    ΔT = 0.05 / α_lin
    mdt = [MCMadd(:T_K, ΔT)]

    # Geometry: a 0.5 m straight followed by a transverse JumpTo to
    # (1.0, 0, 0.5) with outgoing tangent (1, 0, 0). The connector is
    # ~1.6 m long for a chord of ~1.0 m → ~60% slack to absorb expansion.
    destination = (1.0, 0.0, 0.5)

    spec = PathSpecBuilder()
    straight!(spec; length = 0.5)
    jumpto!(spec; destination = destination, tangent = (1.0, 0.0, 0.0))
    baseline = build(spec)

    spec = PathSpecBuilder()
    straight!(spec; length = 0.5, meta = mdt)
    jumpto!(spec; destination = destination, tangent = (1.0, 0.0, 0.0), meta = mdt)
    modified = Fiber(spec;
        cross_section = _DEMO2_MODIFY_XS, T_ref_K = _DEMO2_T_REF).path

    out = joinpath(output_dir, "modify-jumpto-anchor-thermal-2d.html")
    path = _modify_overlay_svg(out, desc_short;
        baseline = baseline, modified = modified)
    return (path = path, desc = desc_long)
end

# =====================================================================
# 3D GROUP — Plotly scenes
# =====================================================================
#
# Renderer + scenes for the 3D demos. Output: standalone HTML files
# loading plotly.js from CDN. Scene functions live below the renderer.
# ---------------------------------------------------------------------
# 3D Plotly renderer
# ---------------------------------------------------------------------

function _jump_row_html(output::AbstractString, title::AbstractString;
                        variants::Vector,
                        variant_spacing::Float64 = 2.5,
                        z_label_offset::Float64 = 1.6)
    js_num(x)   = isnan(x) ? "NaN" : string(Float64(x))
    js_arr(xs)  = "[" * join(js_num.(xs), ",") * "]"
    js_strarr(xs) = "[" * join(("\"" * replace(string(x), "\"" => "\\\"") * "\"" for x in xs), ",") * "]"

    trace_strs  = String[]
    label_xs    = Float64[]; label_ys = Float64[]; label_zs = Float64[]
    label_texts = String[]
    start_xs    = Float64[]; start_ys = Float64[]; start_zs = Float64[]

    for (k, (label, build_fn)) in enumerate(variants)
        path, red_spec = build_fn()
        red_indices = red_spec isa Integer ? Int[red_spec] : Int.(collect(red_spec))
        dx = (k - 1) * variant_spacing
        n_segs = length(_all_placed(path))
        for i in 1:n_segs
            s     = _sample_segment_xyz(path, i)
            xs    = s.x .+ dx
            color = i in red_indices ? "#e55" : "#6c6"
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

        push!(label_xs, dx)
        push!(label_ys, 0.0)
        push!(label_zs, z_label_offset)
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
# 3D demos — individual functions, one per HTML file
# =====================================================================

function demo_jumpby_delta_axial(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpby-delta-axial.html"),
        "JumpBy(delta=(0,0,d)) — axial sweep";
        variants = [
            ("d = 0.3", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.0, 0.0, 0.3))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("d = 0.6", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.0, 0.0, 0.6))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("d = 1.0", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.0, 0.0, 1.0))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 2.2,
    )
end

function demo_jumpby_delta_transverse(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpby-delta-transverse.html"),
        "JumpBy(delta=(d,0,0.5)) — transverse sweep";
        variants = [
            ("d = 0.0", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.0, 0.0, 0.5))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("d = 0.2", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.2, 0.0, 0.5))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("d = 0.5", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta  = (0.5, 0.0, 0.5))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 1.8,
    )
end

function demo_jumpby_tangent_out(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpby-tangent-out.html"),
        "JumpBy — outgoing tangent (local frame)";
        variants = [
            ("tangent = (+1,0,1)/√2", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length  = 1.0)
                (build(spec), 2)
            end),
            ("tangent = (0,0,1)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (0.0, 0.0, 1.0))
                straight!(spec; length  = 1.0)
                (build(spec), 2)
            end),
            ("tangent = (-1,0,1)/√2", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length  = 1.0)
                jumpby!(spec;   delta   = (0.4, 0.0, 0.4),
                                tangent = (-1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length  = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 1.8,
    )
end

function demo_jumpby_curvature_out(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpby-curvature-out.html"),
        "JumpBy — outgoing curvature (G2 knob, local frame)";
        variants = [
            ("κ_out = 0", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta         = (0.5, 0.0, 0.5),
                                tangent       = (1.0, 0.0, 1.0) ./ sqrt(2),
                                curvature_out = (0.0, 0.0, 0.0))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("κ_out = (0,+2,0)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta         = (0.5, 0.0, 0.5),
                                tangent       = (1.0, 0.0, 1.0) ./ sqrt(2),
                                curvature_out = (0.0, 2.0, 0.0))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("κ_out = (0,-2,0)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length = 1.0)
                jumpby!(spec;   delta         = (0.5, 0.0, 0.5),
                                tangent       = (1.0, 0.0, 1.0) ./ sqrt(2),
                                curvature_out = (0.0, -2.0, 0.0))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 1.8,
    )
end

function demo_jumpto_destination(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpto-destination.html"),
        "JumpTo(destination=(0,0,z)) — axial sweep";
        variants = [
            ("z = 1.3", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.0, 0.0, 1.3))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("z = 1.6", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.0, 0.0, 1.6))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("z = 2.0", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.0, 0.0, 2.0))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 2.6,
    )
end

function demo_jumpto_destination_transverse(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpto-destination-transverse.html"),
        "JumpTo(destination=(x,0,1.5)) — transverse sweep";
        variants = [
            ("x = 0.0", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.0, 0.0, 1.5))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("x = 0.3", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.3, 0.0, 1.5))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("x = 0.6", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.6, 0.0, 1.5))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 2.4,
    )
end

function demo_jumpto_tangent_global(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpto-tangent-global.html"),
        "JumpTo — outgoing tangent (GLOBAL frame)";
        variants = [
            ("tangent = (+1,0,0)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (1.0, 0.0, 0.0))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("tangent = (0,0,+1)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (0.0, 0.0, 1.0))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("tangent = (-1,0,1)/√2", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (-1/sqrt(2), 0.0, 1/sqrt(2)))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
        ],
        z_label_offset = 2.4,
    )
end

function demo_jumpto_curvature_global(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpto-curvature-global.html"),
        "JumpTo — outgoing curvature (GLOBAL frame)";
        variants = [
            ("κ_out = 0", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (0.0, 0.0, 1.0),
                                curvature_out = (0.0, 0.0, 0.0))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("κ_out = (10,0,0)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (0.0, 0.0, 1.0),
                                curvature_out = (10.0, 0.0, 0.0))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
            ("κ_out = (-10,0,0)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                jumpto!(spec;   destination = (0.5, 0.0, 1.5),
                                tangent     = (0.0, 0.0, 1.0),
                                curvature_out = (-10.0, 0.0, 0.0))
                straight!(spec; length      = 1.0)
                (build(spec), 2)
            end),
        ],
        variant_spacing = 1.5,
        z_label_offset = 2.4,
    )
end

function demo_jumpby_after_bend(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpby-after-bend.html"),
        "JumpBy after 90° bend — delta in ROTATED local frame";
        variants = [
            ("delta = (0,0,0.5)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length     = 1.0)
                bend!(spec;     radius     = 0.5, angle = π/2, axis_angle = 0.0)
                jumpby!(spec;   delta      = (0.0, 0.0, 0.5))
                straight!(spec; length     = 1.0)
                (build(spec), 3)
            end),
            ("delta = (0.3,0,0.5)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length     = 1.0)
                bend!(spec;     radius     = 0.5, angle = π/2, axis_angle = 0.0)
                jumpby!(spec;   delta      = (0.3, 0.0, 0.5))
                straight!(spec; length     = 1.0)
                (build(spec), 3)
            end),
            ("delta = (-0.3,0,0.5)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length     = 1.0)
                bend!(spec;     radius     = 0.5, angle = π/2, axis_angle = 0.0)
                jumpby!(spec;   delta      = (-0.3, 0.0, 0.5))
                straight!(spec; length     = 1.0)
                (build(spec), 3)
            end),
        ],
        variant_spacing = 3.5,
        z_label_offset = 2.0,
    )
end

function demo_jumpto_after_bend(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpto-after-bend.html"),
        "JumpTo after 90° bend — destination in GLOBAL frame";
        variants = [
            ("dest = (1.0, 0, 1.5)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                bend!(spec;     radius      = 0.5, angle = π/2, axis_angle = 0.0)
                jumpto!(spec;   destination = (1.0, 0.0, 1.5))
                straight!(spec; length      = 1.0)
                (build(spec), 3)
            end),
            ("dest = (1.3, 0, 1.5)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                bend!(spec;     radius      = 0.5, angle = π/2, axis_angle = 0.0)
                jumpto!(spec;   destination = (1.3, 0.0, 1.5))
                straight!(spec; length      = 1.0)
                (build(spec), 3)
            end),
            ("dest = (0.5, 0, 2.0)", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length      = 1.0)
                bend!(spec;     radius      = 0.5, angle = π/2, axis_angle = 0.0)
                jumpto!(spec;   destination = (0.5, 0.0, 2.0))
                straight!(spec; length      = 1.0)
                (build(spec), 3)
            end),
        ],
        variant_spacing = 3.5,
        z_label_offset = 2.4,
    )
end

function demo_jumpby_g2_inheritance(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpby-g2-inheritance.html"),
        "JumpBy after bend — G2 inheritance of incoming κ from prior bend";
        variants = [
            ("R_bend = 0.30 (κ_in ≈ 3.33)", () -> begin
                spec = PathSpecBuilder()
                bend!(spec;   radius = 0.30, angle = π/4, axis_angle = 0.0)
                jumpby!(spec; delta  = (0.3, 0.0, 0.3))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("R_bend = 0.50 (κ_in = 2.00)", () -> begin
                spec = PathSpecBuilder()
                bend!(spec;   radius = 0.50, angle = π/4, axis_angle = 0.0)
                jumpby!(spec; delta  = (0.3, 0.0, 0.3))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
            ("R_bend = 0.80 (κ_in = 1.25)", () -> begin
                spec = PathSpecBuilder()
                bend!(spec;   radius = 0.80, angle = π/4, axis_angle = 0.0)
                jumpby!(spec; delta  = (0.3, 0.0, 0.3))
                straight!(spec; length = 1.0)
                (build(spec), 2)
            end),
        ],
        variant_spacing = 3.0,
        z_label_offset = 1.6,
    )
end

function demo_jumpto_routing(;
    output_dir::AbstractString = joinpath(@__DIR__, "..", "..", "output"),
)
    return _jump_row_html(
        joinpath(output_dir, "jumpto-routing.html"),
        "JumpTo routing — T1/T2/T3 composite (anti-parallel tangents)";
        variants = [
            ("composite", () -> begin
                spec = PathSpecBuilder()
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (1.0, 0.0, 1.0),
                                tangent         = (0.0, 0.0, -1.0),
                                min_bend_radius = 0.1)
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (2.0, 0.0, 0.0),
                                tangent         = (0.0, 0.0, 1.0),
                                min_bend_radius = 0.1)
                straight!(spec; length          = 1.0)
                jumpto!(spec;   destination     = (3.0, 0.0, 1.0),
                                tangent         = (0.0, 0.0, -1.0),
                                min_bend_radius = 0.1)
                (build(spec), [2, 4, 6])
            end),
        ],
        variant_spacing = 0.0,
        z_label_offset = 1.6,
    )
end

# ---------------------------------------------------------------------
# Monolithic index entries
# ---------------------------------------------------------------------

const DEMO2_INDEX = [
    (group = "2D", fn = demo_jumpby_2d_tangent_out, kwargs = (;),
     desc = "2D JumpBy outgoing tangent sweep (local frame), rendered as inline SVG."),
    (group = "2D", fn = demo_jumpto_2d_destination, kwargs = (;),
     desc = "2D JumpTo transverse destination sweep, rendered as inline SVG."),
    (group = "2D", fn = demo_jumpto_2d_tangent_global, kwargs = (;),
     desc = "2D JumpTo outgoing tangent sweep (global frame), rendered as inline SVG."),
    (group = "2D", fn = demo_jumpto_2d_routing, kwargs = (;),
     desc = "2D JumpTo T1/T2/T3 composite routing with anti-parallel tangents, rendered as inline SVG."),
    (group = "2D", fn = demo_jumpto_2d_min_radius, kwargs = (;),
     desc = "2D sweep of the JumpTo min_bend_radius parameter through its feasibility threshold (~0.5 m). Infeasible variants are trapped at build time and rendered with the missing segment omitted; the variant label is annotated with '(infeasible)'."),
    # These three demos provide their `desc` inline (returned in the
    # NamedTuple alongside the html path), so the narrative lives next
    # to the implementation it describes.
    (group = "2D-modify", fn = demo_modify_jumpby_drift_2d,          kwargs = (;)),
    (group = "2D-modify", fn = demo_modify_jumpto_anchor_2d,         kwargs = (;)),
    (group = "2D-modify", fn = demo_modify_jumpto_anchor_thermal_2d, kwargs = (;)),
    (group = "3D", fn = demo_jumpby_delta_axial, kwargs = (;),
     desc = "3D Plotly: JumpBy axial delta sweep. Connector in red, fixed segments in green."),
    (group = "3D", fn = demo_jumpby_delta_transverse, kwargs = (;),
     desc = "3D Plotly: JumpBy transverse delta sweep. Connector in red, fixed segments in green."),
    (group = "3D", fn = demo_jumpby_tangent_out, kwargs = (;),
     desc = "3D Plotly: JumpBy outgoing tangent sweep (local frame)."),
    (group = "3D", fn = demo_jumpby_curvature_out, kwargs = (;),
     desc = "3D Plotly: JumpBy outgoing curvature sweep (G2 knob, local frame)."),
    (group = "3D", fn = demo_jumpto_destination, kwargs = (;),
     desc = "3D Plotly: JumpTo axial destination sweep."),
    (group = "3D", fn = demo_jumpto_destination_transverse, kwargs = (;),
     desc = "3D Plotly: JumpTo transverse destination sweep."),
    (group = "3D", fn = demo_jumpto_tangent_global, kwargs = (;),
     desc = "3D Plotly: JumpTo outgoing tangent sweep (global frame)."),
    (group = "3D", fn = demo_jumpto_curvature_global, kwargs = (;),
     desc = "3D Plotly: JumpTo outgoing curvature sweep (global frame)."),
    (group = "3D", fn = demo_jumpby_after_bend, kwargs = (;),
     desc = "3D Plotly: JumpBy after 90° bend — delta expressed in the rotated local frame."),
    (group = "3D", fn = demo_jumpto_after_bend, kwargs = (;),
     desc = "3D Plotly: JumpTo after 90° bend — destination in the global frame."),
    (group = "3D", fn = demo_jumpby_g2_inheritance, kwargs = (;),
     desc = "3D Plotly: JumpBy G2 inheritance of incoming curvature from a prior bend."),
    (group = "3D", fn = demo_jumpto_routing, kwargs = (;),
     desc = "3D Plotly: JumpTo T1/T2/T3 composite routing with anti-parallel tangents."),
]

"""
    demo2_all(; index_output)

Run every demo in `DEMO2_INDEX` and write `demo-index.html` linking each
output file with a short description.
"""
const _DEMO2_GROUP_TITLES = Dict(
    "2D"        => "2D scenes (inline SVG)",
    "2D-modify" => "meta and JumpTo interplay (2D)",
    "3D"        => "3D scenes (Plotly)",
)

function demo2_entries()
    entries = Tuple{String, String, String, String}[]

    for d in DEMO2_INDEX
        println("[ demo2 ] $(d.fn)")
        result = d.fn(; d.kwargs...)
        desc = _demo_result_desc(result, d)
        for path in _demo_html_paths(result)
            push!(entries, (d.group, basename(path), path, desc))
        end
    end
    return entries
end

function demo2_all(; index_output::AbstractString = DEMO_MONOLITHIC_INDEX_OUTPUT)
    return _write_demo_index(
        [(title = "JumpBy / JumpTo demos",
          entries = demo2_entries(),
          group_titles = _DEMO2_GROUP_TITLES)];
        index_output,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo2_all()
end
