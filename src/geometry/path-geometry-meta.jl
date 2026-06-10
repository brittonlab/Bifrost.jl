# Concrete `AbstractMeta` subtypes (`Nickname`, `MCMadd`, `MCMmul`) for
# per-segment and per-Subpath annotations; path-geometry.jl defines only the
# abstract slot and the `meta` fields. This file contains no sampling or
# interpretation logic — that belongs with whichever layer acts on the
# annotation; foreign meta is carried through inertly (see the
# Bifrost.PathGeometry module docstring in src/Bifrost.jl).
#
# `AbstractMeta` is defined in path-geometry.jl, which includes this file. When
# loaded as part of the Bifrost package, path-geometry.jl is already in scope.

"""
    Nickname(label)

Attach a human-readable label to a path segment or Subpath.

Plotting and demo code use this metadata for display labels. The geometry
builder stores it in the `meta` bag and otherwise leaves interpretation to
consuming layers.
"""
struct Nickname <: AbstractMeta
    label::String
end

"""
    MCMadd(symbol, distribution)

Attach an additive Monte Carlo perturbation to `symbol`.

Consumers combine matching entries as `baseline * product(mul) + sum(add)`;
`MCMadd` contributes to the additive sum. `distribution` may be any object the
consumer knows how to sample or apply, including a scalar or `Particles`.
"""
struct MCMadd{D} <: AbstractMeta
    symbol::Symbol
    distribution::D
end

"""
    MCMmul(symbol, distribution)

Attach a multiplicative Monte Carlo perturbation to `symbol`.

Consumers combine matching entries as `baseline * product(mul) + sum(add)`;
`MCMmul` contributes a direct scale factor to the multiplicative product (so
`MCMmul(:length, 0.5)` halves `length`, and `MCMmul(:length, -0.4)` flips the
sign and shortens). The additive vs. multiplicative distinction is encoded in
the type itself, so consumers dispatch on it rather than looking up a mode
table.
"""
struct MCMmul{D} <: AbstractMeta
    symbol::Symbol
    distribution::D
end

"""
    segment_nickname(seg) → Union{Nothing,String}

Return the first `Nickname` label attached to `seg` via its meta vector, or
`nothing` if none is present.
"""
function segment_nickname(seg)
    for m in segment_meta(seg)
        m isa Nickname && return m.label
    end
    return nothing
end

"""
    MCMcombine(baseline, meta_or_seg, sym::Symbol) → perturbed

Combine every `MCMmul(sym, d)` and `MCMadd(sym, d)` in a meta vector (or on a
segment, via its `meta`) and apply them to `baseline` in a uniform order:

    perturbed = baseline * Π(d_mul) + Σ d_add

All multiplicative factors are applied first (as direct scale factors), then
all additive offsets are summed on top. Non-matching-symbol entries are
ignored. When no matching entries are present, `baseline` is returned
unchanged (same object, no arithmetic performed).

Accepts either a meta vector (e.g. a Subpath's terminal-connector meta) or any
segment with a `meta` field. Intended as the shared composition helper that
consumers (shrinkage, thermal, future generator passes) call so every layer
applies MCM perturbations the same way.
"""
function MCMcombine(baseline, meta::AbstractVector{<:AbstractMeta}, sym::Symbol)
    mul = nothing
    add = nothing
    for m in meta
        if m isa MCMmul && m.symbol === sym
            mul = isnothing(mul) ? m.distribution : mul * m.distribution
        elseif m isa MCMadd && m.symbol === sym
            add = isnothing(add) ? m.distribution : add + m.distribution
        end
    end
    isnothing(mul) && isnothing(add) && return baseline
    scaled = isnothing(mul) ? baseline : baseline * mul
    return isnothing(add) ? scaled : scaled + add
end

MCMcombine(baseline, seg, sym::Symbol) = MCMcombine(baseline, segment_meta(seg), sym)

"""
    _meta_without(meta_or_seg, sym::Symbol) -> Vector{AbstractMeta}

Return a fresh meta vector (copied from a meta vector or a segment's `meta`) with
every `MCMadd(sym, …)` and `MCMmul(sym, …)` removed. Other meta (including MCM
entries for different symbols) is preserved in order. Generic: it names no specific
symbol, so a consuming layer can use it to strip an annotation it has just
interpreted (passing, e.g., its own thermal symbol) without the geometry layer
knowing what that symbol means.
"""
function _meta_without(meta::AbstractVector{<:AbstractMeta}, sym::Symbol)
    out = AbstractMeta[]
    for m in meta
        ((m isa MCMadd || m isa MCMmul) && m.symbol === sym) && continue
        push!(out, m)
    end
    return out
end

_meta_without(seg, sym::Symbol) = _meta_without(segment_meta(seg), sym)
