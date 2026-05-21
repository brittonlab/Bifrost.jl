"""
fiber-path-meta.jl

Concrete `AbstractMeta` subtypes for per-segment annotations. `path-geometry.jl`
defines only the abstract slot (`AbstractMeta`) and the `meta` field on each
segment; the concrete vocabulary lives here so path-geometry stays ignorant of
what annotations mean.

## Vocabulary

- `Nickname(label)` — a human-readable label for the segment (used by plotting
  helpers).
- `MCMadd(symbol, distribution)` — an **additive** MCM perturbation:
  consumers apply `baseline + sample` (e.g. `T_K = T_ref + ΔT`).
- `MCMmul(symbol, distribution)` — a **multiplicative** MCM perturbation:
  consumers apply `baseline * sample` (direct scale factor, so e.g.
  `MCMmul(:length, 0.5)` halves `length`; `MCMmul(:length, -0.4)` flips the
  sign and shortens).

`symbol` is a `Symbol` (e.g. `:T_K`). `distribution` is any object a consumer
knows how to sample from (scalar, `Particles`, `Distributions.*`). The
additive vs. multiplicative distinction is encoded in the type itself, so
consumers dispatch on it rather than looking up a mode table.

This file contains no sampling or interpretation logic — that belongs with
whichever layer acts on the annotation.
"""

struct Nickname <: AbstractMeta
    label::String
end

struct MCMadd{D} <: AbstractMeta
    symbol::Symbol
    distribution::D
end

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
    MCMcombine(baseline, seg, sym::Symbol) → perturbed

Combine every `MCMmul(sym, d)` and `MCMadd(sym, d)` on segment `seg` and
apply them to `baseline` in a uniform order:

    perturbed = baseline * Π(d_mul) + Σ d_add

All multiplicative factors are applied first (as direct scale factors), then
all additive offsets are summed on top. Non-matching-symbol entries are
ignored. When no matching entries are present, `baseline` is returned
unchanged (same object, no arithmetic performed).

Intended as the shared composition helper that consumers (shrinkage, future
generator passes) call so every layer applies MCM perturbations the same way.
"""
function MCMcombine(baseline, seg, sym::Symbol)
    mul = nothing
    add = nothing
    for m in segment_meta(seg)
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
