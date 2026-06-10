# Material-agnostic perturbation mechanism for the geometry layer:
# `build(sub; perturb=true)` applies field-level MCM via `_apply_field_mcm`,
# and `_scale_length_fields` is the reusable isotropic-scaling transform that
# consuming layers (the fiber assembly) drive. Nothing here references
# temperature or material properties.
#
# Segment types (StraightSegment, …, JumpBy), the meta vocabulary (MCMcombine,
# segment_meta), and AbstractPathSegment are all defined in path-geometry.jl /
# path-geometry-meta.jl, which are in scope when this file is included.

"""
    _length_fields(seg) -> Tuple{Vararg{Symbol}}

The length-dimensioned fields of `seg` that scale with an isotropic factor `τ`.
Dimensionless fields (angles, turns) are omitted. The inverse-length `twist`
rate is handled separately by [`_scale_length_fields`](@ref) (divided, not
multiplied). Each concrete `AbstractPathSegment` must declare this; the fallback
errors so a newly added segment type cannot be silently mis-scaled.
"""
_length_fields(::StraightSegment) = (:length,)
_length_fields(::BendSegment)     = (:radius,)
_length_fields(::CatenarySegment) = (:a, :length)
_length_fields(::HelixSegment)    = (:radius, :pitch)
_length_fields(seg::AbstractPathSegment) = error(
    "perturb: $(typeof(seg)) does not declare _length_fields; add a method listing " *
    "its length-dimensioned fields so isotropic scaling is well-defined.")

# Reconstruct a segment of the same type from its non-`meta` field values (in
# declaration order) plus a meta vector, using the *outer* (UnionAll) constructor
# so the element type re-promotes — this lifts a `Float64` segment to `Particles`
# when a scaled value is `Particles`.
_reconstruct_segment(seg::AbstractPathSegment, vals, new_meta) =
    (Base.typename(typeof(seg)).wrapper)(vals...; meta = new_meta)

# Scale a `TwistRate` (rad/m) inversely by `factor`. Mechanical twist is an
# inverse-length rate: under an isotropic length scaling by `factor`, the total
# turns `∫τ_m ds` are conserved, so the rate divides by `factor`. A function rate
# `τ_m(s_local)` is also reparametrized onto the stretched local arc length:
# `g(s) = τ_m(s/factor)/factor`, so `∫₀^{factor·L} g = ∫₀^L τ_m` (turns conserved).
_scale_inverse_twist_rate(::Nothing, _factor) = nothing
_scale_inverse_twist_rate(rate::Real, factor) = rate / factor
_scale_inverse_twist_rate(rate::Function, factor) = s -> rate(s / factor) / factor

"""
    _scale_length_fields(seg, factor, new_meta) -> AbstractPathSegment

Return a copy of `seg` with every length-dimensioned field (see
[`_length_fields`](@ref)) multiplied by `factor`, the inverse-length `twist` rate
divided by `factor` (conserving total turns; see [`_scale_inverse_twist_rate`](@ref)),
all other fields preserved, and `meta` replaced by `new_meta`. Broadcasting
(`.* factor`) covers both scalar and tuple-valued length fields. Material- and
temperature-agnostic: the caller chooses `factor`.
"""
function _scale_length_fields(seg::AbstractPathSegment, factor, new_meta)
    lf   = _length_fields(seg)
    keep = filter(!=(:meta), fieldnames(typeof(seg)))
    vals = map(keep) do f
        if f === :twist
            # inverse-length rate: divide so total turns are conserved
            _scale_inverse_twist_rate(getfield(seg, f), factor)
        elseif f in lf
            # length-dimensioned field: multiply
            getfield(seg, f) .* factor
        else
            getfield(seg, f)
        end
    end
    return _reconstruct_segment(seg, vals, new_meta)
end

# JumpBy is resolved to a QuinticConnector at placement time and has a non-uniform
# constructor (only `delta` is positional), so it does not participate in field
# scaling here: pass it through, applying only the requested meta replacement. Its
# `twist` is also left unscaled — its length is not scaled here, so its turns are
# not redistributed.
_scale_length_fields(seg::JumpBy, _factor, new_meta) =
    JumpBy(seg.delta; tangent_out = seg.tangent_out, curvature_out = seg.curvature_out,
           min_bend_radius = seg.min_bend_radius, twist = seg.twist, meta = new_meta)

"""
    _apply_field_mcm(seg) -> AbstractPathSegment

Apply every `MCMadd`/`MCMmul` whose symbol names one of `seg`'s own fields to that
field via [`MCMcombine`](@ref), returning a reconstructed segment. Meta whose
symbol is not a field of `seg` is ignored (carried through unchanged), so foreign
annotations are inert here. Segments without meta are returned unchanged.
"""
function _apply_field_mcm(seg::AbstractPathSegment)
    isempty(segment_meta(seg)) && return seg
    keep = filter(!=(:meta), fieldnames(typeof(seg)))
    vals = map(f -> MCMcombine(getfield(seg, f), seg, f), keep)
    return _reconstruct_segment(seg, vals, seg.meta)
end

# JumpBy carries no field-level MCM and is resolved at placement; pass through.
_apply_field_mcm(seg::JumpBy) = seg
