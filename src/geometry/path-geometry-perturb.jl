"""
path-geometry-perturb.jl

Material-agnostic perturbation mechanism for the geometry layer.

`build(sub; perturb=true)` applies the meta the geometry layer can interpret on
its own — `MCMadd`/`MCMmul` whose symbol names one of a segment's *own* struct
fields (`:length`, `:radius`, …) — via [`_apply_field_mcm`](@ref). Any meta this
layer does not recognize (e.g. a foreign thermal annotation) is carried through
untouched; interpretation of foreign meta is the consuming layer's job.

The companion transform [`_scale_length_fields`](@ref) multiplies a segment's
length-dimensioned fields by a scalar factor. It is the reusable mechanism that a
consuming layer (the fiber assembly) uses to apply isotropic expansion, supplying
both the factor and the replacement meta. Nothing here references temperature or
material properties.

Which fields are length-dimensioned is declared per segment type by
[`_length_fields`](@ref); a new `AbstractPathSegment` that omits it errors loudly
when perturbed, so the extension point is self-documenting.
"""

# Segment types (StraightSegment, …, JumpBy), the meta vocabulary (MCMcombine,
# segment_meta), and AbstractPathSegment are all defined in path-geometry.jl /
# path-geometry-meta.jl, which are in scope when this file is included.

"""
    _length_fields(seg) -> Tuple{Vararg{Symbol}}

The length-dimensioned fields of `seg` that scale with an isotropic factor `τ`.
Dimensionless fields (angles, turns) are omitted. Each concrete
`AbstractPathSegment` must declare this; the fallback errors so a newly added
segment type cannot be silently mis-scaled.
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

"""
    _scale_length_fields(seg, factor, new_meta) -> AbstractPathSegment

Return a copy of `seg` with every length-dimensioned field (see
[`_length_fields`](@ref)) multiplied by `factor`, all other fields preserved, and
`meta` replaced by `new_meta`. Broadcasting (`.* factor`) covers both scalar
fields and tuple-valued fields. Material- and temperature-agnostic: the caller
chooses `factor`.
"""
function _scale_length_fields(seg::AbstractPathSegment, factor, new_meta)
    lf   = _length_fields(seg)
    keep = filter(!=(:meta), fieldnames(typeof(seg)))
    vals = map(f -> f in lf ? getfield(seg, f) .* factor : getfield(seg, f), keep)
    return _reconstruct_segment(seg, vals, new_meta)
end

# JumpBy is resolved to a QuinticConnector at placement time and has a non-uniform
# constructor (only `delta` is positional), so it does not participate in field
# scaling here: pass it through, applying only the requested meta replacement.
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
