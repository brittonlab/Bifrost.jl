# FiberCrossSection (from FiberCS) and the Subpath geometry API (from
# PathGeometry) are brought into scope by the FiberPath submodule in Bifrost.jl.

"""
Fiber assembly on top of `path-geometry.jl`.

High-level authoring happens in `path-geometry.jl`:
- build a Subpath spec with `SubpathBuilder` (sealed by `start!` and `jumpto!`)
- bind it to a cross section with `Fiber(sb; cross_section, T_ref_K)` (a
  `SubpathBuilder`, `Subpath`, or `Vector{Subpath}`); the `Fiber` constructor
  builds the geometry once, interpreting thermal `:T_K` meta and applying
  field-level `MCMadd`/`MCMmul`

`Fiber` is the compiled query object consumed downstream by `path-integral.jl`.
It owns:
- the immutable built `SubpathBuilt` or `PathBuilt`
- the `FiberCrossSection`
- a single reference temperature `T_ref_K` (reference for path geometry and
  cross-section dimensions)
- the fiber domain `[s_start, s_end]` — a `SubpathBuilt`'s domain starts at
  0 and runs to `arc_length(path)`

Operating wavelength `λ_m` is NOT stored on `Fiber`; it is an argument to
`generator_K` / `generator_Kω` (and to `propagate_fiber` in `path-integral.jl`),
so the same `Fiber` can be queried at multiple wavelengths. Temperature is
fixed at `T_ref_K` for all queries.

# ----------------------------
# Example Use
# ----------------------------

xs = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

sb = SubpathBuilder(); start!(sb)
straight!(sb; length = 5.0)
bend!(sb;
    radius = 4.458, angle = π / 2, axis_angle = 0.0,
    meta = [
        Nickname("90° bend"),
        MCMadd(:T_K, Normal(0.0, 2.0)),   # +ΔT_K ~ N(0, 2 K) on this segment
    ],
)
straight!(sb; length = 8.0)
# Seal the Subpath at the natural exit point; tests/demos commonly use a
# helper to compute the natural exit.
jumpto!(sb; point = (..., ..., ...))

# Pass the builder directly: Fiber builds the geometry and applies meta
# (here the bend's MCMadd(:T_K, …) thermal annotation).
fiber = Fiber(sb; cross_section = xs, T_ref_K = 297.15)

# Operating wavelength is supplied per query; temperature is f.T_ref_K.
K  = generator_K(fiber, 1550e-9)
Kω = generator_Kω(fiber, 1550e-9)
"""

if !isdefined(Main, :DEFAULT_T_REF_K)
    const DEFAULT_T_REF_K = 297.15
end

# Path-backed fibers use the path's local normal/binormal frame. The bend
# axis is the curvature normal, so the local transverse bend components are
# (κ, 0) in that frame. Frame rotation enters through the path spinning rate.
function bend_components(path::Union{SubpathBuilt, PathBuilt}, s::Real)
    κ = curvature(path, s)
    if κ == zero(κ)
        z = zero(κ)
        return (kx = z, ky = z, k2 = z)
    end
    z = zero(κ)
    return (kx = κ, ky = z, k2 = κ * κ)
end

struct Fiber{P,T,S,L}
    path::P
    cross_section::FiberCrossSection
    T_ref_K::T
    s_start::S
    s_end::S
    # Per-segment local temperature derived from `:T_K` (fiber-only): `nothing`
    # when no segment carries `:T_K` (then `local_temperature ≡ T_ref_K`), else a
    # `(breaks, vals)` pair giving the temperature on each placed segment in
    # global arc-length order. See `local_temperature`.
    local_T::L
end

function Fiber(
    path::Union{SubpathBuilt, PathBuilt};
    cross_section::FiberCrossSection,
    T_ref_K = DEFAULT_T_REF_K,
    local_T = nothing,
)
    s_start_val = 0.0
    s_end_val   = Float64(_qc_nominalize(arc_length(path)))
    s_start, s_end = promote(s_start_val, s_end_val)
    @assert s_end > s_start "Fiber requires s_end > s_start"
    return Fiber{typeof(path),typeof(T_ref_K),typeof(s_start),typeof(local_T)}(
        path,
        cross_section,
        T_ref_K,
        s_start,
        s_end,
        local_T,
    )
end

# ----------------------------
# Thermal (`:T_K`) interpretation — fiber-only
# ----------------------------
#
# `:T_K` is a foreign meta to the geometry layer (it cannot be resolved without a
# material). The fiber is its sole interpreter: it converts `:T_K` (a temperature
# excursion ΔT) into an isotropic length scaling `τ = 1 + α_lin·ΔT` using
# `α_lin = cte(cladding_material, T_ref_K)`, bakes that into the affected
# segments, strips `:T_K`, and lets the geometry `build(...; perturb=true)` apply
# any remaining field-level MCM. This is the *only* place `:T_K` is named.

# ΔT for a segment from its additive `:T_K` meta, or `nothing` when it carries
# none (the additive combine of 0.0 returns the unchanged baseline).
_segment_delta_T(seg) = (Δ = MCMcombine(0.0, seg, :T_K); Δ === 0.0 ? nothing : Δ)

# ΔT for the terminal connector from the seal's `:T_K` meta, or `nothing` when
# the seal carries none.
_seal_delta_T(sub::Subpath) =
    (Δ = MCMcombine(0.0, sub.jumpto_meta, :T_K); Δ === 0.0 ? nothing : Δ)

# The terminal connector supports only `MCMadd(:T_K, …)` (thermal expansion,
# issue #33). Any other MCMadd/MCMmul — field-level perturbation, or a
# multiplicative `:T_K` — has no effect on a solved connector, so reject it loudly
# rather than silently ignore. (This lives in the fiber: distinguishing the
# supported `:T_K` from field MCM requires naming `:T_K`, which the geometry layer
# must not do.)
function _validate_seal_meta(sub::Subpath)
    for m in sub.jumpto_meta
        bad = (m isa MCMmul) || (m isa MCMadd && m.symbol !== :T_K)
        bad && throw(ArgumentError(
            "jumpto!: the terminal connector supports only MCMadd(:T_K, …) thermal " *
            "meta (issue #33); got $(nameof(typeof(m)))(:$(m.symbol)). Field-level " *
            "MCMadd/MCMmul is not applied to a solved connector."))
    end
    return nothing
end

# Resolve a Subpath's `:T_K` meta into geometry: scale each thermal interior
# segment's length-fields by τ and strip its `:T_K`. If the terminal `jumpto!`
# connector carries `:T_K`, also compute its thermal target arc length (issue
# #33): the nominal connector length L0 scaled by τ_seal, re-solved to the fixed
# endpoint by `build`. Returns the resolved Subpath, that target length (or
# `nothing`), and the per-placed-segment ΔT excursions in placed order (interior
# segments then terminal connector) — `nothing` when nothing is thermal — so the
# fiber can record each segment's local temperature `T_ref_K + ΔT`.
function _resolve_thermal_subpath(sub::Subpath, cross_section::FiberCrossSection, T_ref_K)
    _validate_seal_meta(sub)   # reject unsupported MCM on the terminal connector
    seal_ΔT     = _seal_delta_T(sub)
    interior_TK = any(seg -> _segment_delta_T(seg) !== nothing, sub.segments)

    # Skip cte (and any thermal work) entirely when nothing is thermal, so a
    # non-thermal fiber on a cladding with no defined CTE still builds.
    (interior_TK || seal_ΔT !== nothing) || return (sub, nothing, nothing)

    # Per-placed-segment ΔT, aligned with the built placed segments (interior in
    # authored order, then the terminal connector).
    deltaT = Any[_segment_delta_T(seg) for seg in sub.segments]
    push!(deltaT, seal_ΔT)

    α_lin = cte(cross_section.cladding_material, T_ref_K)
    new_segments = AbstractPathSegment[
        let ΔT = _segment_delta_T(seg)
            ΔT === nothing ? seg :
                _scale_length_fields(seg, 1 + α_lin * ΔT, _meta_without(seg, :T_K))
        end
        for seg in sub.segments
    ]

    # Issue #33: terminal connector thermal target = τ_seal · L0, where L0 is the
    # nominal connector length (solved without :T_K). `build` re-solves to the
    # fixed `jumpto_point` with this arc length.
    jumpto_target_length = nothing
    if seal_ΔT !== nothing
        L0 = Float64(_qc_nominalize(
            arc_length(build(sub; perturb = false).jumpto_quintic_connector)))
        jumpto_target_length = (1 + α_lin * seal_ΔT) * L0
    end

    resolved = Subpath(
        sub.meta, sub.start_point, sub.start_outgoing_tangent,
        sub.start_outgoing_curvature, new_segments, sub.jumpto_point,
        sub.jumpto_incoming_tangent, sub.jumpto_incoming_curvature,
        sub.jumpto_min_bend_radius, _meta_without(sub.jumpto_meta, :T_K),
        sub.jumpto_twist,
        sub.jumpto_natural, sub.jumpto_natural_extra,
        sub.spin_rate, sub._spin_phi_at_s0,
    )
    return (resolved, jumpto_target_length, deltaT)
end

# Assemble the per-segment local-temperature table from a list of
# `(SubpathBuilt, global_offset)` pairs and the matching per-Subpath ΔT vectors
# (each `nothing` or aligned with that Subpath's placed segments). Returns
# `nothing` when no Subpath is thermal (so `local_temperature ≡ T_ref_K`), else a
# `(breaks, vals)` pair: `breaks[i]` is the global arc length at the end of
# placed segment `i`, `vals[i]` its local temperature. `vals` keeps each entry
# as-stored (no coercion) so a `Particles` ΔT propagates.
function _build_local_T(pairs, deltaTs, T_ref_K)
    any(!isnothing, deltaTs) || return nothing
    breaks = Float64[]
    vals = Any[]
    for k in eachindex(pairs)
        b, off = pairs[k]
        dT = deltaTs[k]
        placed = _all_placed_segs(b)
        for j in eachindex(placed)
            ps = placed[j]
            seg_end = off + Float64(_qc_nominalize(ps.s_offset_eff)) +
                            Float64(_qc_nominalize(arc_length(ps.segment)))
            δ = dT === nothing ? nothing : dT[j]
            push!(breaks, seg_end)
            push!(vals, δ === nothing ? T_ref_K : T_ref_K + δ)
        end
    end
    return (breaks, vals)
end

function _build_perturbed(sub::Subpath, cross_section::FiberCrossSection, T_ref_K)
    resolved, target, deltaT = _resolve_thermal_subpath(sub, cross_section, T_ref_K)
    built = build(resolved; perturb = true, jumpto_target_length = target)
    local_T = _build_local_T(((built, 0.0),), (deltaT,), T_ref_K)
    return built, local_T
end

function _build_perturbed(subs::Vector{Subpath}, cross_section::FiberCrossSection, T_ref_K)
    isempty(subs) && throw(ArgumentError("Fiber: at least one Subpath required"))
    # Build in order so `spin_rate = :inherit` resolves against the prior
    # thermal+perturbed built Subpath before this one is built.
    builts = Vector{SubpathBuilt}(undef, length(subs))
    deltaTs = Vector{Any}(undef, length(subs))
    for i in eachindex(subs)
        sub = i == 1 ? subs[i] :
              PathGeometry._resolve_inherited_spin(subs[i], builts[i-1])
        resolved, target, dT = _resolve_thermal_subpath(sub, cross_section, T_ref_K)
        builts[i] = build(resolved; perturb = true, jumpto_target_length = target)
        deltaTs[i] = dT
    end
    path = build(builts)
    offs = s_offsets(path)
    pairs = Tuple{SubpathBuilt,Float64}[(path.subpaths[i], offs[i])
                                        for i in eachindex(path.subpaths)]
    local_T = _build_local_T(pairs, deltaTs, T_ref_K)
    return path, local_T
end

"""
    Fiber(spec; cross_section, T_ref_K=DEFAULT_T_REF_K) -> Fiber

Build a fiber from authored geometry, applying perturbation meta during the build.
`spec` may be a `SubpathBuilder`, a `Subpath`, a `Vector{Subpath}`, or a
`Vector{SubpathBuilder}` (each builder is frozen to a `Subpath` first, so
thermal handling is identical to the `Vector{Subpath}` path).

Thermal `:T_K` annotations are resolved here using
`α_lin = cte(cross_section.cladding_material, T_ref_K)` — each thermal segment's
length-dimensioned fields are scaled by `1 + α_lin·ΔT`. Field-level `MCMadd`/`MCMmul`
are then applied by the geometry build. If the terminal `jumpto!` connector carries
`:T_K`, it thermally expands too (issue #33): its arc length scales by `τ_seal` while
still landing at the fixed `jumpto_point`. The geometry is built exactly once.
"""
Fiber(spec::SubpathBuilder; cross_section::FiberCrossSection, T_ref_K = DEFAULT_T_REF_K) =
    Fiber(Subpath(spec); cross_section = cross_section, T_ref_K = T_ref_K)

function Fiber(spec::Subpath; cross_section::FiberCrossSection, T_ref_K = DEFAULT_T_REF_K)
    built, local_T = _build_perturbed(spec, cross_section, T_ref_K)
    return Fiber(built; cross_section = cross_section, T_ref_K = T_ref_K, local_T = local_T)
end

function Fiber(spec::Vector{Subpath}; cross_section::FiberCrossSection,
               T_ref_K = DEFAULT_T_REF_K)
    built, local_T = _build_perturbed(spec, cross_section, T_ref_K)
    return Fiber(built; cross_section = cross_section, T_ref_K = T_ref_K, local_T = local_T)
end

Fiber(spec::Vector{SubpathBuilder}; cross_section::FiberCrossSection,
      T_ref_K = DEFAULT_T_REF_K) =
    Fiber(Subpath[Subpath(b) for b in spec];
          cross_section = cross_section, T_ref_K = T_ref_K)

fiber_path(f::Fiber) = f.path

"""
    local_temperature(f::Fiber, s) -> T

Local temperature (K) at fiber arc length `s`: `T_ref_K + ΔT` for the placed
segment containing `s`, where `ΔT` is that segment's `:T_K` excursion. A fiber
with no `:T_K` returns `T_ref_K` everywhere. The cross-section birefringences
are evaluated at this temperature (asymmetric thermal stress ∝ `|T_soft − T|`,
and the indices shift with `T`), so a `:T_K` segment's optical response — not
only its length — reflects the excursion. `ΔT` may be `Particles`; the returned
value carries it into the MCM-safe cross-section `T_K` slot.
"""
function local_temperature(f::Fiber, s::Real)
    lt = f.local_T
    lt === nothing && return f.T_ref_K
    breaks, vals = lt
    @inbounds for i in eachindex(breaks)
        s <= breaks[i] + 1e-9 && return vals[i]
    end
    return vals[end]
end

# ----------------------------
# Generator K(s) and Curvature Kω(s)
# ----------------------------

zero_generator() = zeros(ComplexF64, 2, 2)

"""
    linear_birefringence_generator(Δβ, c2φ, s2φ) -> 2×2 matrix

Local Jones generator for a linear retarder with retardance per unit length
`Δβ` and eigen-axes oriented at angle `φ` in the propagation frame, encoded
via `c2φ = cos(2φ)` and `s2φ = sin(2φ)`. Traceless, anti-Hermitian times `i`.
Shared by the bend, core-ellipticity, and asymmetric-thermal-stress generators.
"""
linear_birefringence_generator(Δβ, c2φ, s2φ) = [
     0.5im * Δβ * c2φ    0.5im * Δβ * s2φ
     0.5im * Δβ * s2φ   -0.5im * Δβ * c2φ
]

"""
    circular_birefringence_generator(Δβc) -> 2×2 matrix

Local Jones generator for circular birefringence (optical activity) with
rotation rate `Δβc`. Real antisymmetric — a pure SO(2) rotation generator.
Used by the mechanical-twist generator.
"""
circular_birefringence_generator(Δβc) = [
     zero(Δβc)   -0.5 * Δβc
     0.5 * Δβc    zero(Δβc)
]

# Orientation of the bend (curvature-direction) birefringence axis in the
# parallel-transport (Bishop) propagation frame. The Frenet curvature direction
# rotates relative to the Bishop frame at the geometric torsion rate, so its
# angle is the accumulated `∫₀ˢ τ_g` (zero for any planar path, recovering the
# previous fixed `(κ,0)` behaviour). Returned as `(cos 2φ, sin 2φ)`.
function _bend_axis_c2s2(f::Fiber, s::Real)
    φ = torsion_phase(f.path, s)
    return (cos(2 * φ), sin(2 * φ))
end

function bend_generator_K(f::Fiber, s::Real, λ_m::Real)
    κ = curvature(f.path, s)
    κ == zero(κ) && return zero_generator()
    T = local_temperature(f, s)
    R = inv(κ)
    Δβb = bending_birefringence(f.cross_section, λ_m, T; bend_radius_m = R)
    c2φ, s2φ = _bend_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβb, c2φ, s2φ)
end

function bend_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    κ = curvature(f.path, s)
    κ == zero(κ) && return zero_generator()
    T = local_temperature(f, s)
    R = inv(κ)
    Δβbω = bending_birefringence(
        WithDerivative(),
        f.cross_section,
        λ_m,
        T;
        bend_radius_m = R
    ).dω
    c2φ, s2φ = _bend_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβbω, c2φ, s2φ)
end

# Mechanical twist (τ_m) is the *only* source of circular birefringence
# (photoelastic optical activity). Geometric torsion and material spin rotate
# linear axes instead and are handled by the bend and ellipticity generators.
function twist_generator_K(f::Fiber, s::Real, λ_m::Real)
    τm = twist_rate(f.path, s)
    τm == zero(τm) && return zero_generator()
    T = local_temperature(f, s)
    Δβc = twisting_birefringence(f.cross_section, λ_m, T; twist_rate_rad_per_m = τm)
    return circular_birefringence_generator(Δβc)
end

function twist_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    τm = twist_rate(f.path, s)
    τm == zero(τm) && return zero_generator()
    T = local_temperature(f, s)
    Δβcω = twisting_birefringence(
        WithDerivative(),
        f.cross_section,
        λ_m,
        T;
        twist_rate_rad_per_m = τm
    ).dω
    return circular_birefringence_generator(Δβcω)
end

# Orientation of the intrinsic-linear (core ellipticity + asymmetric thermal
# stress) birefringence axes in the Bishop frame: the frozen ellipse angle plus
# the material rotation from spin and mechanical twist. `(cos 2φ, sin 2φ)`.
function _intrinsic_axis_c2s2(f::Fiber, s::Real)
    φ = f.cross_section.ellipticity_axis_angle +
        spin_phase(f.path, s) + twist_phase(f.path, s)
    return (cos(2 * φ), sin(2 * φ))
end

# Core ellipticity and asymmetric thermal stress share the ellipse eigen-axes,
# so their magnitudes add before the (single) linear generator. A circular core
# (axis ratio 1) contributes nothing — guarded so a circular fiber pays no cost
# and is bit-for-bit unchanged from the pre-ellipticity model.
function ellipticity_generator_K(f::Fiber, s::Real, λ_m::Real)
    xs = f.cross_section
    xs.ellipticity_axis_ratio == one(xs.ellipticity_axis_ratio) && return zero_generator()
    T = local_temperature(f, s)
    Δβ = core_noncircularity_birefringence(xs, λ_m, T) +
         asymmetric_thermal_stress_birefringence(xs, λ_m, T)
    c2φ, s2φ = _intrinsic_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβ, c2φ, s2φ)
end

function ellipticity_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    xs = f.cross_section
    xs.ellipticity_axis_ratio == one(xs.ellipticity_axis_ratio) && return zero_generator()
    T = local_temperature(f, s)
    Δβω = core_noncircularity_birefringence(WithDerivative(), xs, λ_m, T).dω +
          asymmetric_thermal_stress_birefringence(WithDerivative(), xs, λ_m, T).dω
    c2φ, s2φ = _intrinsic_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβω, c2φ, s2φ)
end

fiber_breakpoints(f::Fiber) = breakpoints(f.path)

# ----------------------------
# Generator dispatch on cross-section type
# ----------------------------
#
# The cross-section layer selects the local physics: `generator_K(f, xs, λ_m)`
# dispatches on `typeof(xs)`. The 2-argument forms are convenience wrappers that
# dispatch on `f.cross_section`, so callers may write `generator_K(fiber, λ_m)`.

# Generic fallbacks — a concrete cross-section type must provide the assembly.
generator_K(f::Fiber, xs::FiberCrossSection, λ_m::Real) =
    error("generator_K not implemented for $(typeof(f)) and $(typeof(xs))")
generator_Kω(f::Fiber, xs::FiberCrossSection, λ_m::Real) =
    error("generator_Kω not implemented for $(typeof(f)) and $(typeof(xs))")

"""
    generator_K(fiber, λ_m) -> (s -> 2×2 ComplexF64)

Return a closure that evaluates the local Jones generator `K(s)` at the given
operating wavelength `λ_m` (metres), dispatching on the fiber's cross section.
Temperature is `fiber.T_ref_K`.
"""
generator_K(f::Fiber, λ_m::Real) = generator_K(f, f.cross_section, λ_m)

"""
    generator_Kω(fiber, λ_m) -> (s -> 2×2 ComplexF64)

Frequency-derivative counterpart of `generator_K`.
"""
generator_Kω(f::Fiber, λ_m::Real) = generator_Kω(f, f.cross_section, λ_m)

# ----------------------------
# Step-index fiber generators
# ----------------------------

function generator_K(f::Fiber, xs::StepIndexCrossSection, λ_m::Real)
    return function (s::Real)
        return bend_generator_K(f, s, λ_m) +        # linear: bending / tension
               twist_generator_K(f, s, λ_m) +       # circular: mechanical twist
               ellipticity_generator_K(f, s, λ_m)   # linear: ellipticity + stress
    end
end

function generator_Kω(f::Fiber, xs::StepIndexCrossSection, λ_m::Real)
    return function (s::Real)
        return bend_generator_Kω(f, s, λ_m) +
               twist_generator_Kω(f, s, λ_m) +
               ellipticity_generator_Kω(f, s, λ_m)
    end
end

# ----------------------------
# Graded-index fiber generators (not yet modeled)
# ----------------------------

generator_K(f::Fiber, xs::GradedIndexCrossSection, λ_m::Real) =
    error("generator_K not implemented for $(typeof(f)) and $(typeof(xs))")
generator_Kω(f::Fiber, xs::GradedIndexCrossSection, λ_m::Real) =
    error("generator_Kω not implemented for $(typeof(f)) and $(typeof(xs))")

# ----------------------------
# Fiber diagnostics for plotting
# ----------------------------

function bend_geometry(f::Fiber, s::Real)
    curv = bend_components(f.path, s)
    kx = curv.kx
    ky = curv.ky
    k2 = kx * kx + ky * ky
    if k2 == 0.0
        return (Rb = Inf, theta_b = 0.0, kx = 0.0, ky = 0.0, k2 = 0.0)
    end

    return (Rb = inv(sqrt(k2)), theta_b = atan(ky, kx), kx = kx, ky = ky, k2 = k2)
end
