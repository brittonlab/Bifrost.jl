# FiberCrossSection (from FiberCS) and the Subpath geometry API (from
# PathGeometry) are brought into scope by the FiberPath submodule in Bifrost.jl.

# Default reference temperature (K), ≈ 24 °C.
if !isdefined(Main, :DEFAULT_T_REF_K)
    const DEFAULT_T_REF_K = 297.15
end

"""
    bend_components(path, s) -> (kx, ky, k2)

Return the local transverse bend components of `path` at arc length `s`.

Path-backed fibers use the path's local normal/binormal frame. The bend axis is
the curvature normal, so the components are `(κ, 0)` in that frame with
`k2 = κ²`; frame rotation enters through the path spin rate.
"""
function bend_components(path::Union{SubpathBuilt, PathBuilt}, s::Real)
    κ = curvature(path, s)
    if κ == zero(κ)
        z = zero(κ)
        return (kx = z, ky = z, k2 = z)
    end
    z = zero(κ)
    return (kx = κ, ky = z, k2 = κ * κ)
end

"""
    Fiber

Compiled, path-backed fiber: immutable built geometry bound to a cross section.

A `Fiber` is the query object consumed by the propagation layer. It owns:

- the immutable built path (`SubpathBuilt` or `PathBuilt`),
- the `FiberCrossSection`,
- the reference temperature `T_ref_K` (K), and
- the fiber domain `[s_start, s_end]` in global arc length (a built path's
  domain starts at `0` and runs to `arc_length(path)`).

Operating wavelength `λ_m` is **not** stored: it is a per-query argument to
[`generator_K`](@ref) / [`generator_Kω`](@ref) (and to `propagate_fiber`), so
one `Fiber` can be queried at several wavelengths. `T_ref_K` is the baseline
temperature; per-segment `:T_K` meta overrides it locally (see
[`temperature`](@ref)).
"""
struct Fiber{P,T,S}
    path::P
    cross_section::FiberCrossSection
    T_ref_K::T
    s_start::S
    s_end::S
end

function Fiber(
    path::Union{SubpathBuilt, PathBuilt};
    cross_section::FiberCrossSection,
    T_ref_K = DEFAULT_T_REF_K,
)
    s_start_val = 0.0
    s_end_val   = Float64(_qc_nominalize(arc_length(path)))
    s_start, s_end = promote(s_start_val, s_end_val)
    @assert s_end > s_start "Fiber requires s_end > s_start"
    return Fiber{typeof(path),typeof(T_ref_K),typeof(s_start)}(
        path,
        cross_section,
        T_ref_K,
        s_start,
        s_end,
    )
end

# ----------------------------
# Thermal (`:T_K`) and tension (`:tension`) interpretation — fiber-only
# ----------------------------
# Both are foreign meta to the geometry layer (they cannot be resolved without a
# material); the fiber is their sole interpreter, in
# `_resolve_thermal_and_tension` below. This is the only place they are named.

"""
    _segment_delta_T(seg)

Return the temperature excursion ΔT (K) from a segment's additive `:T_K` meta,
or `nothing` when it carries none (the additive combine of `0.0` returns the
unchanged baseline).
"""
_segment_delta_T(seg) = (Δ = MCMcombine(0.0, seg, :T_K); Δ === 0.0 ? nothing : Δ)

"""
    _seal_delta_T(sub)

Return ΔT (K) for the terminal connector from the seal's `:T_K` meta, or
`nothing` when the seal carries none.
"""
_seal_delta_T(sub::Subpath) =
    (Δ = MCMcombine(0.0, sub.jumpto_meta, :T_K); Δ === 0.0 ? nothing : Δ)

"""
    _segment_tension(seg)

Return the absolute axial tension (N) from a segment's additive `:tension`
meta, or `nothing` when it carries none.
"""
_segment_tension(seg) = (F = MCMcombine(0.0, seg, :tension); F === 0.0 ? nothing : F)

"""
    _seal_tension(sub)

Return the axial tension (N) for the terminal connector from the seal's
`:tension` meta, or `nothing` when the seal carries none.
"""
_seal_tension(sub::Subpath) =
    (F = MCMcombine(0.0, sub.jumpto_meta, :tension); F === 0.0 ? nothing : F)

"""
    _tension_strain(F, r_clad, E)

Return the axial strain `ε = F / (π·r_clad²·E)` that elongates a segment by
`(1 + ε)` under an absolute tension `F` (N). Shares the denominator of
`axial_tension_dω`.
"""
_tension_strain(F, r_clad, E) = F / (π * r_clad^2 * E)

"""
    _validate_seal_meta(sub)

Throw an `ArgumentError` unless the terminal connector carries only
`MCMadd(:T_K, …)` (thermal expansion) or `MCMadd(:tension, …)` (axial
elongation) meta.

Any other `MCMadd`/`MCMmul` — field-level perturbation, or a multiplicative
`:T_K`/`:tension` — has no effect on a solved connector, so it is rejected
loudly rather than silently ignored. This lives in the fiber because
distinguishing the supported additive meta from field MCM requires naming
them, which the geometry layer must not do.
"""
function _validate_seal_meta(sub::Subpath)
    for m in sub.jumpto_meta
        bad = (m isa MCMmul) ||
              (m isa MCMadd && m.symbol !== :T_K && m.symbol !== :tension)
        bad && throw(ArgumentError(
            "jumpto!: the terminal connector supports only MCMadd(:T_K, …) thermal " *
            "or MCMadd(:tension, …) elongation meta; got " *
            "$(nameof(typeof(m)))(:$(m.symbol)). Field-level MCMadd/MCMmul is not " *
            "applied to a solved connector."))
    end
    return nothing
end

"""
    _resolve_thermal_and_tension(sub, cross_section, T_ref_K)
        -> (resolved_subpath, jumpto_target_length)

Resolve a Subpath's `:T_K` and `:tension` meta into geometry.

`:T_K` is a temperature excursion ΔT (K) from `T_ref_K`, converted into the
length scaling `1 + α_lin·ΔT` with `α_lin = cte(cladding_material, T_ref_K)`.
`:tension` is an absolute axial tension `F` (N), converted into the elongation
`1 + ε` with the strain of [`_tension_strain`](@ref), using
`E = youngs_modulus(cladding_material, T_ref_K)` and
`r_clad = cladding_radius(cross_section)` — the cladding convention, for
consistency between the two paths.

Each affected interior segment's length-fields are scaled by the *combined*
factor `τ = (1 + α_lin·ΔT)·(1 + ε)` (the two scalings compose
multiplicatively) while *retaining* the meta: the segment keeps
`:T_K`/`:tension` as records, so [`temperature`](@ref) and [`tension`](@ref)
recover them on demand at query time and no thermal or tension state is stored
on `Fiber`. (Retaining `:T_K` also lets the optical response — not only the
length — reflect the excursion.)

If the terminal `jumpto!` connector carries either meta, the returned target
arc length is the nominal connector length `L0` scaled by the same combined
factor; `build` re-solves the connector to the fixed endpoint with that length.
Otherwise the target is `nothing`. The material `cte`/`youngs_modulus` lookups
are gated behind the presence of the corresponding meta, so a plain fiber on a
cladding with neither property defined still builds.
"""
function _resolve_thermal_and_tension(sub::Subpath, cross_section::FiberCrossSection, T_ref_K)
    _validate_seal_meta(sub)   # reject unsupported MCM on the terminal connector
    seal_ΔT     = _seal_delta_T(sub)
    seal_F      = _seal_tension(sub)
    interior_TK = any(seg -> _segment_delta_T(seg) !== nothing, sub.segments)
    interior_F  = any(seg -> _segment_tension(seg) !== nothing, sub.segments)

    # Skip all material lookups (and any meta work) when nothing is thermal or tensioned. 
    any_thermal = interior_TK || seal_ΔT !== nothing
    any_tension = interior_F || seal_F !== nothing
    (any_thermal || any_tension) || return (sub, nothing)

    # Consult cladding CTE/stiffness lazily — only the one(s) actually needed.
    α_lin  = any_thermal ? cte(cross_section.cladding_material, T_ref_K) : nothing
    E_clad = any_tension ? youngs_modulus(cross_section.cladding_material, T_ref_K) :
             nothing
    r_clad = cladding_radius(cross_section)

    # Combined length factor for one segment from its `:T_K` (ΔT) and `:tension`
    # (F) meta; either may be absent (`nothing` → unit factor).
    _τ(ΔT, F) =
        (ΔT === nothing ? 1 : (1 + α_lin * ΔT)) *
        (F === nothing ? 1 : (1 + _tension_strain(F, r_clad, E_clad)))

    # Scale one segment's length-fields by its combined τ; a segment with neither
    # `:T_K` nor `:tension` passes through unchanged.
    function _scaled(seg)
        ΔT = _segment_delta_T(seg)
        F  = _segment_tension(seg)
        if ΔT === nothing && F === nothing
            return seg
        end
        return _scale_length_fields(seg, _τ(ΔT, F), segment_meta(seg))
    end

    new_segments = AbstractPathSegment[_scaled(seg) for seg in sub.segments]

    # for jumpto the terminal connector target = τ_seal · L0
    # Here,  L0 is the nominal connector length (solved without :T_K/:tension). 
    # `build` re-solves to the fixed `jumpto_point` with this arc length.
    # When the seal expands, the terminal connector elongates by `seal_factor` (its
    # length is re-solved to `jumpto_target_length`), so its twist rate divides by
    # the same factor to conserve connector turns — mirroring the interior segments.
    seal_expands = seal_ΔT !== nothing || seal_F !== nothing
    seal_factor = seal_expands ? _τ(seal_ΔT, seal_F) : 1
    jumpto_target_length = nothing
    if seal_expands
        L0 = Float64(_qc_nominalize(
            arc_length(build(sub; perturb = false).jumpto_quintic_connector)))
        jumpto_target_length = seal_factor * L0
    end

    resolved = Subpath(
        sub.meta, sub.start_point, sub.start_outgoing_tangent,
        sub.start_outgoing_curvature, new_segments, sub.jumpto_point,
        sub.jumpto_incoming_tangent, sub.jumpto_incoming_curvature,
        sub.jumpto_min_bend_radius, sub.jumpto_meta,
        _scale_inverse_twist_rate(sub.jumpto_twist, seal_factor),
        sub.jumpto_natural, sub.jumpto_natural_extra,
        sub.spin_rate, sub._spin_phi_at_s0,
        sub.inherit_start_point, sub.inherit_start_tangent, sub.inherit_start_curvature,
    )
    return (resolved, jumpto_target_length)
end

"""
    _build_perturbed(sub_or_subs, cross_section, T_ref_K)

Build authored geometry once: resolve `:T_K`/`:tension` meta via
[`_resolve_thermal_and_tension`](@ref), then `build(...; perturb = true)` so the
geometry layer applies any remaining field-level MCM. For a `Vector{Subpath}`,
Subpaths are built in order so `:inherit` start state and spin resolve against
the prior thermal+perturbed built Subpath.
"""
function _build_perturbed(sub::Subpath, cross_section::FiberCrossSection, T_ref_K)
    resolved, target = _resolve_thermal_and_tension(sub, cross_section, T_ref_K)
    return build(resolved; perturb = true, jumpto_target_length = target)
end

function _build_perturbed(subs::Vector{Subpath}, cross_section::FiberCrossSection, T_ref_K)
    isempty(subs) && throw(ArgumentError("Fiber: at least one Subpath required"))
    # Build in order so `spin_rate = :inherit` resolves against the prior
    # thermal+perturbed built Subpath before this one is built.
    builts = Vector{SubpathBuilt}(undef, length(subs))
    for i in eachindex(subs)
        sub = subs[i]
        if i > 1
            # Resolve start-state then spin inheritance against the prior
            # thermal+perturbed built Subpath, mirroring build(::Vector{Subpath}).
            sub = PathGeometry._resolve_inherited_start(sub, builts[i-1])
            sub = PathGeometry._resolve_inherited_spin(sub, builts[i-1])
        end
        resolved, target = _resolve_thermal_and_tension(sub, cross_section, T_ref_K)
        builts[i] = build(resolved; perturb = true, jumpto_target_length = target)
    end
    return build(builts)
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
`:T_K`, it thermally expands too: its arc length scales by `τ_seal` while
still landing at the fixed `jumpto_point`. The geometry is built exactly once.

Per-segment axial tension `:tension` (absolute Newtons, e.g.
`MCMadd(:tension, 0.5)`) is interpreted here too and plays a dual
role exactly like `:T_K`: it elongates the segment by `(1 + ε)` with the axial
strain `ε = F / (π·r_clad²·E)` (`E = youngs_modulus(cladding_material, T_ref_K)`,
`r_clad = cladding_radius(cross_section)`) *and* it sets the segment's axial-tension
photoelastic birefringence (a linear birefringence on the bend eigen-axis, so it
is zero on a straight segment; recovered downstream via `tension(f, s)`). Thermal
and tension length scalings compose multiplicatively. A terminal `jumpto!` carrying
`:tension` elongates the connector the same way. The `:tension` meta is left on the
segment for on-demand recovery, as `:T_K` is.

# Examples
```julia
xs = StepIndexCrossSection(
    SilicaGermaniaGlass(0.036),
    SilicaGermaniaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

sb = SubpathBuilder()
start!(sb)
straight!(sb; length = 5.0)
bend!(sb;
    radius = 4.458, angle = π / 2, axis_angle = 0.0,
    meta = [
        Nickname("90° bend"),
        MCMadd(:T_K, Normal(0.0, 2.0)),   # +ΔT_K ~ N(0, 2 K) on this segment
    ],
)
straight!(sb; length = 8.0)
seal!(sb)   # end at the natural exit

# Pass the builder directly: Fiber builds the geometry and applies meta
# (here the bend's MCMadd(:T_K, …) thermal annotation).
fiber = Fiber(sb; cross_section = xs, T_ref_K = 297.15)

# Operating wavelength is supplied per query; temperature comes from
# temperature(fiber, s).
K  = generator_K(fiber, 1550e-9)
Kω = generator_Kω(fiber, 1550e-9)
```
"""
Fiber(spec::SubpathBuilder; cross_section::FiberCrossSection, T_ref_K = DEFAULT_T_REF_K) =
    Fiber(Subpath(spec); cross_section = cross_section, T_ref_K = T_ref_K)

function Fiber(spec::Subpath; cross_section::FiberCrossSection, T_ref_K = DEFAULT_T_REF_K)
    built = _build_perturbed(spec, cross_section, T_ref_K)
    return Fiber(built; cross_section = cross_section, T_ref_K = T_ref_K)
end

function Fiber(spec::Vector{Subpath}; cross_section::FiberCrossSection,
               T_ref_K = DEFAULT_T_REF_K)
    built = _build_perturbed(spec, cross_section, T_ref_K)
    return Fiber(built; cross_section = cross_section, T_ref_K = T_ref_K)
end

Fiber(spec::Vector{SubpathBuilder}; cross_section::FiberCrossSection,
      T_ref_K = DEFAULT_T_REF_K) =
    Fiber(Subpath[Subpath(b) for b in spec];
          cross_section = cross_section, T_ref_K = T_ref_K)

"""
    fiber_path(f::Fiber)

Return the built path (`SubpathBuilt` or `PathBuilt`) backing `f`.
"""
fiber_path(f::Fiber) = f.path

"""
    temperature(f::Fiber, s) -> T

Return the temperature (K) at fiber arc length `s`, derived on demand as
`T_ref_K + ΔT`, where `ΔT` is the `:T_K` excursion carried by the segment
containing `s`.

A fiber with no `:T_K` returns `T_ref_K` everywhere — `MCMcombine(0.0, seg,
:T_K)` collapses to `0.0` when the segment carries none, so no sentinel is
needed, and no thermal state is stored on `Fiber`. The cross-section
birefringences are evaluated at this temperature (asymmetric thermal stress
∝ `|T_soft − T|`, and the indices shift with `T`), so a `:T_K` segment's
optical response — not only its length — reflects the excursion. `ΔT` may be
`Particles`; the returned value carries it into the MCM-safe cross-section
`T_K` slot.
"""
function temperature(f::Fiber, s::Real)
    return f.T_ref_K + MCMcombine(0.0, local_segment(f.path, s), :T_K)
end

"""
    tension(f::Fiber, s) -> F

Return the axial tension (N) at fiber arc length `s`, derived on demand from
the `:tension` meta carried by the segment containing `s`; `0` where none is
present.

`MCMcombine(0.0, seg, :tension)` collapses to `0.0` when the segment carries no
`:tension`, so no sentinel is needed, and no tension state is stored on
`Fiber`. Consumed by [`tension_generator_K`](@ref) (axial-tension photoelastic
birefringence). May be `Particles`.
"""
tension(f::Fiber, s::Real) = MCMcombine(0.0, local_segment(f.path, s), :tension)

# ----------------------------
# Generator K(s) and Curvature Kω(s)
# ----------------------------

"""
    zero_generator()

Return a `2×2` zero `ComplexF64` matrix, the null local Jones generator
contribution.
"""
zero_generator() = zeros(ComplexF64, 2, 2)

"""
    linear_birefringence_generator(Δβ, c2φ, s2φ) -> 2×2 matrix

Return the local Jones generator for a linear retarder with retardance per unit length
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

Return the local Jones generator for circular birefringence (optical activity) with
rotation rate `Δβc`. Real antisymmetric — a pure SO(2) rotation generator.
Used by the mechanical-twist generator.
"""
circular_birefringence_generator(Δβc) = [
     zero(Δβc)   -0.5 * Δβc
     0.5 * Δβc    zero(Δβc)
]

"""
    _bend_axis_c2s2(f, s) -> (c2φ, s2φ)

Return `(cos 2φ, sin 2φ)` for the orientation `φ` of the bend
(curvature-direction) birefringence axis in the parallel-transport (Bishop)
propagation frame.

The Frenet curvature direction rotates relative to the Bishop frame at the
geometric torsion rate, so `φ` is the accumulated `∫₀ˢ τ_g` (zero for any
planar path, recovering fixed `(κ, 0)` behaviour).
"""
function _bend_axis_c2s2(f::Fiber, s::Real)
    φ = torsion_phase(f.path, s)
    return (cos(2 * φ), sin(2 * φ))
end

"""
    bend_generator_K(f::Fiber, s, λ_m)

Return the bend contribution to the local Jones generator `K(s)`: a linear
birefringence of magnitude `bending_birefringence` on the curvature eigen-axis
(oriented by [`_bend_axis_c2s2`](@ref)). Zero on a straight section.
"""
function bend_generator_K(f::Fiber, s::Real, λ_m::Real)
    κ = curvature(f.path, s)
    κ == zero(κ) && return zero_generator()
    T = temperature(f, s)
    R = inv(κ)
    Δβb = bending_birefringence(f.cross_section, λ_m, T; bend_radius_m = R)
    c2φ, s2φ = _bend_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβb, c2φ, s2φ)
end

"""
    bend_generator_Kω(f::Fiber, s, λ_m)

Frequency-derivative counterpart of [`bend_generator_K`](@ref).
"""
function bend_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    κ = curvature(f.path, s)
    κ == zero(κ) && return zero_generator()
    T = temperature(f, s)
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

"""
    twist_generator_K(f::Fiber, s, λ_m)

Return the mechanical-twist contribution to the local Jones generator `K(s)`: a
circular birefringence of magnitude `twisting_birefringence`.

Mechanical twist `τ_m` is the *only* source of circular birefringence
(photoelastic optical activity); geometric torsion and material spin rotate
linear axes instead and are handled by the bend and ellipticity generators.
"""
function twist_generator_K(f::Fiber, s::Real, λ_m::Real)
    τm = twist_rate(f.path, s)
    τm == zero(τm) && return zero_generator()
    T = temperature(f, s)
    Δβc = twisting_birefringence(f.cross_section, λ_m, T; twist_rate_rad_per_m = τm)
    return circular_birefringence_generator(Δβc)
end

"""
    twist_generator_Kω(f::Fiber, s, λ_m)

Frequency-derivative counterpart of [`twist_generator_K`](@ref).
"""
function twist_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    τm = twist_rate(f.path, s)
    τm == zero(τm) && return zero_generator()
    T = temperature(f, s)
    Δβcω = twisting_birefringence(
        WithDerivative(),
        f.cross_section,
        λ_m,
        T;
        twist_rate_rad_per_m = τm
    ).dω
    return circular_birefringence_generator(Δβcω)
end

"""
    _intrinsic_axis_c2s2(f, s) -> (c2φ, s2φ)

Return `(cos 2φ, sin 2φ)` for the orientation `φ` of the intrinsic-linear (core
ellipticity + asymmetric thermal stress) birefringence axes in the Bishop
frame: the frozen ellipse angle plus the material rotation from spin and
mechanical twist.
"""
function _intrinsic_axis_c2s2(f::Fiber, s::Real)
    φ = f.cross_section.ellipticity_axis_angle +
        spin_phase(f.path, s) + twist_phase(f.path, s)
    return (cos(2 * φ), sin(2 * φ))
end

"""
    ellipticity_generator_K(f::Fiber, s, λ_m)

Return the intrinsic core-ellipticity contribution to the local Jones generator
`K(s)`.

Core ellipticity and asymmetric thermal stress share the ellipse eigen-axes
(oriented by [`_intrinsic_axis_c2s2`](@ref)), so their magnitudes add before
the single linear generator. A circular core (axis ratio `1`) contributes
nothing — guarded so a circular fiber pays no cost.
"""
function ellipticity_generator_K(f::Fiber, s::Real, λ_m::Real)
    xs = f.cross_section
    xs.ellipticity_axis_ratio == one(xs.ellipticity_axis_ratio) && return zero_generator()
    T = temperature(f, s)
    Δβ = core_noncircularity_birefringence(xs, λ_m, T) +
         asymmetric_thermal_stress_birefringence(xs, λ_m, T)
    c2φ, s2φ = _intrinsic_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβ, c2φ, s2φ)
end

"""
    ellipticity_generator_Kω(f::Fiber, s, λ_m)

Frequency-derivative counterpart of [`ellipticity_generator_K`](@ref).
"""
function ellipticity_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    xs = f.cross_section
    xs.ellipticity_axis_ratio == one(xs.ellipticity_axis_ratio) && return zero_generator()
    T = temperature(f, s)
    Δβω = core_noncircularity_birefringence(WithDerivative(), xs, λ_m, T).dω +
          asymmetric_thermal_stress_birefringence(WithDerivative(), xs, λ_m, T).dω
    c2φ, s2φ = _intrinsic_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβω, c2φ, s2φ)
end

"""
    tension_generator_K(f::Fiber, s, λ_m)

Return the axial-tension photoelastic contribution to the local Jones generator
`K(s)`, an additive term separate from bending.

Like bending, it is a *linear* birefringence on the bend eigen-axis and
vanishes on a straight segment (the response ∝ `1/R`), so a tensioned but
unbent segment contributes nothing. Temperature enters through the
photoelastic/stiffness constants exactly as for the bend term.
"""
function tension_generator_K(f::Fiber, s::Real, λ_m::Real)
    F = tension(f, s)
    F == zero(F) && return zero_generator()
    κ = curvature(f.path, s)
    κ == zero(κ) && return zero_generator()   # tension birefringence ∝ 1/R
    T = temperature(f, s)
    Δβ = axial_tension_birefringence(f.cross_section, λ_m, T;
                                     bend_radius_m = inv(κ), axial_tension_N = F)
    c2φ, s2φ = _bend_axis_c2s2(f, s)           # shares the bend eigen-axis
    return linear_birefringence_generator(Δβ, c2φ, s2φ)
end

"""
    tension_generator_Kω(f::Fiber, s, λ_m)

Frequency-derivative counterpart of [`tension_generator_K`](@ref).
"""
function tension_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    F = tension(f, s)
    F == zero(F) && return zero_generator()
    κ = curvature(f.path, s)
    κ == zero(κ) && return zero_generator()
    T = temperature(f, s)
    Δβω = axial_tension_birefringence(
        WithDerivative(),
        f.cross_section,
        λ_m,
        T;
        bend_radius_m   = inv(κ),
        axial_tension_N = F,
    ).dω
    c2φ, s2φ = _bend_axis_c2s2(f, s)
    return linear_birefringence_generator(Δβω, c2φ, s2φ)
end

"""
    fiber_breakpoints(f::Fiber)

Return the globally merged path breakpoints of `f`; the propagator never steps
across one.
"""
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
The local temperature comes from [`temperature`](@ref)`(fiber, s)`.
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
        return bend_generator_K(f, s, λ_m) +        # linear: bending
               twist_generator_K(f, s, λ_m) +       # circular: mechanical twist
               ellipticity_generator_K(f, s, λ_m) + # linear: ellipticity + stress
               tension_generator_K(f, s, λ_m)       # linear: axial tension
    end
end

function generator_Kω(f::Fiber, xs::StepIndexCrossSection, λ_m::Real)
    return function (s::Real)
        return bend_generator_Kω(f, s, λ_m) +
               twist_generator_Kω(f, s, λ_m) +
               ellipticity_generator_Kω(f, s, λ_m) +
               tension_generator_Kω(f, s, λ_m)
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

"""
    bend_geometry(f::Fiber, s) -> (Rb, theta_b, kx, ky, k2)

Return local bend diagnostics for plotting: bend radius `Rb` (m, `Inf` when
straight), bend-axis angle `theta_b` (rad), and the transverse curvature
components from [`bend_components`](@ref).
"""
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
