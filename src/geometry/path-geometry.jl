"""
Construct and query three-dimensional smooth space curves (fiber paths).

This file defines a small authoring DSL for building piecewise paths from
local-frame segment primitives, compiling them to an immutable global-frame
form, and querying differential geometry (position, tangent/normal/binormal,
curvature, torsion) and material spinning along arc length.

# Three layers

- `SubpathBuilder` (mutable) — authoring target for the bang-DSL below. A
  Subpath specification begins with `start!(...)` and ends with a seal:
  `jumpto!(...)` (bend toward a global target) or `seal!(...)` (end at the
  natural exit, no bending).
- `Subpath` (immutable) — frozen snapshot of user-supplied data only.
- `SubpathBuilt` (immutable) — derived layout: `subpath` + `placed_segments`
  + `jumpto_quintic_connector` + `resolved_spinning`.
- `PathBuilt` (immutable) — ordered container of `SubpathBuilt`s.

`build(builder_or_subpath) → SubpathBuilt` runs the placement loop on one
Subpath. `build(::Vector{Subpath})` or `build(::Vector{SubpathBuilt})`
produces a `PathBuilt`.

# Authoring

Sliding-frame interior segments + a global terminal jumpto:

    sb = SubpathBuilder()
    start!(sb; point=(0,0,0), outgoing_tangent=(0,0,1))
    straight!(sb; length=1.0)
    bend!(sb; radius=0.05, angle=π/2)
    jumpto!(sb; point=(0.05, 0.0, 1.05))     # seals the Subpath
    sub = Subpath(sb)
    sb_built = build(sub)

To end a Subpath exactly as authored (no terminal connector bending), seal with
`seal!` instead of `jumpto!`:

    seal!(sb)             # end at the natural exit
    seal!(sb; extra=0.02) # ...plus a 2 cm straight lead-out

Interior `jumpby!` is supported as a relative jump within a Subpath. Each
Subpath must call `start!` before any segment and seal exactly once at the end
(`jumpto!` or `seal!`).

# Spinning

Spinning is attached as per-segment meta via `Spinning <: AbstractMeta`.
A `Spinning` placed in a segment's `meta` vector starts a spinning run at that
segment's start and continues until the next `Spinning`-bearing segment, or
until the Subpath ends. The terminal `jumpto!`/`seal!` connector carries `meta`
too and is treated as the Subpath's last segment, so a `Spinning` on the seal
starts a run over the connector just like one on any interior segment.

`rate` may be a `Real` (constant rad/m) or a `Function` `rate(s_local)`. A
`Spinning` with `is_continuous = true` placed on the first spinning anchor of a
Subpath leaves the resolution pending; `build(::Vector{SubpathBuilt})` then
inherits `phi_0` from the prior Subpath's terminal phase.

# Independence

Each `Subpath` and `SubpathBuilt` is fully independent of all others. The
Subpath holds its `start_point` and `jumpto_point` as the only globally-anchored
values. `build(::Vector{SubpathBuilt})` is the first ordering-aware layer; it
checks endpoint conformity between adjacent Subpaths and resolves any
cross-Subpath spinning continuity.

# Interface

    arc_length(seg_or_path)
    arc_length(path, s1, s2)
    curvature(seg_or_path, s)
    geometric_torsion(seg_or_path, s)
    spinning_rate(path, s)
    position(path, s)
    tangent(path, s)
    normal(path, s)
    binormal(path, s)
    frame(path, s)
    breakpoints(path)
    sample(path, s_values)
    sample_uniform(path; n)
"""

using LinearAlgebra
using QuadGK

# -----------------------------------------------------------------------
# Abstract segment type
# -----------------------------------------------------------------------

"""
    AbstractPathSegment

Supertype for all local-frame path segment primitives (straight, bend,
catenary, helix, and the resolved quintic connector).

Each concrete segment defines its geometry in its own local frame, starting at
the origin with tangent along local `+z`. The `build` placement loop rotates
and chains these local descriptions into the global frame.

# Implementation

A concrete segment must implement the following, where local arc length is
`s ∈ [0, arc_length(seg)]`. The return element type `T` follows the segment's
own type parameter (`Float64` for deterministic segments, `Particles` for
MCM-valued fields):

| Method | Returns |
| --- | --- |
| `arc_length(seg)` | `T` |
| `curvature(seg, s)` | `T` (κ, 1/m) |
| `geometric_torsion(seg, s)` | `T` (τ_geom, rad/m) |
| `position_local(seg, s)` | length-3 `Vector{T}` |
| `tangent_local(seg, s)` | unit length-3 `Vector{T}` |
| `normal_local(seg, s)` | unit length-3 `Vector{T}` |
| `binormal_local(seg, s)` | unit length-3 `Vector{T}` |
| `end_position_local(seg)` | length-3 `Vector{T}` |
| `end_frame_local(seg)` | `(T, N, B)`, each length-3 `Vector{T}` |

`JumpBy` resolves at `build` time into a parametric [`QuinticConnector`](@ref)
that supports MCM `Particles` via nominalized-scalar branching.
"""
abstract type AbstractPathSegment end

# -----------------------------------------------------------------------
# AbstractMeta — per-segment annotations
# -----------------------------------------------------------------------

"""
    AbstractMeta

Supertype for per-segment annotations carried in a segment's
`meta::Vector{AbstractMeta}` bag.

`path-geometry.jl` is deliberately ignorant of what the bag contains;
downstream layers define their own `AbstractMeta` subtypes (see
`path-geometry-meta.jl`) and decide how to act on them. The only `AbstractMeta`
this file interprets is [`Spinning`](@ref).
"""
abstract type AbstractMeta end

"""
    segment_meta(seg::AbstractPathSegment) -> Vector{AbstractMeta}

Return the segment's `meta` annotation bag, or an empty vector if it has none.
"""
segment_meta(seg::AbstractPathSegment) =
    :meta ∈ fieldnames(typeof(seg)) ? seg.meta : AbstractMeta[]

# -----------------------------------------------------------------------
# Spinning  (per-segment material-spinning annotation)
# -----------------------------------------------------------------------

"""
    Spinning(; rate, phi_0 = 0.0, is_continuous = false)

Per-segment material-spinning annotation with a constant or function-valued
rate.

Attached to a segment's `meta` vector. The spinning run begins at the host
segment's start (effective arc length `s_offset_eff`) and extends until the
next segment carrying a `Spinning`, or until the path ends.

- `rate` may be a `Real` (constant rad/m) or a `Function` `rate(s_local)` of
  run-local arc length where `s_local = 0` at the start of the spinning run.
- `phi_0` is the absolute initial phase (rad) at the start of the run when
  `is_continuous = false`.
- `is_continuous = true` means the resolver computes `phi_0` from the prior
  spinning run's accumulated phase. In that case `phi_0` must be left at its
  default `0.0`.

At most one `Spinning` may appear in any single segment's `meta`. The first
`Spinning` on a path must have `is_continuous = false`.
"""
struct Spinning <: AbstractMeta
    rate::Union{Float64, Function}
    phi_0::Float64
    is_continuous::Bool

    function Spinning(; rate, phi_0 = 0.0, is_continuous::Bool = false)
        if is_continuous && phi_0 != 0.0
            throw(ArgumentError(
                "Spinning: do not specify phi_0 when is_continuous=true (phase is carried over from prior run)"))
        end
        r = rate isa Function ? rate : Float64(rate)
        new(r, Float64(phi_0), is_continuous)
    end
end

"""
    ResolvedSpinningRate(s_eff_start, s_eff_end, rate, phi_0)

A single spinning run resolved to absolute path coordinates. `rate` is called
(when a `Function`) with run-local arc length `s_local = s - s_eff_start`.
`phi_0` is the absolute phase (rad) at `s_eff_start`.
"""
struct ResolvedSpinningRate
    s_eff_start::Float64
    s_eff_end::Float64
    rate::Union{Float64, Function}
    phi_0::Float64
end

# -----------------------------------------------------------------------
# Quadrature helper
# -----------------------------------------------------------------------

"""
    _integrate_rate(rate, a, b; rtol=1e-8, atol=0.0)

Integrate a spinning `rate` over the run-local interval `[a, b]`.

Constant rates take the analytic branch; `Function` rates use QuadGK adaptive
Gauss–Kronrod, which subdivides automatically for oscillatory integrands.
"""
_integrate_rate(rate::Float64, a::Float64, b::Float64;
                rtol::Float64 = 1e-8, atol::Float64 = 0.0) = rate * (b - a)

function _integrate_rate(rate::Function, a::Float64, b::Float64;
                         rtol::Float64 = 1e-8, atol::Float64 = 0.0)
    val, _err = QuadGK.quadgk(rate, a, b; rtol = rtol, atol = atol)
    return val
end

# -----------------------------------------------------------------------
# StraightSegment
# -----------------------------------------------------------------------

"""
    StraightSegment{T} <: AbstractPathSegment

Straight line segment of signed `length` along the local tangent.

A negative `length` walks backward along the local tangent: `arc_length` is
`|length|`, position advances as `sign(length)·s`, and the end-frame
tangent/normal carry the same sign so downstream segments inherit a consistent,
right-handed frame matching the direction of motion. Using `sign` (not a
conditional) keeps this branch-free and MCM/`Particles`-compatible — `sign`
broadcasts elementwise over a `Particles` length.
"""
struct StraightSegment{T} <: AbstractPathSegment
    length::T
    meta::Vector{AbstractMeta}

    """
        StraightSegment(length; meta=AbstractMeta[]) -> StraightSegment

    Construct a straight segment of signed `length`. See
    [`AbstractPathSegment`](@ref).
    """
    function StraightSegment(length; meta = AbstractMeta[])
        new{typeof(length)}(length, Vector{AbstractMeta}(meta))
    end
end

"""
    arc_length(seg::AbstractPathSegment) -> T

Return the segment's total arc length (the local-frame interface method; see
[`AbstractPathSegment`](@ref)).
"""
arc_length(seg::StraightSegment)         = abs(seg.length)

"""
    curvature(seg::AbstractPathSegment, s) -> T

Return the curvature κ (1/m) at local arc length `s`. See
[`AbstractPathSegment`](@ref).
"""
curvature(seg::StraightSegment, _)       = zero(seg.length)

"""
    geometric_torsion(seg::AbstractPathSegment, s) -> T

Return the geometric torsion τ (rad/m) at local arc length `s`. See
[`AbstractPathSegment`](@ref).
"""
geometric_torsion(seg::StraightSegment, _) = zero(seg.length)

"""
    position_local(seg::AbstractPathSegment, s) -> Vector{T}

Return the position at local arc length `s`, in the segment's local frame. See
[`AbstractPathSegment`](@ref).
"""
position_local(seg::StraightSegment, s)   = [zero(s), zero(s), sign(seg.length) * s]

"""
    tangent_local(seg::AbstractPathSegment, s) -> Vector{T}

Return the unit tangent at local arc length `s`, in the local frame. See
[`AbstractPathSegment`](@ref).
"""
tangent_local(seg::StraightSegment, _)    = [zero(seg.length), zero(seg.length), sign(seg.length)]

"""
    normal_local(seg::AbstractPathSegment, s) -> Vector{T}

Return the unit normal at local arc length `s`, in the local frame. See
[`AbstractPathSegment`](@ref).
"""
normal_local(seg::StraightSegment, _)     = [sign(seg.length), zero(seg.length), zero(seg.length)]

"""
    binormal_local(seg::AbstractPathSegment, s) -> Vector{T}

Return the unit binormal at local arc length `s`, in the local frame. See
[`AbstractPathSegment`](@ref).
"""
binormal_local(seg::StraightSegment, _)   = [zero(seg.length), one(seg.length), zero(seg.length)]

"""
    end_position_local(seg::AbstractPathSegment) -> Vector{T}

Return the segment's endpoint position in its local frame. See
[`AbstractPathSegment`](@ref).
"""
end_position_local(seg::StraightSegment)  = [zero(seg.length), zero(seg.length), seg.length]

"""
    end_frame_local(seg::AbstractPathSegment) -> (T, N, B)

Return the `(tangent, normal, binormal)` frame at the segment's end, in its
local frame. See [`AbstractPathSegment`](@ref).
"""
function end_frame_local(seg::StraightSegment)
    sgn = sign(seg.length)
    T = [zero(seg.length), zero(seg.length), sgn]
    N = [sgn, zero(seg.length), zero(seg.length)]
    B = [zero(seg.length), one(seg.length), zero(seg.length)]
    return (T, N, B)
end

# -----------------------------------------------------------------------
# BendSegment  (circular arc)
# -----------------------------------------------------------------------

"""
    BendSegment(radius, angle, axis_angle)

Circular arc of radius `radius` (m) sweeping `angle` (rad) in the plane whose
inward normal is at `axis_angle` (rad) from the local N-axis.

In the local frame (local z = incoming tangent, local x = incoming normal,
local y = incoming binormal), the inward normal direction is:
    n̂ = cos(axis_angle)·x̂ + sin(axis_angle)·ŷ

Curvature κ = 1 / radius.
"""
struct BendSegment{T} <: AbstractPathSegment
    radius::T
    angle::T       # total angle swept (rad)
    axis_angle::T  # orientation of inward normal in transverse plane (rad)
    meta::Vector{AbstractMeta}

    function BendSegment(radius, angle, axis_angle = 0.0;
                         meta = AbstractMeta[])
        @assert radius > 0 "BendSegment: radius must be positive"
        r, a, x = promote(radius, angle, axis_angle)
        new{typeof(r)}(r, a, x, Vector{AbstractMeta}(meta))
    end
end

arc_length(seg::BendSegment)         = seg.radius * abs(seg.angle)
curvature(seg::BendSegment, _)       = one(seg.radius) / seg.radius
geometric_torsion(seg::BendSegment, _) = zero(seg.radius)

function position_local(seg::BendSegment, s)
    R   = seg.radius
    θ   = s / R
    φ   = seg.axis_angle
    n̂   = [cos(φ), sin(φ), zero(φ)]
    return R * (1 - cos(θ)) * n̂ + [zero(R), zero(R), R * sin(θ)]
end

function tangent_local(seg::BendSegment, s)
    R = seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]
end

function normal_local(seg::BendSegment, s)
    R = seg.radius
    θ = s / R
    φ = seg.axis_angle
    return [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
end

function binormal_local(seg::BendSegment, _)
    φ = seg.axis_angle
    return [-sin(φ), cos(φ), zero(φ)]   # constant for circular arc (zero torsion)
end

function end_position_local(seg::BendSegment)
    R = seg.radius
    θ = seg.angle
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), zero(φ)]
    return R * (1 - cos(θ)) * n̂ + [zero(R), zero(R), R * sin(θ)]
end

function end_frame_local(seg::BendSegment)
    θ = seg.angle
    φ = seg.axis_angle
    T = [sin(θ) * cos(φ), sin(θ) * sin(φ),  cos(θ)]
    N = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]
    B = [-sin(φ), cos(φ), zero(φ)]
    return (T, N, B)
end

# -----------------------------------------------------------------------
# CatenarySegment
# -----------------------------------------------------------------------

"""
    CatenarySegment(a, length, axis_angle)

A catenary curve in the plane whose horizontal direction is at `axis_angle`
from the local N-axis.  `a` (m) is the catenary parameter (a = T₀/(ρg),
the ratio of horizontal tension to weight per unit length).  The curve
starts with tangent along the local z-axis (vertical catenary vertex) and
curves horizontally.  Curvature κ(s) = a / (a² + s²).
"""
struct CatenarySegment{T} <: AbstractPathSegment
    a::T
    length::T
    axis_angle::T
    meta::Vector{AbstractMeta}

    function CatenarySegment(a, length, axis_angle = 0.0;
                             meta = AbstractMeta[])
        @assert a > 0      "CatenarySegment: a must be positive"
        @assert length > 0 "CatenarySegment: length must be positive"
        av, L, x = promote(a, length, axis_angle)
        new{typeof(av)}(av, L, x, Vector{AbstractMeta}(meta))
    end
end

arc_length(seg::CatenarySegment)         = seg.length
geometric_torsion(seg::CatenarySegment, _) = zero(seg.a)

function curvature(seg::CatenarySegment, s)
    a = seg.a
    return a / (a^2 + s^2)
end

function position_local(seg::CatenarySegment, s)
    a = seg.a
    φ = seg.axis_angle
    n̂ = [cos(φ), sin(φ), zero(φ)]
    horiz = a * (sqrt(1 + (s / a)^2) - 1)
    vert  = a * asinh(s / a)
    return horiz * n̂ + [zero(a), zero(a), vert]
end

function tangent_local(seg::CatenarySegment, s)
    a = seg.a
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [(s / a) / q * cos(φ), (s / a) / q * sin(φ), one(q) / q]
end

function normal_local(seg::CatenarySegment, s)
    # N = dT/ds / |dT/ds|, derived analytically: N = [n̂_horiz/q, -s/a/q] normalised
    a = seg.a
    φ = seg.axis_angle
    q = sqrt(1 + (s / a)^2)
    return [cos(φ) / q, sin(φ) / q, -(s / a) / q]
end

function binormal_local(seg::CatenarySegment, s)
    return cross(tangent_local(seg, s), normal_local(seg, s))
end

function end_position_local(seg::CatenarySegment)
    return position_local(seg, arc_length(seg))
end

function end_frame_local(seg::CatenarySegment)
    s = arc_length(seg)
    return (tangent_local(seg, s), normal_local(seg, s), binormal_local(seg, s))
end

# -----------------------------------------------------------------------
# HelixSegment
# -----------------------------------------------------------------------

"""
    HelixSegment(radius, pitch, turns, axis_angle)

A helix whose entry tangent is ẑ (the incoming sliding-frame tangent), ensuring
continuity with the prior segment.  `axis_angle` (rad) selects which transverse
direction n̂ = cos(axis_angle)·x̂ + sin(axis_angle)·ŷ the helix curves toward.

The helix axis â = (h·ẑ + R·n̂) / ℓ' is tilted from the transverse plane by
arctan(h/R) toward ẑ, where h = pitch/(2π) and ℓ' = √(R²+h²).  This tilt is
the geometric consequence of demanding tangent(0) = ẑ: a zero-pitch helix
reduces to a circular arc (BendSegment) in the n̂ direction.

    κ(s) = R / (R² + h²)          (constant)
    τ_geom(s) = h / (R² + h²)     (constant)
    arc_length = turns · 2π · √(R² + h²)

Local-frame basis vectors:
    n̂  = [cos(axis_angle), sin(axis_angle), 0]     (toward helix axis)
    r̂₀ = [-sin(axis_angle), cos(axis_angle), 0]    (outward radial at s=0)
    ê_φ = (R·ẑ - h·n̂) / ℓ'                        (tangential at s=0, ⊥ axis)
"""
struct HelixSegment{T} <: AbstractPathSegment
    radius::T
    pitch::T
    turns::T
    axis_angle::T
    meta::Vector{AbstractMeta}

    function HelixSegment(radius, pitch, turns,
                          axis_angle = 0.0; meta = AbstractMeta[])
        @assert radius > 0 "HelixSegment: radius must be positive"
        @assert turns  > 0 "HelixSegment: turns must be positive"
        r, p, n, x = promote(radius, pitch, turns, axis_angle)
        new{typeof(r)}(r, p, n, x, Vector{AbstractMeta}(meta))
    end
end

"""
    _helix_h(seg::HelixSegment)

Return the helix axial advance per radian, `h = pitch / 2π`.
"""
function _helix_h(seg::HelixSegment)
    seg.pitch / (2π)          # axial advance per radian
end

function arc_length(seg::HelixSegment)
    h = _helix_h(seg)
    return seg.turns * 2π * sqrt(seg.radius^2 + h^2)
end

function curvature(seg::HelixSegment, _)
    R = seg.radius
    h = _helix_h(seg)
    return R / (R^2 + h^2)
end

function geometric_torsion(seg::HelixSegment, _)
    R = seg.radius
    h = _helix_h(seg)
    return h / (R^2 + h^2)
end

"""
    _helix_basis(seg::HelixSegment)

Return the helix geometric basis `(R, h, ℓ′, â, r̂₀, ê_φ)`: radius, axial
advance per radian, arc length per radian, axis unit vector, outward radial at
`s = 0`, and tangential direction at `s = 0`.
"""
function _helix_basis(seg::HelixSegment)
    R  = seg.radius
    h  = _helix_h(seg)
    φ_a = seg.axis_angle
    ℓ′  = sqrt(R^2 + h^2)
    z0 = zero(R); z1 = one(R)
    ẑ   = [z0, z0, z1]
    n̂   = [cos(φ_a), sin(φ_a), z0]
    â   = (h .* ẑ .+ R .* n̂) ./ ℓ′
    r̂₀  = [-sin(φ_a), cos(φ_a), z0]   # ẑ × n̂
    ê_φ = (R .* ẑ .- h .* n̂) ./ ℓ′  # â × r̂₀, tangential at s=0
    return R, h, ℓ′, â, r̂₀, ê_φ
end

function position_local(seg::HelixSegment, s)
    R, h, ℓ′, â, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    return h .* â .* φ .+ R .* (cos(φ) - 1) .* r̂₀ .+ R .* sin(φ) .* ê_φ
end

function tangent_local(seg::HelixSegment, s)
    R, h, ℓ′, â, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    return (h .* â .- R .* sin(φ) .* r̂₀ .+ R .* cos(φ) .* ê_φ) ./ ℓ′
end

function normal_local(seg::HelixSegment, s)
    _, _, ℓ′, _, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    return -(cos(φ) .* r̂₀ .+ sin(φ) .* ê_φ)
end

function binormal_local(seg::HelixSegment, s)
    R, h, ℓ′, â, r̂₀, ê_φ = _helix_basis(seg)
    φ = s / ℓ′
    # B = T × N; expanding in {â, r̂₀, ê_φ} orthonormal frame:
    # B = (R·â + h·sin φ·r̂₀ - h·cos φ·ê_φ) / ℓ′
    return (R .* â .+ h .* sin(φ) .* r̂₀ .- h .* cos(φ) .* ê_φ) ./ ℓ′
end

function end_position_local(seg::HelixSegment)
    return position_local(seg, arc_length(seg))
end

function end_frame_local(seg::HelixSegment)
    s = arc_length(seg)
    return (tangent_local(seg, s), normal_local(seg, s), binormal_local(seg, s))
end

"""
    JumpBy <: AbstractPathSegment

Relative jump connector resolved at build time into a quintic G2 Hermite curve.

Connect the current position to `current_position + delta`. The incoming
tangent and curvature vector are inherited
from the prior segment automatically. `tangent_out` is the desired outgoing
tangent direction (relative-frame); if nothing, falls back to the chord
direction. `curvature_out` is the desired outgoing curvature vector dT/ds in
the relative frame; defaults to zero (G1 outgoing match).

`min_bend_radius` (metres, nothing = unconstrained) sets a lower bound on the
radius of curvature of the connector. The handle scale λ is extended beyond
the chord default (with bisection refinement) to keep κ ≤ 1/R_min.
"""
struct JumpBy <: AbstractPathSegment
    delta::NTuple{3, Float64}
    tangent_out::Union{Nothing, NTuple{3, Float64}}
    curvature_out::Union{Nothing, NTuple{3, Float64}}
    min_bend_radius::Union{Nothing, Float64}
    meta::Vector{AbstractMeta}

    """
        JumpBy(delta; tangent_out=nothing, curvature_out=nothing,
                      min_bend_radius=nothing, meta=AbstractMeta[]) -> JumpBy

    Construct a relative jump connector to `current_position + delta`, resolved
    at `build` time. See the type docstring above for the full argument
    semantics.
    """
    function JumpBy(delta; tangent_out = nothing, curvature_out = nothing,
                    min_bend_radius = nothing, meta = AbstractMeta[])
        d = (Float64(delta[1]), Float64(delta[2]), Float64(delta[3]))
        t = isnothing(tangent_out) ? nothing :
            (Float64(tangent_out[1]), Float64(tangent_out[2]), Float64(tangent_out[3]))
        k = isnothing(curvature_out) ? nothing :
            (Float64(curvature_out[1]), Float64(curvature_out[2]), Float64(curvature_out[3]))
        r = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
        new(d, t, k, r, Vector{AbstractMeta}(meta))
    end
end

# JumpBy is context-dependent: geometry is resolved at build() time into a
# QuinticConnector. Calling these methods on the raw struct is unsupported.
arc_length(::JumpBy)             = error("JumpBy: call build() to resolve jump geometry")
curvature(::JumpBy, ::Real)      = error("JumpBy: call build() to resolve jump geometry")
geometric_torsion(::JumpBy, ::Real) = 0.0
position_local(::JumpBy, ::Real) = error("JumpBy: call build() to resolve jump geometry")
tangent_local(::JumpBy, ::Real)  = error("JumpBy: call build() to resolve jump geometry")
normal_local(::JumpBy, ::Real)   = error("JumpBy: call build() to resolve jump geometry")
binormal_local(::JumpBy, ::Real) = error("JumpBy: call build() to resolve jump geometry")
end_position_local(::JumpBy)     = error("JumpBy: call build() to resolve jump geometry")
end_frame_local(::JumpBy)        = error("JumpBy: call build() to resolve jump geometry")

# -----------------------------------------------------------------------
# QuinticConnector  (resolved form of JumpBy and Subpath terminal connectors)
# -----------------------------------------------------------------------

include(joinpath(@__DIR__, "path-geometry-connector.jl"))

# Concrete meta vocabulary lives in path-geometry-meta.jl. It is part of the
# geometry layer (Spinning, Nickname, MCMadd, MCMmul) and makes no reference to
# fiber. Included here so Subpath constructors can reference MCMadd/MCMmul for
# validation.
include(joinpath(@__DIR__, "path-geometry-meta.jl"))

"""
    _resolve_at_placement(seg, pos, frame_mat, K_in_global)

Resolve a `JumpBy` to a `QuinticConnector` at `build` time; pass other segment
types through unchanged.

`K_in_global` is the prior segment's terminal curvature vector in the global
frame, computed by `build`. It is rotated into the new segment's local frame
here, mirroring the tangent-rotation pattern.
"""
function _resolve_at_placement(seg::JumpBy, pos::AbstractVector, frame_mat::AbstractMatrix,
                                K_in_global::AbstractVector)
    p1_local  = collect(seg.delta)
    chord     = norm(p1_local)
    t_hat_out = isnothing(seg.tangent_out) ?
        (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
        normalize(collect(seg.tangent_out))
    K0_local = frame_mat' * K_in_global
    K1_local = isnothing(seg.curvature_out) ? zeros(eltype(K0_local), 3) :
                                              collect(seg.curvature_out)
    return _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                                    min_bend_radius = seg.min_bend_radius,
                                    meta = seg.meta)
end

# Default fallback for non-JumpBy segments (ignores K_in_global).
_resolve_at_placement(seg::AbstractPathSegment, ::AbstractVector, ::AbstractMatrix,
                      ::AbstractVector) = seg

# -----------------------------------------------------------------------
# PlacedSegment  (derived layout — lives inside PathSpecCached)
# -----------------------------------------------------------------------

"""
    PlacedSegment

A segment together with its placement in the global frame.

Records the segment, its cumulative arc-length offset at its start
(`s_offset_eff`), its global start `origin`, and the 3×3 `frame` whose columns
are `[N_global | B_global | T_global]` and which transforms local vectors to
global as `v_g = frame * v_l`.
"""
struct PlacedSegment
    segment::AbstractPathSegment
    s_offset_eff::Real          # cumulative arc-length at segment start
    origin::AbstractVector      # global start position (length 3)
    frame::AbstractMatrix       # 3×3, columns [N_global | B_global | T_global]
                                # transforms local vectors → global: v_g = frame * v_l
end

# -----------------------------------------------------------------------
# SubpathBuilder (mutable authoring) → Subpath (immutable spec)
# -----------------------------------------------------------------------

"""
    SubpathBuilder()

Mutable Subpath authoring target. Lifecycle:

1. `start!(builder; ...)` seals the start state. Must be called before any
   interior segment.
2. Append interior segments via `straight!`, `bend!`, `helix!`, `catenary!`,
   `jumpby!`.
3. Seal the end with `jumpto!(builder; ...)` (bend toward a global target) or
   `seal!(builder; ...)` (end at the natural exit, no bending). After this no
   further segments may be added.

`Subpath(builder)` (or `build(builder)`) requires both the start and the end
seal to have been called.
"""
mutable struct SubpathBuilder
    meta::Vector{AbstractMeta}
    start_point::Union{Nothing, NTuple{3, Float64}}
    start_outgoing_tangent::Union{Nothing, NTuple{3, Float64}}
    start_outgoing_curvature::Union{Nothing, NTuple{3, Float64}}
    segments::Vector{AbstractPathSegment}
    jumpto_point::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    # Meta attached to the terminal `jumpto!` connector. The geometry layer does
    # not interpret it; a consuming layer (the fiber) may — e.g. to thermally
    # expand the connector via build's `jumpto_target_length`.
    jumpto_meta::Vector{AbstractMeta}
    # Natural seal (set by seal!): seal at the natural exit with no connector
    # bending. `jumpto_natural_extra` is an optional straight lead-out length.
    jumpto_natural::Bool
    jumpto_natural_extra::Float64
    # `:inherit` start-state intent (set by start!): when a flag is true the
    # corresponding start_* field is a placeholder, resolved at vector build time
    # from the predecessor Subpath's endpoint. The start_* tuples stay concrete
    # Float64 so all downstream readers remain type-stable.
    inherit_start_point::Bool
    inherit_start_tangent::Bool
    inherit_start_curvature::Bool

    SubpathBuilder(; meta::AbstractVector{<:AbstractMeta} = AbstractMeta[]) =
        new(Vector{AbstractMeta}(meta),
            nothing, nothing, nothing,
            AbstractPathSegment[],
            nothing, nothing, nothing, nothing, AbstractMeta[],
            false, 0.0,
            false, false, false)
end

"""
    _check_started(b::SubpathBuilder)

Throw if `start!` has not yet been called on the builder.
"""
_check_started(b::SubpathBuilder) =
    isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: call start!() before adding segments"))

"""
    _check_unsealed(b::SubpathBuilder)

Throw if the builder has already been sealed by `jumpto!` or `seal!`.
"""
_check_unsealed(b::SubpathBuilder) =
    (!isnothing(b.jumpto_point) || b.jumpto_natural) &&
        throw(ArgumentError("SubpathBuilder: subpath already sealed by jumpto!() or seal!()"))

"""
    start!(builder; point=(0,0,0), outgoing_tangent=(0,0,1),
                    outgoing_curvature=(0,0,0))
    start!(builder, :inherit)

Seal the Subpath start state. Throws if `start!` has already been called or
if any interior segment has already been appended.

Any of `point`, `outgoing_tangent`, `outgoing_curvature` may be the symbol
`:inherit`, meaning "take this field from the previous Subpath's endpoint". The
inheritance is resolved at vector build time (`build([...])`) from the prior
built Subpath, so `:inherit` is only valid for a non-first Subpath in a list;
`build` errors otherwise. The positional form `start!(builder, :inherit)` is
shorthand for inheriting all three fields. See [`build`](@ref).
"""
function start!(b::SubpathBuilder;
                point = (0.0, 0.0, 0.0),
                outgoing_tangent = (0.0, 0.0, 1.0),
                outgoing_curvature = (0.0, 0.0, 0.0))
    !isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: start!() already called"))
    !isempty(b.segments) &&
        throw(ArgumentError("SubpathBuilder: start!() must be called before any segments"))
    b.inherit_start_point      = _is_inherit(point,              "point")
    b.inherit_start_tangent    = _is_inherit(outgoing_tangent,   "outgoing_tangent")
    b.inherit_start_curvature  = _is_inherit(outgoing_curvature, "outgoing_curvature")
    # Inherited fields hold a concrete placeholder; the real value is filled in at
    # vector build time, so start_* stays type-stable NTuple{3,Float64}.
    b.start_point              = b.inherit_start_point     ? (0.0, 0.0, 0.0) :
        (Float64(point[1]),              Float64(point[2]),              Float64(point[3]))
    b.start_outgoing_tangent   = b.inherit_start_tangent   ? (0.0, 0.0, 1.0) :
        (Float64(outgoing_tangent[1]),   Float64(outgoing_tangent[2]),   Float64(outgoing_tangent[3]))
    b.start_outgoing_curvature = b.inherit_start_curvature ? (0.0, 0.0, 0.0) :
        (Float64(outgoing_curvature[1]), Float64(outgoing_curvature[2]), Float64(outgoing_curvature[3]))
    return b
end

# Positional shorthand: start!(builder, :inherit) inherits all three start
# fields from the previous Subpath. Only :inherit is accepted.
function start!(b::SubpathBuilder, mode::Symbol)
    mode === :inherit ||
        throw(ArgumentError("start!: positional argument must be :inherit; got :$mode"))
    return start!(b; point = :inherit, outgoing_tangent = :inherit,
                  outgoing_curvature = :inherit)
end

"""
    _is_inherit(value, field_name) -> Bool

Return `true` if `value` is the `:inherit` sentinel, `false` if it is an ordinary
tuple. Any other symbol is rejected with a pointed error naming `field_name`.
"""
_is_inherit(value, ::AbstractString) = false
function _is_inherit(value::Symbol, field_name::AbstractString)
    value === :inherit ||
        throw(ArgumentError("start!: $field_name accepts a tuple or :inherit; got :$value"))
    return true
end

"""
    straight!(builder; length, meta=AbstractMeta[]) -> builder

Append a [`StraightSegment`](@ref) of signed `length` to `builder`.

See also [`bend!`](@ref), [`helix!`](@ref), [`catenary!`](@ref),
[`jumpby!`](@ref), [`seal!`](@ref), [`jumpto!`](@ref).
"""
function straight!(spec::SubpathBuilder; length,
                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, StraightSegment(length; meta))
    return spec
end

"""
    bend!(builder; radius, angle, axis_angle=0.0, meta=AbstractMeta[]) -> builder

Append a [`BendSegment`](@ref) (circular arc) to `builder`. `axis_angle`
orients the bend plane in the transverse frame.

See also [`straight!`](@ref), [`helix!`](@ref), [`catenary!`](@ref),
[`jumpby!`](@ref).
"""
function bend!(spec::SubpathBuilder; radius::Real, angle::Real, axis_angle::Real = 0.0,
               meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, BendSegment(radius, angle, axis_angle; meta))
    return spec
end

"""
    helix!(builder; radius, pitch, turns, axis_angle=0.0, meta=AbstractMeta[]) -> builder

Append a [`HelixSegment`](@ref) to `builder`. `axis_angle` selects the
transverse direction the helix curves toward.

See also [`straight!`](@ref), [`bend!`](@ref), [`catenary!`](@ref),
[`jumpby!`](@ref).
"""
function helix!(spec::SubpathBuilder; radius::Real, pitch::Real, turns::Real,
                axis_angle::Real = 0.0, meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, HelixSegment(radius, pitch, turns, axis_angle; meta))
    return spec
end

"""
    catenary!(builder; a, length, axis_angle=0.0, meta=AbstractMeta[]) -> builder

Append a [`CatenarySegment`](@ref) to `builder`. `a` is the catenary parameter
and `axis_angle` orients the curve in the transverse frame.

See also [`straight!`](@ref), [`bend!`](@ref), [`helix!`](@ref),
[`jumpby!`](@ref).
"""
function catenary!(spec::SubpathBuilder; a::Real, length::Real, axis_angle::Real = 0.0,
                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, CatenarySegment(a, length, axis_angle; meta))
    return spec
end

"""
    jumpby!(builder; delta, tangent=nothing, curvature_out=nothing,
                     min_bend_radius=nothing, meta=AbstractMeta[]) -> builder

Append a relative [`JumpBy`](@ref) connector that advances the path by `delta`,
resolved into a quintic G2 connector at `build` time.

See also [`jumpto!`](@ref) (the terminal global-target form) and the segment
builders [`straight!`](@ref), [`bend!`](@ref), [`helix!`](@ref),
[`catenary!`](@ref).
"""
function jumpby!(spec::SubpathBuilder; delta, tangent = nothing,
                 curvature_out = nothing, min_bend_radius = nothing,
                 meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    _check_started(spec); _check_unsealed(spec)
    push!(spec.segments, JumpBy(delta; tangent_out = tangent,
                                curvature_out = curvature_out,
                                min_bend_radius, meta))
    return spec
end

"""
    jumpto!(builder; point, incoming_tangent=nothing, incoming_curvature=nothing,
                     min_bend_radius=nothing, meta=AbstractMeta[])

Seal the Subpath end by bending toward a global target `point`, storing the
terminal connector spec on the builder.

`meta` attaches annotations to the terminal connector, which is placed and
queried exactly like any interior segment: a `Spinning` here starts a spinning
run over the connector (see [`Spinning`](@ref)). Meta the geometry layer does
not recognize is carried through for a consuming layer (e.g. the fiber reads a
thermal annotation and expands the connector to land at `point` with a scaled
arc length — see `build`'s `jumpto_target_length`).

Throws if `start!` was not called, or if the builder is already sealed. See
also [`seal!`](@ref) to end at the natural exit without bending.
"""
function jumpto!(b::SubpathBuilder;
                 point,
                 incoming_tangent = nothing,
                 incoming_curvature = nothing,
                 min_bend_radius = nothing,
                 meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: call start!() before jumpto!()"))
    (!isnothing(b.jumpto_point) || b.jumpto_natural) &&
        throw(ArgumentError("SubpathBuilder: subpath already sealed"))
    b.jumpto_point = (Float64(point[1]), Float64(point[2]), Float64(point[3]))
    b.jumpto_incoming_tangent   = isnothing(incoming_tangent)   ? nothing :
        (Float64(incoming_tangent[1]),   Float64(incoming_tangent[2]),   Float64(incoming_tangent[3]))
    b.jumpto_incoming_curvature = isnothing(incoming_curvature) ? nothing :
        (Float64(incoming_curvature[1]), Float64(incoming_curvature[2]), Float64(incoming_curvature[3]))
    b.jumpto_min_bend_radius    = isnothing(min_bend_radius) ? nothing : Float64(min_bend_radius)
    b.jumpto_meta = Vector{AbstractMeta}(meta)
    return b
end

"""
    seal!(builder; extra=0.0, meta=AbstractMeta[])

Seal the Subpath at its **natural exit** — the position and tangent at the end
of the last interior segment — with no terminal connector bending. This is the
seal to use when a path should end exactly as authored, rather than being bent
toward a global target with `jumpto!`.

`extra > 0` appends a straight lead-out of `extra` meters along the natural exit
tangent. `extra == 0` (default) produces a true zero-length terminal connector.

`meta` attaches annotations to the terminal connector, which is placed and
queried like any interior segment: a `Spinning` here starts a spinning run over
the lead-out (see [`Spinning`](@ref)). A zero-length seal still carries the meta,
but its run is zero-length.

Like `jumpto!`, this seals the builder: no further segments may be added, and a
builder may be sealed exactly once (by either `seal!` or `jumpto!`). Throws if
`start!` was not called, if the builder is already sealed, or if `extra < 0`.
"""
function seal!(b::SubpathBuilder; extra::Real = 0.0,
               meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    isnothing(b.start_point) &&
        throw(ArgumentError("SubpathBuilder: call start!() before seal!()"))
    (!isnothing(b.jumpto_point) || b.jumpto_natural) &&
        throw(ArgumentError("SubpathBuilder: subpath already sealed"))
    extra < 0.0 &&
        throw(ArgumentError("SubpathBuilder: seal!() extra must be non-negative; got $extra"))
    b.jumpto_natural = true
    b.jumpto_natural_extra = Float64(extra)
    b.jumpto_meta = Vector{AbstractMeta}(meta)
    return b
end

# -----------------------------------------------------------------------
# Subpath (immutable user-supplied snapshot)
# -----------------------------------------------------------------------

"""
    Subpath

Immutable Subpath specification — user-supplied data only. Geometry queries on
a `Subpath` (arc_length, curvature, position, ...) throw with "call build()
first"; build the Subpath into a `SubpathBuilt` to query.
"""
struct Subpath
    meta::Vector{AbstractMeta}
    start_point::NTuple{3, Float64}
    start_outgoing_tangent::NTuple{3, Float64}
    start_outgoing_curvature::NTuple{3, Float64}
    segments::Vector{AbstractPathSegment}
    jumpto_point::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_tangent::Union{Nothing, NTuple{3, Float64}}
    jumpto_incoming_curvature::Union{Nothing, NTuple{3, Float64}}
    jumpto_min_bend_radius::Union{Nothing, Float64}
    # Meta on the terminal connector, carried verbatim and interpreted only by
    # consuming layers (e.g. the fiber's thermal expansion).
    jumpto_meta::Vector{AbstractMeta}
    # Natural seal (set by seal!): when true, jumpto_point is nothing and the
    # terminal connector is built directly at the natural exit (see build).
    jumpto_natural::Bool
    jumpto_natural_extra::Float64
    # `:inherit` start-state intent: when a flag is true the corresponding
    # start_* field is a placeholder, resolved at vector build time from the
    # predecessor Subpath's endpoint. The start_* tuples stay concrete so the
    # standalone build and fiber reconstruction remain type-stable.
    inherit_start_point::Bool
    inherit_start_tangent::Bool
    inherit_start_curvature::Bool
end

"""
    Subpath(b::SubpathBuilder) -> Subpath

Freeze a started-and-sealed builder into an immutable [`Subpath`](@ref) spec.

Throws if `start!` was not called, or if the builder was not sealed by
`jumpto!` or `seal!`.
"""
function Subpath(b::SubpathBuilder)
    isnothing(b.start_point) &&
        throw(ArgumentError("Subpath: builder has no start; call start!() first"))
    (isnothing(b.jumpto_point) && !b.jumpto_natural) &&
        throw(ArgumentError("Subpath: builder is not sealed; call jumpto!() or seal!() " *
                            "before constructing"))
    return Subpath(deepcopy(b.meta),
                   b.start_point::NTuple{3, Float64},
                   b.start_outgoing_tangent::NTuple{3, Float64},
                   b.start_outgoing_curvature::NTuple{3, Float64},
                   deepcopy(b.segments),
                   b.jumpto_point,
                   b.jumpto_incoming_tangent,
                   b.jumpto_incoming_curvature,
                   b.jumpto_min_bend_radius,
                   deepcopy(b.jumpto_meta),
                   b.jumpto_natural,
                   b.jumpto_natural_extra,
                   b.inherit_start_point,
                   b.inherit_start_tangent,
                   b.inherit_start_curvature)
end

# Geometry queries on an unbuilt Subpath fail loudly.
arc_length(::Subpath)               = error("Subpath: call build(subpath) before querying arc_length")
arc_length(::Subpath, ::Real, ::Real) = error("Subpath: call build(subpath) before querying arc_length")
curvature(::Subpath, ::Real)        = error("Subpath: call build(subpath) before querying curvature")
geometric_torsion(::Subpath, ::Real) = error("Subpath: call build(subpath) before querying geometric_torsion")
spinning_rate(::Subpath, ::Real)   = error("Subpath: call build(subpath) before querying spinning_rate")
position(::Subpath, ::Real)         = error("Subpath: call build(subpath) before querying position")
tangent(::Subpath, ::Real)          = error("Subpath: call build(subpath) before querying tangent")
normal(::Subpath, ::Real)           = error("Subpath: call build(subpath) before querying normal")
binormal(::Subpath, ::Real)         = error("Subpath: call build(subpath) before querying binormal")
frame(::Subpath, ::Real)            = error("Subpath: call build(subpath) before querying frame")

# -----------------------------------------------------------------------
# SubpathBuilt and PathBuilt
# -----------------------------------------------------------------------

"""
    SubpathBuilt

Built form of a `Subpath`. Contains:

- `subpath` — the source spec.
- `placed_segments` — the placed **interior** segments only.
- `jumpto_quintic_connector` — the resolved terminal connector segment.
- `jumpto_placed` — the `PlacedSegment` wrapper for the terminal connector
  (s offset, global origin, global frame at its start). Carried so query
  functions can treat the terminal connector uniformly with interior
  segments.
- `resolved_spinning` — spinning runs resolved within this Subpath.
- `pending_continuous_first_spinning` — true if the first spinning anchor was
  deferred for `build(::Vector{SubpathBuilt})` to fill in.

Local arc length runs from `0` to `s_end(::SubpathBuilt)` (computed on demand
from the placed segments + terminal connector).
"""
struct SubpathBuilt
    subpath::Subpath
    placed_segments::Vector{PlacedSegment}        # interior only
    jumpto_quintic_connector::QuinticConnector    # terminal connector
    jumpto_placed::PlacedSegment                  # placement of the terminal connector
    resolved_spinning::Vector{ResolvedSpinningRate}
    pending_continuous_first_spinning::Bool
end

"""
    PathBuilt(subpaths::Vector{SubpathBuilt})

Ordered container of independent `SubpathBuilt`s. Global s offsets and the
total `s_end` are computed on demand to avoid consistency hazards.
"""
struct PathBuilt
    subpaths::Vector{SubpathBuilt}
end

# -----------------------------------------------------------------------
# build()
# -----------------------------------------------------------------------

"""
    _safe_normalize(v) -> AbstractVector

Return `v` scaled to unit length. Assumes `v` is non-zero.
"""
_safe_normalize(v::AbstractVector) = v ./ sqrt(sum(abs2, v))

"""
    _initial_frame_from_tangent(T_tuple) -> (T, N, B)

Derive an orthonormal start frame from a tangent direction.

Picks `N` orthogonal to `T` via Gram-Schmidt against the world axis least
aligned with `T` (guaranteeing a non-degenerate cross product), then
`B = T × N`. Throws if the tangent is zero.
"""
function _initial_frame_from_tangent(T_tuple::NTuple{3, Float64})
    T = collect(T_tuple)
    Tn = sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    Tn > 0.0 || throw(ArgumentError("start_outgoing_tangent must be non-zero"))
    T = T ./ Tn
    # Pick the world axis least aligned with T to guarantee a non-degenerate cross.
    aT = (abs(T[1]), abs(T[2]), abs(T[3]))
    ref = if aT[1] <= aT[2] && aT[1] <= aT[3]
        [1.0, 0.0, 0.0]
    elseif aT[2] <= aT[3]
        [0.0, 1.0, 0.0]
    else
        [0.0, 0.0, 1.0]
    end
    Nraw = ref - dot(ref, T) .* T
    N = _safe_normalize(Nraw)
    B = cross(T, N)
    return (T, N, B)
end

"""
    build(builder::SubpathBuilder; perturb=false, jumpto_target_length=nothing) -> SubpathBuilt

Compile a started-and-sealed `SubpathBuilder` into a [`SubpathBuilt`](@ref) by
freezing it to a [`Subpath`](@ref) first. See [`build(::Subpath)`](@ref) for
`perturb` and `jumpto_target_length`.
"""
build(b::SubpathBuilder; perturb::Bool = false, jumpto_target_length = nothing) =
    build(Subpath(b); perturb, jumpto_target_length)

"""
    build(sub::Subpath; perturb=false, jumpto_target_length=nothing) -> SubpathBuilt

Resolve a `Subpath` into the global coordinate system, including all
geometry-related metas, returning a [`SubpathBuilt`](@ref).

With `perturb=false` (default) the geometry is nominal: all meta is stored but
none is applied. With `perturb=true` each authored segment's field-level
`MCMadd`/`MCMmul` are applied via [`_apply_field_mcm`](@ref) before placement;
meta this layer does not recognize (symbols that are not a segment's own fields)
is carried through untouched.

`jumpto_target_length`, when given, constrains the terminal `jumpto!` connector's
arc length to that value while still landing at the fixed `jumpto_point` (a plain
length, layer-agnostic). It is supplied by a consuming layer — the fiber uses it
to thermally expand the connector (issue #33). It combines with
`jumpto_min_bend_radius`: when a target length is set the solver picks the handle
scale by arc length and validates peak curvature against the radius limit.

# Extended help

Every segment defines its geometry in its own local frame, where it always
starts at the origin heading along local `+z`. Each segment type answers two
questions in local coordinates:

1. `end_position_local(seg)` — where does this segment end, relative to its own
   start?
2. `end_frame_local(seg)` — what is the `(T, N, B)` frame at its end, in its own
   local axes?

The placement loop rotates those local answers into global coordinates and
chains them (the "frame propagation" steps below). The frame matrix (columns
`[N | B | T]`) is the rotation from local to global, so `frame * v_local`
converts any local vector to global.

Note: this propagated frame is a parallel-transport-style frame (continuous
across joints), not the textbook Frenet frame, which flips discontinuously at
inflection points where curvature → 0. For a fiber path the continuous frame is
the desired one, as done here.
"""
function build(sub::Subpath; perturb::Bool = false, jumpto_target_length = nothing)::SubpathBuilt
    (sub.inherit_start_point || sub.inherit_start_tangent || sub.inherit_start_curvature) &&
        throw(ArgumentError("build(Subpath): start state uses :inherit but there is no " *
                            "predecessor; :inherit is only valid for a non-first Subpath " *
                            "in build([...])"))
    pos = collect(sub.start_point)
    T_frame, N_frame, B_frame = _initial_frame_from_tangent(sub.start_outgoing_tangent)
    K_in_global = collect(sub.start_outgoing_curvature)

    # Under perturb, apply each authored segment's field-level MCM before
    # placement; meta whose symbol is not one of a segment's own fields (foreign
    # annotations) is carried through untouched. Nominal build leaves segments
    # exactly as authored.
    segs = perturb ? map(_apply_field_mcm, sub.segments) : sub.segments

    s_eff  = 0.0
    placed = PlacedSegment[]  # in global coordinate system

    for seg_orig in segs
        frame      = hcat(N_frame, B_frame, T_frame)   # columns: [N | B | T]
        seg_placed = _resolve_at_placement(seg_orig, pos, frame, K_in_global)
        push!(placed, PlacedSegment(seg_placed, s_eff, copy(pos), copy(frame)))

        # Advance position and frame
        # Note: Each of Straight, Bend, Catenary, Helix and QuniticConnector has its own 
        # end_position_local() and end_frame_local() methods.
        pos_end_local         = end_position_local(seg_placed)
        (T_end_l, N_end_l, _) = end_frame_local(seg_placed)

        # Frame propagation
        pos     = pos + frame * pos_end_local
        T_frame = _safe_normalize(frame * T_end_l)
        N_end_g = frame * N_end_l
        N_frame = _safe_normalize(N_end_g - dot(N_end_g, T_frame) * T_frame)
        B_frame = cross(T_frame, N_frame)

        # Update K_in_global for the next segment: κ at end times the unit
        # normal at the end, both expressed globally.
        L_seg = arc_length(seg_placed)
        κ_end = curvature(seg_placed, L_seg)
        K_in_global = κ_end .* N_frame

        s_eff += L_seg
    end

    # Resolve the terminal connector inline.
    frame = hcat(N_frame, B_frame, T_frame)
    if sub.jumpto_natural
        # Natural seal: no target to reach. Build the terminal connector
        # directly at the natural exit, bypassing the quintic solver (a
        # coincident-endpoint solve is ill-conditioned). `extra > 0` gives a
        # straight lead-out along the local tangent ẑ; `extra == 0` gives a
        # true zero-length connector. Both are degenerate-safe in the
        # QuinticConnector query path (zero-speed → local ẑ tangent).
        connector = _build_straight_connector(sub.jumpto_natural_extra,
                                              eltype(K_in_global);
                                              meta = sub.jumpto_meta)
    else
        # Destination is global (sub.jumpto_point); transform to local frame to
        # call _build_quintic_connector.
        p1_local  = frame' * (collect(sub.jumpto_point) .- pos)
        chord     = norm(p1_local)
        t_hat_out = isnothing(sub.jumpto_incoming_tangent) ?
            (chord > 1e-15 ? p1_local ./ chord : [0.0, 0.0, 1.0]) :
            _safe_normalize(frame' * _safe_normalize(collect(sub.jumpto_incoming_tangent)))
        K0_local = frame' * K_in_global
        K1_local = isnothing(sub.jumpto_incoming_curvature) ?
            zeros(eltype(K0_local), 3) :
            frame' * collect(sub.jumpto_incoming_curvature)
        # The caller may constrain the terminal connector's arc length (the fiber
        # supplies this to thermally expand the connector — issue #33).
        # `min_bend_radius` is always honored: with a target set the
        # solver picks the handle by arc length and validates the radius limit
        # post-hoc; with no target it drives the handle selection.
        connector = _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                                             min_bend_radius    = sub.jumpto_min_bend_radius,
                                             target_path_length = jumpto_target_length,
                                             meta               = sub.jumpto_meta)
    end
    # PlacedSegment wrapper for the terminal connector: anchor at the
    # position/frame at the end of the interior segments. Stored alongside
    # the connector so query functions treat it like any other placed segment.
    jumpto_placed = PlacedSegment(connector, s_eff, copy(pos), copy(frame))
    L_conn = arc_length(connector)
    s_end_eff = s_eff + L_conn

    # Resolve spinnings local to this Subpath. The terminal connector is included
    # as the last placed segment so a `Spinning` on `jumpto!`/`seal!` anchors a run
    # exactly like any interior segment. If the first spinning anchor has
    # is_continuous=true it is left pending for the PathBuilt build to fix up.
    resolved, pending = _resolve_spinning_subpath_local(
        vcat(placed, jumpto_placed), Float64(_qc_nominalize(s_end_eff)))
    return SubpathBuilt(sub, placed, connector, jumpto_placed, resolved, pending)
end

# -----------------------------------------------------------------------
# Spinning resolution (per Subpath)
# -----------------------------------------------------------------------
# Walk placed segments in order (the terminal connector is the last one),
# collect Spinning anchors from each segment's meta, and emit one
# ResolvedSpinningRate per anchor. Each run extends from the anchor's segment
# start to the next anchor's segment start, or to the end of the Subpath.
# is_continuous=true on a non-first anchor takes its phi_0 from the prior run's
# accumulated phase. is_continuous=true on the *first* anchor of a Subpath
# leaves resolution pending — `build(::Vector{SubpathBuilt})` fixes it up by
# inheriting from the prior Subpath.

"""
    _collect_spinning_anchors(placed) -> Vector{Tuple{Float64, Spinning}}

Collect `(s_offset_eff, Spinning)` anchors from the placed segments' meta, in
order. The terminal connector is an ordinary entry in `placed` (see
[`_all_placed_segs`](@ref)) and is treated no differently from an interior
segment.

Throws if any single segment carries more than one `Spinning`.
"""
function _collect_spinning_anchors(placed::Vector{PlacedSegment})
    anchors = Tuple{Float64, Spinning}[]
    for ps in placed
        spinnings_here = Spinning[]
        for m in segment_meta(ps.segment)
            m isa Spinning && push!(spinnings_here, m)
        end
        if length(spinnings_here) > 1
            throw(ArgumentError(
                "Subpath build: segment at s_offset_eff = $(ps.s_offset_eff) carries " *
                "$(length(spinnings_here)) Spinning meta entries; at most one is permitted"))
        end
        if !isempty(spinnings_here)
            push!(anchors, (Float64(_qc_nominalize(ps.s_offset_eff)), spinnings_here[1]))
        end
    end
    return anchors
end

"""
    _resolve_spinning_subpath_local(placed, s_end)
        -> (Vector{ResolvedSpinningRate}, pending_first::Bool)

Resolve spinning runs for a single Subpath. `placed` is the full ordered list of
placed segments including the terminal connector (see [`_all_placed_segs`](@ref)).

Each run extends from its anchor's segment start to the next anchor's start, or
to `s_end`. A non-first `is_continuous` anchor inherits `phi_0` from the prior
run's accumulated phase. If the *first* anchor is `is_continuous`, its run is
emitted with `phi_0 = NaN` as a sentinel and `pending_first = true`; the
`PathBuilt`-level resolver later replaces it with the inherited value.
"""
function _resolve_spinning_subpath_local(placed::Vector{PlacedSegment},
                                       s_end::Float64)
    anchors = _collect_spinning_anchors(placed)
    isempty(anchors) && return (ResolvedSpinningRate[], false)

    n = length(anchors)
    out = Vector{ResolvedSpinningRate}(undef, n)
    prev_phi_0 = 0.0
    prev_run_length = 0.0
    prev_rate::Union{Float64, Function} = 0.0
    pending_first = false

    for i in 1:n
        s_start_i, tw = anchors[i]
        s_run_end = (i < n) ? anchors[i + 1][1] : s_end

        if tw.is_continuous
            if i == 1
                # Defer to PathBuilt-level resolution.
                pending_first = true
                phi_0 = NaN
            else
                phi_0 = prev_phi_0 + _integrate_rate(prev_rate, 0.0, prev_run_length)
            end
        else
            phi_0 = tw.phi_0
        end

        out[i] = ResolvedSpinningRate(s_start_i, s_run_end, tw.rate, phi_0)

        prev_phi_0 = phi_0
        prev_run_length = s_run_end - s_start_i
        prev_rate = tw.rate
    end

    return (out, pending_first)
end

"""
    _resolve_pending_continuous_spinning(builts) -> Vector{SubpathBuilt}

Fix up cross-Subpath continuous spinning phase across an ordered `builts` list.

Where a Subpath has `pending_continuous_first_spinning = true`, look up the
prior Subpath's terminal spinning phase and rebuild its `resolved_spinning`
list with the inherited `phi_0`. The first Subpath cannot be pending; throw if
it is.
"""
function _resolve_pending_continuous_spinning(builts::Vector{SubpathBuilt})
    n = length(builts)
    out = Vector{SubpathBuilt}(undef, n)
    for i in 1:n
        b = builts[i]
        if !b.pending_continuous_first_spinning
            out[i] = b
            continue
        end
        if i == 1
            throw(ArgumentError(
                "PathBuilt: first Subpath has Spinning(is_continuous=true) on its first " *
                "anchor, but there is no prior Subpath to inherit phase from"))
        end
        prev = out[i - 1]
        # Compute the prior Subpath's terminal phase: the last resolved spinning
        # run's phi_0 plus integral over its length.
        if isempty(prev.resolved_spinning)
            throw(ArgumentError(
                "PathBuilt: Subpath $i has pending continuous first spinning but " *
                "Subpath $(i-1) has no spinning runs to inherit from"))
        end
        last_run = prev.resolved_spinning[end]
        run_len  = last_run.s_eff_end - last_run.s_eff_start
        phi_end  = last_run.phi_0 + _integrate_rate(last_run.rate, 0.0, run_len)

        # Rebuild this Subpath's resolved_spinning list, replacing the first run's phi_0.
        rs_old = b.resolved_spinning
        rs_new = Vector{ResolvedSpinningRate}(undef, length(rs_old))
        # First run: inherit phi_end. Then propagate forward through subsequent
        # runs that were resolved with prev_phi_0 originating from the (NaN)
        # first run — re-run the phi_0 chain from this corrected start.
        prev_phi_0 = phi_end
        prev_run_length = 0.0
        prev_rate::Union{Float64, Function} = 0.0
        for k in eachindex(rs_old)
            r = rs_old[k]
            if k == 1
                phi_0 = phi_end
            else
                # Re-derive: was tw.is_continuous? We don't have the raw Spinning
                # object handy, so detect by NaN propagation in the original
                # resolution. The simple rule: if the original phi_0 isn't NaN
                # and was deterministic (from tw.phi_0), keep it. If it was
                # derived from prev_phi_0 chain (k > 1 and is_continuous), it
                # would have been computed without the corrected base.
                # To be safe, recompute: phi_0 = prev_phi_0 + ∫ prev_rate over prev_run_length
                # only if the original is "close" to that expression. We
                # cannot disambiguate without the source Spinning meta, so we
                # walk the segments again below.
                phi_0 = r.phi_0
            end
            rs_new[k] = ResolvedSpinningRate(r.s_eff_start, r.s_eff_end, r.rate, phi_0)
            prev_phi_0 = phi_0
            prev_run_length = r.s_eff_end - r.s_eff_start
            prev_rate = r.rate
        end

        # If the first run's phi_0 was NaN and there are downstream is_continuous
        # runs that inherited from it, those would also have been computed with
        # prev_phi_0=NaN and need recomputation. Re-resolve from raw anchors with
        # phi_end as the seed.
        if any(r -> isnan(r.phi_0), rs_old) && length(rs_old) >= 2
            anchors = _collect_spinning_anchors(_all_placed_segs(b))
            # Recompute with first phi_0 = phi_end (forced is_continuous semantics)
            n_a = length(anchors)
            rs_new = Vector{ResolvedSpinningRate}(undef, n_a)
            prev_phi_0 = phi_end
            prev_run_length = 0.0
            prev_rate = 0.0
            # s_end for the runs is the same as the original
            s_end = isempty(rs_old) ? 0.0 : rs_old[end].s_eff_end
            for k in 1:n_a
                s_start_k, tw = anchors[k]
                s_run_end = (k < n_a) ? anchors[k + 1][1] : s_end
                if tw.is_continuous
                    if k == 1
                        phi_0 = phi_end
                    else
                        phi_0 = prev_phi_0 + _integrate_rate(prev_rate, 0.0, prev_run_length)
                    end
                else
                    phi_0 = tw.phi_0
                end
                rs_new[k] = ResolvedSpinningRate(s_start_k, s_run_end, tw.rate, phi_0)
                prev_phi_0 = phi_0
                prev_run_length = s_run_end - s_start_k
                prev_rate = tw.rate
            end
        end

        out[i] = SubpathBuilt(b.subpath, b.placed_segments, b.jumpto_quintic_connector,
                              b.jumpto_placed, rs_new, false)
    end
    return out
end

# -----------------------------------------------------------------------
# build(::Vector{Subpath} | ::Vector{SubpathBuilt} | ::SubpathBuilt) → PathBuilt
# -----------------------------------------------------------------------

"""
    _tuple_isapprox(a, b; atol=1e-9, rtol=1e-9) -> Bool

Return whether two length-3 tuples are elementwise approximately equal (used
for endpoint conformity checks).
"""
_tuple_isapprox(a::NTuple{3, Float64}, b::NTuple{3, Float64};
                atol::Float64 = 1e-9, rtol::Float64 = 1e-9) =
    all(isapprox(a[i], b[i]; atol = atol, rtol = rtol) for i in 1:3)

"""
    _subpath_endpoint_state(prev_built) -> (point, tangent, curvature)

Return the endpoint state of a built Subpath as
`(point::NTuple{3,Float64}, tangent::Union{Nothing,NTuple{3,Float64}},
curvature::Union{Nothing,NTuple{3,Float64}})`.

For a naturally-sealed predecessor (`jumpto_natural`), the point and tangent are
read from the built geometry (nominalized) and there is no declared end curvature
(`nothing`). For a `jumpto!`-sealed predecessor they come from the declared
`jumpto_point` / `jumpto_incoming_tangent` / `jumpto_incoming_curvature` (the
latter two may be `nothing`).
"""
function _subpath_endpoint_state(prev_built::SubpathBuilt)
    prev = prev_built.subpath
    if prev.jumpto_natural
        s_e        = Float64(_qc_nominalize(s_end(prev_built)))
        ep_point   = Tuple(Float64.(position(prev_built, s_e)))::NTuple{3, Float64}
        ep_tangent = Tuple(Float64.(tangent(prev_built, s_e)))::NTuple{3, Float64}
        return (ep_point, ep_tangent, nothing)
    else
        return (prev.jumpto_point::NTuple{3, Float64},
                prev.jumpto_incoming_tangent,
                prev.jumpto_incoming_curvature)
    end
end

"""
    _resolve_inherited_start(cur, prev_built) -> Subpath

Resolve a Subpath's `:inherit` start fields from the predecessor's endpoint.
Returns `cur` unchanged when no inherit flag is set. Otherwise builds a new
`Subpath` with each flagged field replaced by a concrete value:

- `point` ← the predecessor endpoint point;
- `outgoing_tangent` ← the declared incoming tangent, or — when the predecessor
  was `jumpto!`-sealed with a chord-direction default (`nothing`) — the actual
  exit tangent queried from the built geometry;
- `outgoing_curvature` ← the declared incoming curvature, or `(0,0,0)` when none
  was declared (matching conformity, which only ever checks declared curvature).

The inherit flags are cleared so the resolved Subpath builds standalone.
"""
function _resolve_inherited_start(cur::Subpath, prev_built::SubpathBuilt)
    (cur.inherit_start_point || cur.inherit_start_tangent ||
     cur.inherit_start_curvature) || return cur
    point, tangent, curvature = _subpath_endpoint_state(prev_built)
    new_point = cur.inherit_start_point ? point : cur.start_point
    new_tangent = if cur.inherit_start_tangent
        # Chord-direction default (nothing) is only realized once the connector is
        # built, so query the built exit tangent for a concrete vector.
        isnothing(tangent) ?
            Tuple(Float64.(end_tangent(prev_built)))::NTuple{3, Float64} : tangent
    else
        cur.start_outgoing_tangent
    end
    new_curvature = if cur.inherit_start_curvature
        isnothing(curvature) ? (0.0, 0.0, 0.0) : curvature
    else
        cur.start_outgoing_curvature
    end
    return Subpath(cur.meta, new_point, new_tangent, new_curvature, cur.segments,
                   cur.jumpto_point, cur.jumpto_incoming_tangent,
                   cur.jumpto_incoming_curvature, cur.jumpto_min_bend_radius,
                   cur.jumpto_meta, cur.jumpto_natural, cur.jumpto_natural_extra,
                   false, false, false)
end

"""
    _check_subpath_conformity(prev_built, cur, idx)

Throw unless the endpoint state of `prev_built` matches the start state of the
`cur` Subpath (Subpaths are independent in spec but must stack in order).

`prev_built` is the built form so a naturally-sealed predecessor
(`jumpto_point === nothing`) can have its endpoint position/tangent read from
the built geometry rather than from a (nonexistent) global jumpto spec.
"""
function _check_subpath_conformity(prev_built::SubpathBuilt, cur::Subpath, idx::Int)
    jump_point, jump_tangent, pk = _subpath_endpoint_state(prev_built)
    if !_tuple_isapprox(jump_point, cur.start_point)
        throw(ArgumentError(
            "PathBuilt: Subpath $(idx-1) endpoint $(jump_point) " *
            "does not match Subpath $idx start_point $(cur.start_point)"))
    end
    # jumpto_incoming_tangent default is "chord direction"; nothing to compare
    # if prev has nothing. If both supplied, compare; if one supplied, that's a
    # mismatch.
    pt = jump_tangent
    ct = cur.start_outgoing_tangent
    if !isnothing(pt)
        if !_tuple_isapprox(pt, ct)
            throw(ArgumentError(
                "PathBuilt: Subpath $(idx-1) endpoint tangent $pt " *
                "does not match Subpath $idx start_outgoing_tangent $ct"))
        end
    end
    ck = cur.start_outgoing_curvature
    if !isnothing(pk)
        if !_tuple_isapprox(pk, ck)
            throw(ArgumentError(
                "PathBuilt: Subpath $(idx-1) jumpto_incoming_curvature $pk " *
                "does not match Subpath $idx start_outgoing_curvature $ck"))
        end
    end
    return nothing
end

"""
    build(builts::Vector{SubpathBuilt})   → PathBuilt
    build(subpaths::Vector{Subpath})      → PathBuilt
    build(builders::Vector{SubpathBuilder}) → PathBuilt
    build(spb::SubpathBuilt)              → PathBuilt

Stitch already-built `SubpathBuilt`s into a `PathBuilt`, validating that
adjacent Subpaths' endpoint states agree, and resolving any cross-Subpath
spinning continuity. The vector-of-Subpath form builds each first; the
`Vector{SubpathBuilder}` convenience form freezes each builder to a `Subpath`
before building; the single `SubpathBuilt` form wraps a length-1 PathBuilt.
"""
function build(builts::Vector{SubpathBuilt})
    isempty(builts) && throw(ArgumentError("PathBuilt: at least one SubpathBuilt required"))
    for i in 2:length(builts)
        _check_subpath_conformity(builts[i-1], builts[i].subpath, i)
    end
    fixed = _resolve_pending_continuous_spinning(builts)
    return PathBuilt(fixed)
end

# Build a vector of Subpaths sequentially so that each Subpath's `:inherit` start
# fields can be resolved from the predecessor's already-built endpoint.
# `_resolve_inherited_start` is a no-op when no inherit flag is set.
function build(subpaths::Vector{Subpath}; perturb::Bool = false)
    isempty(subpaths) && throw(ArgumentError("PathBuilt: at least one Subpath required"))
    first_sub = subpaths[1]
    (first_sub.inherit_start_point || first_sub.inherit_start_tangent ||
     first_sub.inherit_start_curvature) &&
        throw(ArgumentError("build([...]): :inherit is not allowed on the first Subpath " *
                            "(no predecessor to inherit from)"))
    builts = Vector{SubpathBuilt}(undef, length(subpaths))
    builts[1] = build(first_sub; perturb)
    for i in 2:lastindex(subpaths)
        resolved = _resolve_inherited_start(subpaths[i], builts[i-1])
        builts[i] = build(resolved; perturb)
    end
    return build(builts)
end

build(builders::Vector{SubpathBuilder}; perturb::Bool = false) =
    build(Subpath[Subpath(b) for b in builders]; perturb = perturb)

build(spb::SubpathBuilt) = build(SubpathBuilt[spb])


# -----------------------------------------------------------------------
# Segment lookup helpers
# -----------------------------------------------------------------------

"""
    _find_placed_segment(b::SubpathBuilt, s) -> (PlacedSegment, s_local)

Locate the placed segment (interior or terminal connector) containing global
arc length `s`, returning it with `s_local` clamped into the segment's domain.
"""
function _find_placed_segment(b::SubpathBuilt, s)
    s_eff = 0.0
    for ps in b.placed_segments
        seg_len = arc_length(ps.segment)
        seg_len_nom = Float64(_qc_nominalize(seg_len))
        s_eff_next  = s_eff + seg_len_nom
        if s <= s_eff_next + 1e-12
            s_local = clamp(s - s_eff, zero(seg_len), seg_len)
            return ps, s_local
        end
        s_eff = s_eff_next
    end
    # Past all interior segments → terminal connector.
    ps_t      = b.jumpto_placed
    seg_len_t = arc_length(ps_t.segment)
    s_local   = clamp(s - s_eff, zero(seg_len_t), seg_len_t)
    return ps_t, s_local
end

"""
    _local_to_global(ps::PlacedSegment, v_local) -> AbstractVector

Rotate a local-frame vector into global coordinates via the placed segment's
frame (`v_g = frame * v_l`).
"""
function _local_to_global(ps::PlacedSegment, v_local::AbstractVector)
    return ps.frame * v_local
end

"""
    s_end(b::SubpathBuilt) -> T

Return the Subpath's total arc length: the sum of the interior segment lengths
plus the terminal connector length (computed on demand).
"""
function s_end(b::SubpathBuilt)
    total = zero(arc_length(b.jumpto_quintic_connector))
    for ps in b.placed_segments
        total = total + arc_length(ps.segment)
    end
    return total + arc_length(b.jumpto_quintic_connector)
end

# -----------------------------------------------------------------------
# Differential geometry interface on SubpathBuilt
# -----------------------------------------------------------------------

"""
    arc_length(path) -> T
    arc_length(path, s1, s2) -> T

Return a `path`'s total arc length, or the arc length between `s1` and `s2`
(which is simply `s2 - s1` in the global arc-length parameterization).
Applies to both `SubpathBuilt` and `PathBuilt`.
"""
arc_length(b::SubpathBuilt) = s_end(b)

function arc_length(::SubpathBuilt, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

"""
    curvature(path, s) -> T

Return the curvature κ (1/m) of a built `path` at global arc length `s`.
"""
function curvature(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return curvature(ps.segment, s_local)
end

"""
    geometric_torsion(path, s) -> T

Return the geometric torsion τ (rad/m) of a built `path` at global arc
length `s`.
"""
function geometric_torsion(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return geometric_torsion(ps.segment, s_local)
end

"""
    spinning_rate(b, s)

Spinning rate (rad/m) at local arc length `s`, summed over all resolved
spinning runs that contain `s`. Runs are disjoint by construction; the sum
exists only as a robustness guard.
"""
function spinning_rate(b::SubpathBuilt, s)
    τ = zero(s isa AbstractFloat ? s : Float64(s))
    for r in b.resolved_spinning
        if r.s_eff_start <= s <= r.s_eff_end
            τ += r.rate isa Function ? r.rate(s - r.s_eff_start) : r.rate
        end
    end
    return τ
end

"""
    position(path, s) -> Vector

Return the global-frame position of a built `path` at arc length `s`.
"""
function position(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return ps.origin + _local_to_global(ps, position_local(ps.segment, s_local))
end

"""
    tangent(path, s) -> Vector

Return the global-frame unit tangent of a built `path` at arc length `s`.
"""
function tangent(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return _local_to_global(ps, tangent_local(ps.segment, s_local))
end

"""
    normal(path, s) -> Vector

Return the global-frame unit normal of a built `path` at arc length `s` (the
parallel-transport frame; see [`build`](@ref)).
"""
function normal(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return _local_to_global(ps, normal_local(ps.segment, s_local))
end

"""
    binormal(path, s) -> Vector

Return the global-frame unit binormal of a built `path` at arc length `s`.
"""
function binormal(b::SubpathBuilt, s::Real)
    ps, s_local = _find_placed_segment(b, s)
    return _local_to_global(ps, binormal_local(ps.segment, s_local))
end

"""
    frame(path, s) -> NamedTuple

Return all differential-geometry quantities of a built `path` at arc length
`s` as a `NamedTuple`: `position`, `tangent`, `normal`, `binormal`,
`curvature`, `geometric_torsion`, and `spinning_rate`.
"""
function frame(b::SubpathBuilt, s::Real)
    T = tangent(b, s)
    N = normal(b, s)
    Bi = binormal(b, s)
    κ = curvature(b, s)
    τ = geometric_torsion(b, s)
    m = spinning_rate(b, s)
    return (; position = position(b, s), tangent = T, normal = N, binormal = Bi,
              curvature = κ, geometric_torsion = τ, spinning_rate = m)
end

# -----------------------------------------------------------------------
# Endpoint access on SubpathBuilt
# -----------------------------------------------------------------------

"""
    start_point(path) -> Vector

Return the global-frame position at a built `path`'s start (`s = 0`).
"""
start_point(b::SubpathBuilt)   = position(b, 0.0)

"""
    end_point(path) -> Vector

Return the global-frame position at a built `path`'s end
(`s = arc_length(path)`).
"""
end_point(b::SubpathBuilt)     = position(b, Float64(_qc_nominalize(s_end(b))))

"""
    start_tangent(path) -> Vector

Return the global-frame unit tangent at a built `path`'s start (`s = 0`).
"""
start_tangent(b::SubpathBuilt) = tangent(b, 0.0)

"""
    end_tangent(path) -> Vector

Return the global-frame unit tangent at a built `path`'s end
(`s = arc_length(path)`).
"""
end_tangent(b::SubpathBuilt)   = tangent(b, Float64(_qc_nominalize(s_end(b))))

# -----------------------------------------------------------------------
# Path measures
# -----------------------------------------------------------------------

"""
    path_length(path) -> T

Return a built `path`'s total arc length. Alias for [`arc_length`](@ref).
"""
path_length(b::SubpathBuilt) = arc_length(b)

"""
    cartesian_distance(path, s1, s2) -> T

Return the straight-line (chord) distance between the points at arc lengths
`s1` and `s2` on a built `path`.
"""
function cartesian_distance(b::SubpathBuilt, s1::Real, s2::Real)
    return norm(position(b, s2) - position(b, s1))
end

"""
    bounding_box(path; n=512) -> (; lo, hi)

Return the axis-aligned bounding box `(lo, hi)` of a built `path`, sampled at
`n` points along its arc length.
"""
function bounding_box(b::SubpathBuilt; n::Int = 512)
    s0 = 0.0
    s1 = Float64(_qc_nominalize(s_end(b)))
    ss = range(s0, s1; length = n)
    pts = [position(b, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

"""
    _all_placed_segs(b::SubpathBuilt) -> Vector{PlacedSegment}

Return the interior placed segments followed by the terminal connector. Used
by the `total_*` iterators.
"""
function _all_placed_segs(b::SubpathBuilt)
    result = Vector{PlacedSegment}(undef, length(b.placed_segments) + 1)
    @inbounds for i in eachindex(b.placed_segments)
        result[i] = b.placed_segments[i]
    end
    result[end] = b.jumpto_placed
    return result
end

"""
    total_turning_angle(path) -> Float64

Return the integral of curvature over a built `path` (total bending angle, in
radians). Applies to `SubpathBuilt` and `PathBuilt`.
"""
function total_turning_angle(b::SubpathBuilt)
    total = 0.0
    for ps in _all_placed_segs(b)
        seg = ps.segment
        if seg isa StraightSegment
            # κ = 0
        elseif seg isa BendSegment
            total += abs(seg.angle)
        elseif seg isa HelixSegment
            total += curvature(seg, 0.0) * arc_length(seg)
        else
            n = 64
            ss = range(0.0, arc_length(seg); length = n + 1)
            h = ss[2] - ss[1]
            total += h * sum(curvature(seg, s) for s in ss)
        end
    end
    return total
end

"""
    total_torsion(path) -> Float64

Return the integral of geometric torsion over a built `path` (in radians).
Applies to `SubpathBuilt` and `PathBuilt`.
"""
function total_torsion(b::SubpathBuilt)
    total = 0.0
    for ps in _all_placed_segs(b)
        seg = ps.segment
        if seg isa HelixSegment
            total += geometric_torsion(seg, 0.0) * arc_length(seg)
        elseif seg isa StraightSegment || seg isa BendSegment || seg isa CatenarySegment
            # τ_geom = 0
        else
            n = 64
            ss = range(0.0, arc_length(seg); length = n + 1)
            h = ss[2] - ss[1]
            total += h * sum(geometric_torsion(seg, s) for s in ss)
        end
    end
    return total
end

"""
    total_spinning(b; s_start, s_end, rtol = 1e-8, atol = 0.0) → Float64

Integrated material spinning ``∫ τ_{\\mathrm{mat}}(s) \\, ds`` over local arc length
from `s_start` to `s_end` (defaults: full Subpath). Both endpoints must lie
in `[0, s_end(b)]`.
"""
function total_spinning(
    b::SubpathBuilt;
    s_start::Real = 0.0,
    s_end::Real   = s_end(b),
    rtol::Real    = 1e-8,
    atol::Real    = 0.0,
)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    if s_lo > s_hi
        throw(ArgumentError(
            "total_spinning: require s_start ≤ s_end; got s_start=$(s_lo), s_end=$(s_hi)"))
    end
    ps0 = 0.0
    ps1 = Float64(_qc_nominalize(arc_length(b)))
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_spinning: require 0 ≤ s ≤ s_end(b) for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs subpath domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for r in b.resolved_spinning
        a = max(r.s_eff_start, s_lo)
        bb = min(r.s_eff_end, s_hi)
        bb <= a && continue
        total += _integrate_rate(r.rate, a - r.s_eff_start, bb - r.s_eff_start;
                                 rtol = rtolf, atol = atolf)
    end
    return total
end

"""
    total_frame_rotation(b; s_start, s_end, rtol = 1e-8, atol = 0.0) → Float64

Total frame rotation `∫ (τ_geom + Ω_spin) ds` over local arc length from
`s_start` to `s_end`.
"""
function total_frame_rotation(
    b::SubpathBuilt;
    s_start::Real = 0.0,
    s_end::Real   = s_end(b),
    rtol::Real    = 1e-8,
    atol::Real    = 0.0,
)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    if s_lo > s_hi
        throw(ArgumentError(
            "total_frame_rotation: require s_start ≤ s_end; got s_start=$(s_lo), s_end=$(s_hi)"))
    end
    ps0 = 0.0
    ps1 = Float64(_qc_nominalize(arc_length(b)))
    if !(ps0 - 1e-12 <= s_lo <= ps1 + 1e-12) || !(ps0 - 1e-12 <= s_hi <= ps1 + 1e-12)
        throw(ArgumentError(
            "total_frame_rotation: require 0 ≤ s ≤ s_end(b) for both endpoints; " *
            "got [$(s_lo), $(s_hi)] m vs subpath domain [$(ps0), $(ps1)] m"))
    end
    s_lo == s_hi && return 0.0

    τ_total = 0.0
    rtolf = Float64(rtol)
    atolf = Float64(atol)
    for ps in _all_placed_segs(b)
        seg = ps.segment
        seg_s_start_nom = Float64(_qc_nominalize(ps.s_offset_eff))
        seg_s_end_nom   = seg_s_start_nom +
                          Float64(_qc_nominalize(arc_length(seg)))
        a_nom = max(seg_s_start_nom, s_lo)
        b_nom = min(seg_s_end_nom,   s_hi)
        b_nom <= a_nom && continue

        L_seg  = arc_length(seg)
        seg_span_nom = seg_s_end_nom - seg_s_start_nom
        frac_a = (a_nom - seg_s_start_nom) / seg_span_nom
        frac_b = (b_nom - seg_s_start_nom) / seg_span_nom
        a_loc = frac_a * L_seg
        b_loc = frac_b * L_seg

        if seg isa StraightSegment || seg isa BendSegment || seg isa CatenarySegment
            # τ_geom = 0
        elseif seg isa HelixSegment
            τ_total += geometric_torsion(seg, zero(L_seg)) * (b_loc - a_loc)
        elseif seg isa QuinticConnector
            # QuinticConnector returns geometric_torsion(::QuinticConnector, ::Real)
            # = 0.0 identically (planar treatment); skip without QuadGK so the
            # MCM Particles arc-length bounds don't trigger a kronrod method error.
        else
            val, _ = QuadGK.quadgk(s -> geometric_torsion(seg, s), a_loc, b_loc;
                                   rtol = rtolf, atol = atolf)
            τ_total += val
        end
    end

    Ω_total = total_spinning(b; s_start = s_lo, s_end = s_hi,
                                   rtol = rtolf, atol = atolf)
    return τ_total + Ω_total
end

"""
    writhe(b; n) → Float64

Writhe of the Subpath: numerical double integral on `n` samples.
"""
function writhe(b::SubpathBuilt; n::Int = 256)
    s0 = 0.0
    s1 = Float64(_qc_nominalize(s_end(b)))
    ss = collect(range(s0, s1; length = n))
    rs = [position(b, s) for s in ss]
    ts = [tangent(b, s)  for s in ss]
    ds = (s1 - s0) / (n - 1)

    Wr = 0.0
    for i in 1:n, j in 1:n
        i == j && continue
        r_ij = rs[i] - rs[j]
        d = norm(r_ij)
        d < 1e-14 * (s1 - s0) && continue
        Wr += dot(cross(ts[i], ts[j]), r_ij) / d^3
    end
    return Wr * ds^2 / (4π)
end

# -----------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------

"""
    Sample

One evaluated point on a path: arc-length coordinate `s` plus all frame
quantities (position, tangent/normal/binormal, curvature, geometric torsion,
spinning rate).

The frame is the propagated (parallel-transport) frame, not the textbook
Frenet frame — see [`build`](@ref).
"""
struct Sample
    s                 :: Real
    position          :: AbstractVector
    tangent           :: AbstractVector
    normal            :: AbstractVector
    binormal          :: AbstractVector
    curvature         :: Real
    geometric_torsion :: Real
    spinning_rate    :: Real
end

"""
    PathSample

Dense samples of a Subpath/PathBuilt over `[s_start, s_end]`.
"""
struct PathSample
    samples :: Vector{Sample}
    s_start :: Float64
    s_end   :: Float64
    n       :: Int
end

"""
    _segment_total_angle(seg) -> Float64

Return the total turning angle (radians) swept by a segment, used to budget
adaptive sample density. Defaults to `0.0` for segments with no net turning.
"""
function _segment_total_angle(seg::BendSegment)
    return abs(seg.angle)
end

function _segment_total_angle(seg::CatenarySegment)
    return atan(arc_length(seg) / seg.a)
end

function _segment_total_angle(seg::HelixSegment)
    return seg.turns * 2π
end

function _segment_total_angle(seg::QuinticConnector)
    L = arc_length(seg)
    L < 1e-15 && return 0.0
    n = 16
    ss = range(0.0, L; length = n + 1)
    h  = L / n
    return h * (sum(curvature(seg, s) for s in ss) - curvature(seg, 0.0)/2 - curvature(seg, L)/2)
end

_segment_total_angle(::AbstractPathSegment) = 0.0

"""
    _budget_scalar(x) -> Float64

Reduce a value to a plain `Float64` for sample-budgeting; for MCM `Particles`
takes the worst-case (`maximum`) so the budget covers the whole ensemble.
"""
_budget_scalar(x::AbstractFloat) = Float64(x)
_budget_scalar(x::Integer) = Float64(x)
function _budget_scalar(x)
    if hasfield(typeof(x), :particles)
        return Float64(maximum(getfield(x, :particles)))
    end
    return Float64(x)
end

"""
    _segment_point_budget(ps, b, s_lo, s_hi, fidelity) -> Int

Return the number of sample points to allocate to placed segment `ps` over the
clamped window `[s_lo, s_hi]`, scaled by `fidelity` and the segment's geometric
turning and spinning angle.
"""
function _segment_point_budget(
    ps::PlacedSegment,
    b::SubpathBuilt,
    s_lo::Float64,
    s_hi::Float64,
    fidelity::Float64,
)
    seg = ps.segment
    seg_s_start = Float64(_qc_nominalize(ps.s_offset_eff))
    seg_s_end   = seg_s_start + Float64(_qc_nominalize(arc_length(seg)))

    a = max(s_lo, seg_s_start)
    bb = min(s_hi, seg_s_end)
    bb <= a && return 2

    seg_len = seg_s_end - seg_s_start
    frac = seg_len > 0.0 ? (bb - a) / seg_len : 1.0

    geom_angle  = _budget_scalar(_segment_total_angle(seg) * frac)
    geom_budget = max(2, ceil(Int, fidelity * geom_angle / (2π) * 32))

    spinning_total  = total_spinning(b; s_start = a, s_end = bb, rtol = 1e-3)
    spinning_angle  = abs(_budget_scalar(spinning_total))
    spinning_budget = max(2, ceil(Int, fidelity * spinning_angle / (2π) * 32))

    return max(geom_budget, spinning_budget)
end

"""
    sample_path(b, s1, s2; fidelity = 1.0) → PathSample

Adaptive sampling of a `SubpathBuilt` over `[s1, s2]`.
"""
function sample_path(b::SubpathBuilt, s1::Real, s2::Real; fidelity::Float64 = 1.0)
    @assert s2 > s1    "sample_path: require s2 > s1"
    @assert fidelity > 0.0 "sample_path: fidelity must be positive"

    s_lo = Float64(_qc_nominalize(s1))
    s_hi = Float64(_qc_nominalize(s2))

    all_s = Float64[]
    for ps in _all_placed_segs(b)
        seg_s_start = Float64(_qc_nominalize(ps.s_offset_eff))
        seg_s_end   = seg_s_start + Float64(_qc_nominalize(arc_length(ps.segment)))

        a = max(s_lo, seg_s_start)
        bb = min(s_hi, seg_s_end)
        bb <= a && continue

        n_seg = _segment_point_budget(ps, b, s_lo, s_hi, fidelity)
        seg_ss = collect(range(a, bb; length = n_seg))

        if isempty(all_s)
            append!(all_s, seg_ss)
        else
            start_idx = (seg_ss[1] ≈ all_s[end]) ? 2 : 1
            append!(all_s, @view seg_ss[start_idx:end])
        end
    end

    if isempty(all_s)
        all_s = [s_lo, s_hi]
    elseif length(all_s) == 1
        push!(all_s, s_hi)
    end

    n = length(all_s)
    samples = Vector{Sample}(undef, n)
    for i in eachindex(all_s)
        fr = frame(b, all_s[i])
        samples[i] = Sample(
            all_s[i],
            fr.position,
            fr.tangent,
            fr.normal,
            fr.binormal,
            fr.curvature,
            fr.geometric_torsion,
            fr.spinning_rate,
        )
    end
    return PathSample(samples, s_lo, s_hi, n)
end

# -----------------------------------------------------------------------
# Breakpoints
# -----------------------------------------------------------------------

"""
    normalize_breakpoints(breakpoints) -> Vector

Return the breakpoints sorted and deduplicated.
"""
function normalize_breakpoints(breakpoints::AbstractVector{<:Real})
    return sort(unique(copy(breakpoints)))
end

"""
    path_segment_breakpoints(b::SubpathBuilt) -> Vector{Float64}

Return the normalized arc-length breakpoints at every segment boundary
(including the terminal connector and the path endpoints).
"""
function path_segment_breakpoints(b::SubpathBuilt)
    points = Float64[0.0]
    for ps in b.placed_segments
        push!(points, Float64(_qc_nominalize(ps.s_offset_eff)))
        push!(points, Float64(_qc_nominalize(
            ps.s_offset_eff + arc_length(ps.segment))))
    end
    push!(points, Float64(_qc_nominalize(b.jumpto_placed.s_offset_eff)))
    push!(points, Float64(_qc_nominalize(arc_length(b))))
    return normalize_breakpoints(points)
end

"""
    path_spinning_breakpoints(b::SubpathBuilt) -> Vector{Float64}

Return the normalized arc-length breakpoints at every spinning-run boundary
(plus the path endpoints).
"""
function path_spinning_breakpoints(b::SubpathBuilt)
    points = Float64[0.0, Float64(_qc_nominalize(arc_length(b)))]
    for r in b.resolved_spinning
        push!(points, r.s_eff_start)
        push!(points, r.s_eff_end)
    end
    return normalize_breakpoints(points)
end

"""
    breakpoints(path) -> Vector{Float64}

Return the merged, normalized arc-length breakpoints of a built `path`: the
union of segment and spinning-run boundaries. The propagator never steps
across one of these.
"""
function breakpoints(b::SubpathBuilt)
    return normalize_breakpoints(vcat(path_segment_breakpoints(b),
                                      path_spinning_breakpoints(b)))
end

"""
    sample(path, s_values) -> Vector{NamedTuple}

Return the [`frame`](@ref) of a built `path` at each arc length in `s_values`.
"""
function sample(b::SubpathBuilt, s_values)
    return [frame(b, s) for s in s_values]
end

"""
    sample_uniform(path; n=256) -> Vector{NamedTuple}

Return `n` frames of a built `path` uniformly spaced over its arc length.
"""
function sample_uniform(b::SubpathBuilt; n::Int = 256)
    ss = range(0.0, Float64(_qc_nominalize(arc_length(b))); length = n)
    return sample(b, ss)
end

# -----------------------------------------------------------------------
# PathBuilt query interface — single _find_subpath glue + @eval forwards
# -----------------------------------------------------------------------

"""
    s_offsets(p::PathBuilt) -> Vector{Float64}

Return the cumulative global arc-length offset at the start of each subpath,
computed on demand.
"""
function s_offsets(p::PathBuilt)
    n = length(p.subpaths)
    offs = Vector{Float64}(undef, n)
    cum = 0.0
    @inbounds for i in 1:n
        offs[i] = cum
        cum += Float64(_qc_nominalize(arc_length(p.subpaths[i])))
    end
    return offs
end

"""
    s_end(p::PathBuilt) -> Float64

Return a `PathBuilt`'s total arc length (the sum of its subpath lengths).
"""
s_end(p::PathBuilt) = sum(Float64(_qc_nominalize(arc_length(b))) for b in p.subpaths;
                          init = 0.0)

arc_length(p::PathBuilt) = s_end(p)

function arc_length(::PathBuilt, s1, s2)
    @assert s2 >= s1 "arc_length: require s2 >= s1"
    return s2 - s1
end

path_length(p::PathBuilt) = arc_length(p)

"""
    _find_subpath(p::PathBuilt, s) -> (SubpathBuilt, s_local)

Locate the subpath of `p` containing global arc length `s`, returning it with
the subpath-local arc length.
"""
function _find_subpath(p::PathBuilt, s)
    n = length(p.subpaths)
    n == 0 && error("PathBuilt is empty")
    offs = s_offsets(p)
    @inbounds for i in 1:n
        local_end = offs[i] + Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        if s <= local_end + 1e-12 || i == n
            return p.subpaths[i], s - offs[i]
        end
    end
    error("PathBuilt: s = $s out of path bounds")
end

# Forward all "point-query at s" methods through _find_subpath.
for f in (:curvature, :geometric_torsion, :spinning_rate,
          :position, :tangent, :normal, :binormal, :frame)
    @eval function $f(p::PathBuilt, s::Real)
        sb, s_local = _find_subpath(p, s)
        return $f(sb, s_local)
    end
end

# Endpoint access
start_point(p::PathBuilt)   = position(p, 0.0)
end_point(p::PathBuilt)     = position(p, s_end(p))
start_tangent(p::PathBuilt) = tangent(p, 0.0)
end_tangent(p::PathBuilt)   = tangent(p, s_end(p))

function cartesian_distance(p::PathBuilt, s1::Real, s2::Real)
    return norm(position(p, s2) - position(p, s1))
end

function bounding_box(p::PathBuilt; n::Int = 512)
    ss = range(0.0, s_end(p); length = n)
    pts = [position(p, s) for s in ss]
    lo = minimum(reduce(hcat, pts); dims = 2) |> vec
    hi = maximum(reduce(hcat, pts); dims = 2) |> vec
    return (; lo, hi)
end

# Aggregate measures: sum across subpaths.
total_turning_angle(p::PathBuilt) = sum(total_turning_angle(b) for b in p.subpaths;
                                        init = 0.0)
total_torsion(p::PathBuilt)       = sum(total_torsion(b)       for b in p.subpaths;
                                        init = 0.0)

function total_spinning(p::PathBuilt;
                              s_start::Real = 0.0,
                              s_end::Real   = s_end(p),
                              rtol::Real    = 1e-8,
                              atol::Real    = 0.0)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    s_lo == s_hi && return 0.0
    total = 0.0
    offs = s_offsets(p)
    for i in eachindex(p.subpaths)
        L = Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        a_local = max(0.0, s_lo - offs[i])
        b_local = min(L,   s_hi - offs[i])
        b_local <= a_local && continue
        total += total_spinning(p.subpaths[i];
                                      s_start = a_local, s_end = b_local,
                                      rtol = rtol, atol = atol)
    end
    return total
end

function total_frame_rotation(p::PathBuilt;
                              s_start::Real = 0.0,
                              s_end::Real   = s_end(p),
                              rtol::Real    = 1e-8,
                              atol::Real    = 0.0)
    s_lo = Float64(_qc_nominalize(s_start))
    s_hi = Float64(_qc_nominalize(s_end))
    s_lo == s_hi && return 0.0
    total = 0.0
    offs = s_offsets(p)
    for i in eachindex(p.subpaths)
        L = Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        a_local = max(0.0, s_lo - offs[i])
        b_local = min(L,   s_hi - offs[i])
        b_local <= a_local && continue
        total += total_frame_rotation(p.subpaths[i];
                                      s_start = a_local, s_end = b_local,
                                      rtol = rtol, atol = atol)
    end
    return total
end

function breakpoints(p::PathBuilt)
    isempty(p.subpaths) && return Float64[]
    offs = s_offsets(p)
    bps = Float64[]
    for i in eachindex(p.subpaths)
        for x in breakpoints(p.subpaths[i])
            push!(bps, x + offs[i])
        end
    end
    return normalize_breakpoints(bps)
end

function sample(p::PathBuilt, s_values)
    return [frame(p, s) for s in s_values]
end

function sample_uniform(p::PathBuilt; n::Int = 256)
    ss = range(0.0, s_end(p); length = n)
    return sample(p, ss)
end

"""
    sample_path(p::PathBuilt, s1, s2; fidelity = 1.0) → PathSample

Adaptive sampling of a `PathBuilt` over `[s1, s2]`. Walks each constituent
`SubpathBuilt`'s clipped interval and concatenates samples with adjacent-
duplicate suppression at Subpath boundaries. Sample positions/tangents/etc.
are taken via `frame(p, s)` so they live in the global frame.
"""
function sample_path(p::PathBuilt, s1::Real, s2::Real; fidelity::Float64 = 1.0)
    @assert s2 > s1 "sample_path: require s2 > s1"
    @assert fidelity > 0.0 "sample_path: fidelity must be positive"

    s_lo = Float64(_qc_nominalize(s1))
    s_hi = Float64(_qc_nominalize(s2))

    offs = s_offsets(p)
    all_s = Float64[]
    for i in eachindex(p.subpaths)
        L = Float64(_qc_nominalize(arc_length(p.subpaths[i])))
        a = max(s_lo, offs[i])
        b = min(s_hi, offs[i] + L)
        b <= a && continue
        ps = sample_path(p.subpaths[i], a - offs[i], b - offs[i]; fidelity = fidelity)
        for sm in ps.samples
            sg = sm.s + offs[i]
            if isempty(all_s) || !(sg ≈ all_s[end])
                push!(all_s, sg)
            end
        end
    end

    if isempty(all_s)
        all_s = [s_lo, s_hi]
    elseif length(all_s) == 1
        push!(all_s, s_hi)
    end

    n = length(all_s)
    samples = Vector{Sample}(undef, n)
    for i in eachindex(all_s)
        fr = frame(p, all_s[i])
        samples[i] = Sample(
            all_s[i],
            fr.position,
            fr.tangent,
            fr.normal,
            fr.binormal,
            fr.curvature,
            fr.geometric_torsion,
            fr.spinning_rate,
        )
    end
    return PathSample(samples, s_lo, s_hi, n)
end
