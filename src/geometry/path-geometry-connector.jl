"""
path-geometry-connector.jl

QuinticConnector{T<:Real}: parametric quintic Hermite curve used as the resolved
form of a `JumpBy` segment and of the terminal `jumpto!` seal.

Matches position + tangent direction + curvature vector at both endpoints (G2).
A single positive scalar λ (handle scale) supplies the parameter-speed degree of
freedom that geometric G2 boundary data does not pin down.

The struct is parametric over `T<:Real` so MCM `Particles` flow through.
Branching predicates (`κ ≤ 1/R_min`, `s_table` bisection, Newton refinement)
are evaluated on nominalized Float64 scalars; the final coefficients are built
from the original `T`-typed inputs so uncertainty propagates.

This file is `include`d into the same lexical scope as `path-geometry.jl`. It
uses `AbstractPathSegment` and `AbstractMeta` from there.
"""

# ---------------------------------------------------------------------------
# Nominalization (deterministic branching under MCM)
# ---------------------------------------------------------------------------

# `Particles{T,N} <: Real` (but not <: AbstractFloat). Dispatch the identity
# branch on the concrete already-nominal types so Particles falls through to
# the pmean lookup below. Listing AbstractFloat first covers Float64/BigFloat;
# Integer covers Int*; the generic Real method catches Particles and other
# custom Real subtypes that need pmean reduction.
_qc_nominalize(x::AbstractFloat) = x
_qc_nominalize(x::Integer) = x
function _qc_nominalize(x::Real)
    if isdefined(Main, :MonteCarloMeasurements)
        M = getfield(Main, :MonteCarloMeasurements)
        if isdefined(M, :pmean)
            f = getfield(M, :pmean)
            applicable(f, x) && return f(x)
        end
    end
    return x
end
# Fallback for anything that isn't a Real (shouldn't normally happen but
# preserves the prior catch-all behavior).
_qc_nominalize(x) = x

# ---------------------------------------------------------------------------
# Gauss–Legendre 4-point nodes/weights for arc-length quadrature
# ---------------------------------------------------------------------------

const _QC_GAUSS4_NODES   = (-0.8611363115940526, -0.3399810435848563,
                              0.3399810435848563,  0.8611363115940526)
const _QC_GAUSS4_WEIGHTS = ( 0.3478548451374538,  0.6521451548625461,
                              0.6521451548625461,  0.3478548451374538)

# ---------------------------------------------------------------------------
# Struct
# ---------------------------------------------------------------------------

"""
    QuinticConnector{T<:Real} <: AbstractPathSegment

Resolved form of a [`JumpBy`](@ref) segment and of the terminal `jumpto!` seal: a
quintic Hermite curve

```
r(u) = a₀ + a₁ u + a₂ u² + a₃ u³ + a₄ u⁴ + a₅ u⁵,   u ∈ [0,1]
```

in 3D that matches position, tangent direction, and curvature vector at both endpoints
(G2 continuity). Implements the [`AbstractPathSegment`](@ref) local-geometry interface.

# Fields

- `a`: `6 × 3` coefficient matrix; row `i` stores the coefficient of `uⁱ⁻¹`.
- `lambda`: chosen handle scale (`Float64` — must be deterministic under MCM inputs).
- `s_table`: arc-length lookup, where `s_table[i]` is the arc length from `u = 0` to
  `u = (i-1)/(n-1)`.
- `e1_table`: `n × 3` parallel-transport lookup; row `i` holds the local
  components of the transported `e1` at `u = (i-1)/(n-1)`, with
  `e1_table[1, :] = x̂` (the entry transverse axis). Built by discrete
  double-reflection transport along the same `u`-grid as `s_table`.
- `twist`: mechanical-twist rate (rad/m) carried from the originating
  `JumpBy`/seal (see [`TwistRate`](@ref)).
- `meta`: per-segment annotation bag (see [`AbstractMeta`](@ref)).
"""
struct QuinticConnector{T<:Real} <: AbstractPathSegment
    a        :: Matrix{T}    # 6 × 3
    lambda   :: Float64
    s_table  :: Vector{T}
    e1_table :: Matrix{T}    # n × 3 transported e1 at the s_table u-nodes
    twist    :: TwistRate    # mechanical-twist rate (rad/m); see `TwistRate`
    meta     :: Vector{AbstractMeta}
end

"""
    QuinticConnector(a, lambda, s_table; twist = nothing, meta = AbstractMeta[])

Construct a [`QuinticConnector`](@ref) from a precomputed coefficient matrix `a`, handle
scale `lambda`, and arc-length table `s_table`. The parallel-transport `e1_table` is
computed here by discrete double-reflection over the `s_table` `u`-grid.

Coerces `lambda` to `Float64` and wraps `meta` in an `AbstractMeta` vector. `twist`
is the connector's mechanical-twist rate (see [`TwistRate`](@ref)). Callers
normally obtain `a`/`s_table` from the internal connector builders rather than calling
this directly.
"""
QuinticConnector(a::Matrix{T}, lambda::Real, s_table::AbstractVector{T};
                 twist = nothing, meta = AbstractMeta[]) where {T<:Real} =
    QuinticConnector{T}(a, Float64(lambda), Vector{T}(s_table),
                        _qc_build_e1_table(a, length(s_table)),
                        twist, Vector{AbstractMeta}(meta))

# ---------------------------------------------------------------------------
# Discrete parallel transport (double-reflection / rotation-minimizing frame)
# ---------------------------------------------------------------------------

"""
    _qc_tangent_at_u(a, u) -> Vector

Unit tangent of the quintic at parameter `u` (not arc length). Zero-speed
points (degenerate/zero-length connectors) map to the local `ẑ`, mirroring
`tangent_local`.
"""
function _qc_tangent_at_u(a::AbstractMatrix, u::Real)
    d1 = _qc_eval_d1(a, u)
    spd2 = d1[1]^2 + d1[2]^2 + d1[3]^2
    if Float64(_qc_nominalize(spd2)) < 1e-30
        T = eltype(a)
        return [zero(T), zero(T), one(T)]
    end
    return d1 ./ sqrt(spd2)
end

"""
    _qc_double_reflect(v, x0, t0, x1, t1) -> Vector

One double-reflection parallel-transport step (Wang–Jüttler–Zheng–Liu
rotation-minimizing frame, doi:10.1145/1330511.1330513): transport the
transverse vector `v` from the point/tangent pair `(x0, t0)` to `(x1, t1)`.
Degenerate steps (coincident points, parallel reflection data) skip the
affected reflection on a nominalized predicate; the algebra itself is built
from the inputs' element type so MCM `Particles` propagate.
"""
function _qc_double_reflect(v::AbstractVector, x0::AbstractVector, t0::AbstractVector,
                            x1::AbstractVector, t1::AbstractVector)
    d = x1 .- x0
    c1 = d[1]^2 + d[2]^2 + d[3]^2
    Float64(_qc_nominalize(c1)) < 1e-30 && return v
    vL = v .- (2 * (d[1]*v[1] + d[2]*v[2] + d[3]*v[3]) / c1) .* d
    tL = t0 .- (2 * (d[1]*t0[1] + d[2]*t0[2] + d[3]*t0[3]) / c1) .* d
    e = t1 .- tL
    c2 = e[1]^2 + e[2]^2 + e[3]^2
    Float64(_qc_nominalize(c2)) < 1e-30 && return vL
    return vL .- (2 * (e[1]*vL[1] + e[2]*vL[2] + e[3]*vL[3]) / c2) .* e
end

"""
    _qc_build_e1_table(a, n) -> Matrix

Build the `n × 3` transported-`e1` lookup for a quintic with coefficients `a`:
row `i` is the parallel transport of the entry `x̂` to `u = (i-1)/(n-1)`,
accumulated by double-reflection steps between consecutive nodes. Exact for
straight and zero-length connectors (identity transport); O(h²)-accurate in
the node spacing otherwise.
"""
function _qc_build_e1_table(a::Matrix{T}, n::Int) where {T<:Real}
    e1 = Matrix{T}(undef, n, 3)
    v = T[one(T), zero(T), zero(T)]
    e1[1, 1] = v[1]; e1[1, 2] = v[2]; e1[1, 3] = v[3]
    x_prev = _qc_eval(a, 0.0)
    t_prev = _qc_tangent_at_u(a, 0.0)
    @inbounds for i in 2:n
        u = (i - 1) / (n - 1)
        x_cur = _qc_eval(a, u)
        t_cur = _qc_tangent_at_u(a, u)
        v = _qc_double_reflect(v, x_prev, t_prev, x_cur, t_cur)
        e1[i, 1] = v[1]; e1[i, 2] = v[2]; e1[i, 3] = v[3]
        x_prev = x_cur
        t_prev = t_cur
    end
    return e1
end

# ---------------------------------------------------------------------------
# Coefficient solve (closed-form quintic Hermite)
# ---------------------------------------------------------------------------

function _qc_coefficients(P0::AbstractVector{T}, V0::AbstractVector{T}, A0::AbstractVector{T},
                          P1::AbstractVector{T}, V1::AbstractVector{T}, A1::AbstractVector{T}) where {T<:Real}
    half = T(1)/T(2)
    a0 = collect(P0)
    a1 = collect(V0)
    a2 = A0 .* half
    b0 = P1 .- (a0 .+ a1 .+ a2)
    b1 = V1 .- (a1 .+ T(2) .* a2)
    b2 = A1 .- (T(2) .* a2)
    a3 =  T(10) .* b0 .- T(4) .* b1 .+ half .* b2
    a4 = -T(15) .* b0 .+ T(7) .* b1 .- b2
    a5 =  T(6)  .* b0 .- T(3) .* b1 .+ half .* b2
    coeffs = Matrix{T}(undef, 6, 3)
    @inbounds for j in 1:3
        coeffs[1,j] = a0[j]
        coeffs[2,j] = a1[j]
        coeffs[3,j] = a2[j]
        coeffs[4,j] = a3[j]
        coeffs[5,j] = a4[j]
        coeffs[6,j] = a5[j]
    end
    return coeffs
end

# ---------------------------------------------------------------------------
# Polynomial evaluators (Horner)
# ---------------------------------------------------------------------------

@inline function _qc_eval(a::AbstractMatrix, u::Real)
    val = collect(a[6, :])
    @inbounds for i in 5:-1:1
        val = val .* u .+ collect(a[i, :])
    end
    return val
end

@inline function _qc_eval_d1(a::AbstractMatrix, u::Real)
    c1 = collect(a[2, :])
    c2 = 2  .* collect(a[3, :])
    c3 = 3  .* collect(a[4, :])
    c4 = 4  .* collect(a[5, :])
    c5 = 5  .* collect(a[6, :])
    val = c5
    @inbounds for c in (c4, c3, c2, c1)
        val = val .* u .+ c
    end
    return val
end

@inline function _qc_eval_d2(a::AbstractMatrix, u::Real)
    c0 = 2  .* collect(a[3, :])
    c1 = 6  .* collect(a[4, :])
    c2 = 12 .* collect(a[5, :])
    c3 = 20 .* collect(a[6, :])
    val = c3
    @inbounds for c in (c2, c1, c0)
        val = val .* u .+ c
    end
    return val
end

@inline function _qc_eval_d3(a::AbstractMatrix, u::Real)
    c0 = 6  .* collect(a[4, :])
    c1 = 24 .* collect(a[5, :])
    c2 = 60 .* collect(a[6, :])
    val = c2
    @inbounds for c in (c1, c0)
        val = val .* u .+ c
    end
    return val
end

# ---------------------------------------------------------------------------
# Speed and arc-length table
# ---------------------------------------------------------------------------

function _qc_speed(a::AbstractMatrix, u::Real)
    d1 = _qc_eval_d1(a, u)
    return sqrt(d1[1]^2 + d1[2]^2 + d1[3]^2)
end

function _qc_build_table(a::Matrix{T}, n::Int) where {T<:Real}
    s = Vector{T}(undef, n)
    s[1] = zero(T)
    @inbounds for i in 2:n
        u0 = (i - 2) / (n - 1)
        u1 = (i - 1) / (n - 1)
        um = (u0 + u1) / 2
        uh = (u1 - u0) / 2
        s[i] = s[i-1] + uh * (
            _QC_GAUSS4_WEIGHTS[1] * _qc_speed(a, um + uh * _QC_GAUSS4_NODES[1]) +
            _QC_GAUSS4_WEIGHTS[2] * _qc_speed(a, um + uh * _QC_GAUSS4_NODES[2]) +
            _QC_GAUSS4_WEIGHTS[3] * _qc_speed(a, um + uh * _QC_GAUSS4_NODES[3]) +
            _QC_GAUSS4_WEIGHTS[4] * _qc_speed(a, um + uh * _QC_GAUSS4_NODES[4]))
    end
    return s
end

# ---------------------------------------------------------------------------
# s ↔ u inversion (deterministic under MCM)
# ---------------------------------------------------------------------------

function _qc_t_from_s(seg::QuinticConnector{T}, s_target::Real) where {T<:Real}
    s_table = seg.s_table
    n = length(s_table)
    L_nom = Float64(_qc_nominalize(s_table[end]))
    s_nom = Float64(_qc_nominalize(s_target))
    L_nom < 1e-15 && return 0.0
    sc = clamp(s_nom, 0.0, L_nom)

    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        s_mid = Float64(_qc_nominalize(s_table[mid]))
        s_mid <= sc ? (lo = mid) : (hi = mid)
    end
    u0 = (lo - 1) / (n - 1)
    u1 = (hi - 1) / (n - 1)
    s0 = Float64(_qc_nominalize(s_table[lo]))
    ds = Float64(_qc_nominalize(s_table[hi])) - s0

    u = ds > 1e-15 ? u0 + (sc - s0) / ds * (u1 - u0) : u0

    a_nom = [Float64(_qc_nominalize(x)) for x in seg.a]
    a_nom_mat = reshape(a_nom, size(seg.a))
    for _ in 1:2
        spd_mid = _qc_speed(a_nom_mat, (u0 + u) / 2)
        spd_mid < 1e-15 && break
        s_est = s0 + (u - u0) * spd_mid
        spd   = _qc_speed(a_nom_mat, u)
        spd   < 1e-15 && break
        u = clamp(u - (s_est - sc) / spd, u0, u1)
    end
    return u
end

# ---------------------------------------------------------------------------
# AbstractPathSegment interface
# ---------------------------------------------------------------------------

arc_length(seg::QuinticConnector) = seg.s_table[end]

function curvature(seg::QuinticConnector, s::Real)
    u  = _qc_t_from_s(seg, s)
    d1 = _qc_eval_d1(seg.a, u)
    d2 = _qc_eval_d2(seg.a, u)
    spd2 = d1[1]^2 + d1[2]^2 + d1[3]^2
    Float64(_qc_nominalize(spd2)) < 1e-30 && return zero(spd2)
    cx = d1[2]*d2[3] - d1[3]*d2[2]
    cy = d1[3]*d2[1] - d1[1]*d2[3]
    cz = d1[1]*d2[2] - d1[2]*d2[1]
    spd = sqrt(spd2)
    return sqrt(cx^2 + cy^2 + cz^2) / (spd * spd2)
end

function geometric_torsion(seg::QuinticConnector, s::Real)
    u  = _qc_t_from_s(seg, s)
    d1 = _qc_eval_d1(seg.a, u)
    d2 = _qc_eval_d2(seg.a, u)
    d3 = _qc_eval_d3(seg.a, u)
    cx = d1[2]*d2[3] - d1[3]*d2[2]
    cy = d1[3]*d2[1] - d1[1]*d2[3]
    cz = d1[1]*d2[2] - d1[2]*d2[1]
    denom = cx^2 + cy^2 + cz^2
    Float64(_qc_nominalize(denom)) < 1e-30 && return zero(denom)
    return (cx*d3[1] + cy*d3[2] + cz*d3[3]) / denom
end

position_local(seg::QuinticConnector, s::Real) =
    _qc_eval(seg.a, _qc_t_from_s(seg, s))

function tangent_local(seg::QuinticConnector, s::Real)
    u  = _qc_t_from_s(seg, s)
    d1 = _qc_eval_d1(seg.a, u)
    spd2 = d1[1]^2 + d1[2]^2 + d1[3]^2
    if Float64(_qc_nominalize(spd2)) < 1e-30
        T = eltype(seg.a)
        return [zero(T), zero(T), one(T)]
    end
    return d1 ./ sqrt(spd2)
end

function _qc_normal_from_derivs(d1::AbstractVector, d2::AbstractVector)
    spd2 = d1[1]^2 + d1[2]^2 + d1[3]^2
    Tc = promote_type(eltype(d1), eltype(d2))
    if Float64(_qc_nominalize(spd2)) < 1e-30
        Tv = Tc[zero(Tc), zero(Tc), one(Tc)]
    else
        Tv = d1 ./ sqrt(spd2)
    end
    dotTd2 = Tv[1]*d2[1] + Tv[2]*d2[2] + Tv[3]*d2[3]
    acc = d2 .- dotTd2 .* Tv
    an2 = acc[1]^2 + acc[2]^2 + acc[3]^2
    if Float64(_qc_nominalize(an2)) < 1e-30
        return Tv, _qc_perp_unit(Tv)
    end
    return Tv, acc ./ sqrt(an2)
end

function _normal_local(seg::QuinticConnector, s::Real)
    u  = _qc_t_from_s(seg, s)
    d1 = _qc_eval_d1(seg.a, u)
    d2 = _qc_eval_d2(seg.a, u)
    _, N = _qc_normal_from_derivs(d1, d2)
    return N
end

function _binormal_local(seg::QuinticConnector, s::Real)
    u  = _qc_t_from_s(seg, s)
    d1 = _qc_eval_d1(seg.a, u)
    d2 = _qc_eval_d2(seg.a, u)
    Tv, Nv = _qc_normal_from_derivs(d1, d2)
    return [Tv[2]*Nv[3] - Tv[3]*Nv[2],
            Tv[3]*Nv[1] - Tv[1]*Nv[3],
            Tv[1]*Nv[2] - Tv[2]*Nv[1]]
end

end_position_local(seg::QuinticConnector) = _qc_eval(seg.a, 1.0)

function _end_frame_local(seg::QuinticConnector)
    d1 = _qc_eval_d1(seg.a, 1.0)
    d2 = _qc_eval_d2(seg.a, 1.0)
    Tv, Nv = _qc_normal_from_derivs(d1, d2)
    Bv = [Tv[2]*Nv[3] - Tv[3]*Nv[2],
          Tv[3]*Nv[1] - Tv[1]*Nv[3],
          Tv[1]*Nv[2] - Tv[2]*Nv[1]]
    return (Tv, Nv, Bv)
end

# Curvature vector k⃗ = dT̂/ds from the raw parametric derivatives:
# k⃗ = (r″ − T̂(T̂·r″)) / |r′|². No κ normalization, so no inflection or
# near-straight degeneracy — zero-speed points return the zero vector.
function _curvature_vector_local(seg::QuinticConnector, s::Real)
    u  = _qc_t_from_s(seg, s)
    d1 = _qc_eval_d1(seg.a, u)
    d2 = _qc_eval_d2(seg.a, u)
    spd2 = d1[1]^2 + d1[2]^2 + d1[3]^2
    if Float64(_qc_nominalize(spd2)) < 1e-30
        return [zero(spd2), zero(spd2), zero(spd2)]
    end
    proj = (d1[1]*d2[1] + d1[2]*d2[2] + d1[3]*d2[3]) / spd2
    acc = d2 .- proj .* d1
    return acc ./ spd2
end

# Parallel transport on the connector: decompose `v` on the entry transverse
# pair (x̂, ŷ), look up the transported pair at the table node at or below the
# query point, and finish with one double-reflection step from the node to the
# exact query point — continuous in `s`, O(h²) in the table spacing.
function _parallel_transport_local(seg::QuinticConnector, v::AbstractVector, s)
    n = length(seg.s_table)
    u = _qc_t_from_s(seg, s)          # deterministic Float64 under MCM
    lo = clamp(floor(Int, u * (n - 1) + 1e-12) + 1, 1, n)
    u_lo = (lo - 1) / (n - 1)
    e1_lo = [seg.e1_table[lo, 1], seg.e1_table[lo, 2], seg.e1_table[lo, 3]]
    t_lo = _qc_tangent_at_u(seg.a, u_lo)
    e2_lo = [t_lo[2]*e1_lo[3] - t_lo[3]*e1_lo[2],
             t_lo[3]*e1_lo[1] - t_lo[1]*e1_lo[3],
             t_lo[1]*e1_lo[2] - t_lo[2]*e1_lo[1]]
    w = v[1] .* e1_lo .+ v[2] .* e2_lo
    x_lo = _qc_eval(seg.a, u_lo)
    x_u = _qc_eval(seg.a, u)
    t_u = _qc_tangent_at_u(seg.a, u)
    return _qc_double_reflect(w, x_lo, t_lo, x_u, t_u)
end

# Unit vector perpendicular to t (arbitrary orientation), parametric in T.
function _qc_perp_unit(t::AbstractVector)
    T = eltype(t)
    t3_nom = abs(Float64(_qc_nominalize(t[3])))
    ref = t3_nom < 0.9 ? T[zero(T), zero(T), one(T)] :
                         T[one(T), zero(T), zero(T)]
    dotrt = ref[1]*t[1] + ref[2]*t[2] + ref[3]*t[3]
    n = ref .- dotrt .* t
    nn2 = n[1]^2 + n[2]^2 + n[3]^2
    if Float64(_qc_nominalize(nn2)) < 1e-24
        return T[one(T), zero(T), zero(T)]
    end
    return n ./ sqrt(nn2)
end

# ---------------------------------------------------------------------------
# Peak curvature (sampled grid; nominalized)
# ---------------------------------------------------------------------------

function _qc_peak_curvature(a::AbstractMatrix; n_check::Int = 128)
    κ_max = 0.0
    @inbounds for i in 0:n_check
        u  = i / n_check
        d1 = _qc_eval_d1(a, u)
        d2 = _qc_eval_d2(a, u)
        spd2 = d1[1]^2 + d1[2]^2 + d1[3]^2
        spd2_nom = Float64(_qc_nominalize(spd2))
        spd2_nom < 1e-30 && continue
        cx = d1[2]*d2[3] - d1[3]*d2[2]
        cy = d1[3]*d2[1] - d1[1]*d2[3]
        cz = d1[1]*d2[2] - d1[2]*d2[1]
        cn2_nom = Float64(_qc_nominalize(cx^2 + cy^2 + cz^2))
        κ = sqrt(cn2_nom) / (sqrt(spd2_nom) * spd2_nom)
        κ > κ_max && (κ_max = κ)
    end
    return κ_max
end

# ---------------------------------------------------------------------------
# Coefficient assembly at a chosen λ
# ---------------------------------------------------------------------------

function _qc_assemble(p1_local::AbstractVector, t_hat_in::AbstractVector,
                       t_hat_out::AbstractVector, K0::AbstractVector,
                       K1::AbstractVector, λ::Real)
    Tc = promote_type(eltype(p1_local), eltype(t_hat_in), eltype(t_hat_out),
                      eltype(K0), eltype(K1), typeof(λ))
    P0 = Tc[zero(Tc), zero(Tc), zero(Tc)]
    P1 = Tc.(p1_local)
    λT = Tc(λ)
    V0 = λT .* Tc.(t_hat_in)
    V1 = λT .* Tc.(t_hat_out)
    A0 = (λT*λT) .* Tc.(K0)
    A1 = (λT*λT) .* Tc.(K1)
    return _qc_coefficients(P0, V0, A0, P1, V1, A1)
end

# Initial λ guess: chord-driven floor with curvature-aware lower bound.
function _qc_initial_lambda(p1_local::AbstractVector, t_hat_in::AbstractVector,
                             t_hat_out::AbstractVector, R_min::Float64)
    chord_nom = sqrt(Float64(_qc_nominalize(
        p1_local[1]^2 + p1_local[2]^2 + p1_local[3]^2)))
    dot_tt = t_hat_in[1]*t_hat_out[1] + t_hat_in[2]*t_hat_out[2] + t_hat_in[3]*t_hat_out[3]
    dθ = acos(clamp(Float64(_qc_nominalize(dot_tt)), -1.0, 1.0))
    return max(chord_nom, R_min * dθ, 1e-6)
end

# Project K onto the plane perpendicular to t (silently — the new segment's
# local frame is constructed so the prior endpoint's K should already be ⊥ to
# ẑ; user-supplied K may have a tangential component which we simply discard).
function _qc_project_perp(K::AbstractVector, t::AbstractVector)
    dot_tk = t[1]*K[1] + t[2]*K[2] + t[3]*K[3]
    return K .- dot_tk .* t
end

# ---------------------------------------------------------------------------
# Straight / zero-length terminal connector (natural seal)
# ---------------------------------------------------------------------------

"""
    _build_straight_connector(extra, T; meta = AbstractMeta[]) → QuinticConnector{T}

Build a terminal connector that travels straight along the local tangent ẑ for
`extra` meters, with no bending. Used by `seal!` (natural seal) to terminate a
Subpath at its natural exit without invoking the quintic solver, whose
coincident-endpoint solve is ill-conditioned.

`extra == 0` yields a true zero-length connector (all coefficients zero); its
query path is degenerate-safe in `QuinticConnector` (zero speed maps to the
local ẑ tangent). `extra > 0` yields an exact straight line `r(u) = (0,0,extra·u)`
with a linear arc-length table.

`meta` is the seal's annotation bag; it is stored on the connector so any
consuming-layer annotation (e.g. a thermal `:T_K`) is carried through like any
other segment's meta.
"""
function _build_straight_connector(extra::Real, ::Type{T};
                                   twist = nothing,
                                   meta::AbstractVector{<:AbstractMeta} = AbstractMeta[]) where {T<:Real}
    L = T(extra)
    a = zeros(T, 6, 3)
    a[2, 3] = L                       # r(u) = (0, 0, L·u): linear in ẑ
    s_table = T[zero(T), L]           # straight ⇒ two table points suffice
    return QuinticConnector(a, 1.0, s_table; twist = twist, meta = meta)
end

# ---------------------------------------------------------------------------
# Build a connector with optional min_bend_radius (λ search)
# ---------------------------------------------------------------------------

"""
    _build_quintic_connector(p1_local, t_hat_out, K0_local, K1_local;
                             min_bend_radius, target_path_length,
                             n_table, meta) → QuinticConnector

Construct a `QuinticConnector` from a Subpath terminal-jumpto resolve or a
`JumpBy` resolve.

`t_hat_in` is the local-frame ẑ = (0,0,1). `K0_local` is the prior segment's
terminal curvature vector projected into the new segment's local frame.
`K1_local` is the user-specified outgoing curvature (zero by default). Both
are silently re-projected onto the plane perpendicular to their tangent.

If `min_bend_radius` is set, performs an exponential-bracket / bisection
search over the handle scale λ until the sampled peak curvature falls below
`1/min_bend_radius`. 

If `target_path_length` is set, instead searches λ
until the connector's arc length matches the target — supplied via
`build(...; jumpto_target_length=…)` (a consuming layer, the fiber, uses it to
thermally expand the terminal connector) while the endpoint stays a
lab-frame invariant.

The two constraints can be
combined: when both are set, λ is chosen by the arc-length search and the
result is validated against `min_bend_radius`. All branching uses
nominalized scalars so the chosen λ is deterministic under MCM inputs; the
final coefficients carry the original `T` element type.
"""
function _build_quintic_connector(p1_local::AbstractVector,
                                  t_hat_out::AbstractVector,
                                  K0_local::AbstractVector,
                                  K1_local::AbstractVector;
                                  min_bend_radius::Union{Nothing, Real} = nothing,
                                  target_path_length::Union{Nothing, Real} = nothing,
                                  n_table::Int = 256,
                                  n_check::Int = 128,
                                  growth::Float64 = 1.5,
                                  max_iter::Int = 30,
                                  bisect_iter::Int = 64,
                                  rel_tol::Float64 = 1e-6,
                                  curvature_tol::Float64 = 1e-8,
                                  twist = nothing,
                                  meta::AbstractVector{<:AbstractMeta} = AbstractMeta[])
    Tc = promote_type(eltype(p1_local), eltype(t_hat_out),
                      eltype(K0_local), eltype(K1_local), Float64)
    t_hat_in = Tc[zero(Tc), zero(Tc), one(Tc)]

    # Project K's onto plane ⊥ tangent. K0 should be ⊥ ẑ by construction;
    # K1 from user may not be exactly ⊥ to t_hat_out.
    K0_perp = _qc_project_perp(Tc.(K0_local), t_hat_in)
    K1_perp = _qc_project_perp(Tc.(K1_local), Tc.(t_hat_out))

    R_min = isnothing(min_bend_radius) ? Inf : Float64(min_bend_radius)
    R_min > 0.0 || throw(ArgumentError("min_bend_radius must be positive"))
    κ_limit = isfinite(R_min) ? 1.0 / R_min : Inf

    if isfinite(κ_limit)
        K0n = sqrt(Float64(_qc_nominalize(K0_perp[1]^2 + K0_perp[2]^2 + K0_perp[3]^2)))
        K1n = sqrt(Float64(_qc_nominalize(K1_perp[1]^2 + K1_perp[2]^2 + K1_perp[3]^2)))
        K0n > κ_limit + curvature_tol && throw(ArgumentError(
            "connector: incoming endpoint curvature ($(round(K0n;digits=3)) m⁻¹) " *
            "exceeds 1/min_bend_radius ($(round(κ_limit;digits=3)) m⁻¹)"))
        K1n > κ_limit + curvature_tol && throw(ArgumentError(
            "connector: outgoing endpoint curvature ($(round(K1n;digits=3)) m⁻¹) " *
            "exceeds 1/min_bend_radius ($(round(κ_limit;digits=3)) m⁻¹)"))
    end

    R_for_init = isfinite(R_min) ? R_min : 0.0
    λ = _qc_initial_lambda(p1_local, t_hat_in, t_hat_out, R_for_init)

    # Length-constrained mode: search λ so arc length ≈ target. λ is monotone
    # in arc length over the typical range (chord-aligned baseline upward). If
    # both target and min_bend_radius are set, pick λ by arc length and verify
    # peak curvature post-hoc.
    if !isnothing(target_path_length)
        target = Float64(_qc_nominalize(target_path_length))
        target > 0.0 || throw(ArgumentError(
            "target_path_length must be positive; got $(target)"))
        chord_nom = sqrt(Float64(_qc_nominalize(
            p1_local[1]^2 + p1_local[2]^2 + p1_local[3]^2)))
        target < chord_nom * (1 - rel_tol) && throw(ArgumentError(
            "target path length infeasible: target=$(target) m is shorter " *
            "than chord ($(round(chord_nom;digits=6)) m)"))

        arc_at = function (λv)
            c = _qc_assemble(p1_local, t_hat_in, t_hat_out, K0_perp, K1_perp, λv)
            s = _qc_build_table(c, n_table)
            return Float64(_qc_nominalize(s[end]))
        end

        arc_init = arc_at(λ)

        local λ_lo::Float64, λ_hi::Float64
        if abs(arc_init - target) <= max(rel_tol * target, 1e-12)
            λ_lo = λ
            λ_hi = λ
        elseif arc_init < target
            # Grow λ until arc length ≥ target.
            λ_lo = λ
            λ_hi = λ
            found = false
            for _ in 1:max_iter
                λ_hi = λ_hi * growth
                if arc_at(λ_hi) >= target
                    found = true
                    break
                end
                λ_lo = λ_hi
            end
            found || throw(ArgumentError(
                "target path length infeasible: arc length did not reach " *
                "$(target) m within λ=$(round(λ_hi;digits=3))"))
        else
            # Shrink λ until arc length ≤ target.
            λ_lo = λ
            λ_hi = λ
            found = false
            for _ in 1:max_iter
                λ_lo = λ_lo / growth
                if arc_at(λ_lo) <= target
                    found = true
                    break
                end
                λ_hi = λ_lo
            end
            found || throw(ArgumentError(
                "target path length infeasible: arc length did not shrink to " *
                "$(target) m down to λ=$(round(λ_lo;digits=6))"))
        end

        # Bisect (arc length is monotone-ish in λ over the bracket).
        for _ in 1:bisect_iter
            (λ_hi - λ_lo) < rel_tol * max(λ_hi, 1e-12) && break
            λ_mid = (λ_lo + λ_hi) / 2
            if arc_at(λ_mid) < target
                λ_lo = λ_mid
            else
                λ_hi = λ_mid
            end
        end
        λ_final = (λ_lo + λ_hi) / 2

        coeffs_final = _qc_assemble(p1_local, t_hat_in, t_hat_out,
                                    K0_perp, K1_perp, λ_final)
        s_table = _qc_build_table(coeffs_final, n_table)

        if isfinite(κ_limit)
            κ_final = _qc_peak_curvature(coeffs_final; n_check = n_check)
            κ_final > κ_limit + curvature_tol && throw(ArgumentError(
                "target_path_length=$(target) m and min_bend_radius=$(R_min) m " *
                "are jointly infeasible: peak curvature " *
                "$(round(κ_final;digits=3)) m⁻¹ exceeds 1/R_min=" *
                "$(round(κ_limit;digits=3)) m⁻¹"))
        end

        return QuinticConnector(coeffs_final, λ_final, s_table; twist = twist, meta = meta)
    end

    # Unconstrained: build at the initial λ.
    if !isfinite(κ_limit)
        coeffs = _qc_assemble(p1_local, t_hat_in, t_hat_out, K0_perp, K1_perp, λ)
        s_table = _qc_build_table(coeffs, n_table)
        return QuinticConnector(coeffs, λ, s_table; twist = twist, meta = meta)
    end

    # Test λ₀.
    coeffs0 = _qc_assemble(p1_local, t_hat_in, t_hat_out, K0_perp, K1_perp, λ)
    κ0 = _qc_peak_curvature(coeffs0; n_check = n_check)
    if κ0 <= κ_limit + curvature_tol
        s_table = _qc_build_table(coeffs0, n_table)
        return QuinticConnector(coeffs0, λ, s_table; twist = twist, meta = meta)
    end

    # Exponential bracket: grow λ until feasible.
    λ_lo = λ
    λ_hi = λ
    found = false
    for _ in 1:max_iter
        λ_hi = λ_hi * growth
        coeffs_hi = _qc_assemble(p1_local, t_hat_in, t_hat_out, K0_perp, K1_perp, λ_hi)
        κ_hi = _qc_peak_curvature(coeffs_hi; n_check = n_check)
        if κ_hi <= κ_limit + curvature_tol
            found = true
            break
        end
        λ_lo = λ_hi
    end
    found || throw(ArgumentError(
        "connector: min_bend_radius=$(R_min) m infeasible; could not bring " *
        "peak curvature below $(round(κ_limit;digits=3)) m⁻¹ within λ=$(round(λ_hi;digits=3))."))

    # Bisect for tighter λ.
    for _ in 1:bisect_iter
        (λ_hi - λ_lo) < rel_tol * λ_hi && break
        λ_mid = (λ_lo + λ_hi) / 2
        coeffs_mid = _qc_assemble(p1_local, t_hat_in, t_hat_out, K0_perp, K1_perp, λ_mid)
        κ_mid = _qc_peak_curvature(coeffs_mid; n_check = n_check)
        if κ_mid <= κ_limit + curvature_tol
            λ_hi = λ_mid
        else
            λ_lo = λ_mid
        end
    end

    coeffs_final = _qc_assemble(p1_local, t_hat_in, t_hat_out, K0_perp, K1_perp, λ_hi)
    s_table = _qc_build_table(coeffs_final, n_table)
    return QuinticConnector(coeffs_final, λ_hi, s_table; twist = twist, meta = meta)
end
