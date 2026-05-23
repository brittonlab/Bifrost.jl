using LinearAlgebra
using Printf


# Code Overview:
# - fiber-path.jl defines the fiber model, source specification, and breakpoint assembly.
# - generator_K and generator_Kω assemble the local Jones dynamics.
# - the numerical propagation stack consists of
#     - exp_jones_generator
#     - exp_midpoint_step
#     - phase_insensitive_error
#     - propagate_interval!
#     - propagate_piecewise
# - the DGD sensitivity stack consists of
#     - exp_sensitivity_midpoint_step
#     - sensitivity_phase_insensitive_error
#     - propagate_interval_sensitivity!
#     - propagate_piecewise_sensitivity
#     - output_dgd
#
# fiber-path-plot.jl sits on top of that stack; it provides visual diagnostics.

if !isdefined(@__MODULE__, :Fiber)
    include(joinpath(@__DIR__, "fiber", "fiber-path.jl"))
end
# NOTE: fiber-path-plot.jl is intentionally not included here. Callers that need
# plotting helpers should include it themselves. This keeps path-integral.jl
# focused on the propagation stack and lets test files guard their own plot
# includes without double-inclusion.

# ----------------------------
# MCM scalar reduction
# ----------------------------
#
# The adaptive step controller needs to make one decision per step for the whole
# ensemble — we cannot take half a step for some particles and a full step for
# others. `scalar_reduce` collapses a possibly-Particles scalar down to a Float64
# using `pmaximum` (conservative: worst-case particle drives step size).
#
# Performance note: `pmean` would be a less conservative compromise and could
# reduce the step count under high-variance inputs; switch if MCM propagations
# become too slow under tight tolerances. Keep `scalar_reduce` as the single
# switch point.

scalar_reduce(x::AbstractFloat) = Float64(x)
scalar_reduce(x::Integer) = Float64(x)
# Catch-all for Particles-like types (ducked via a `.particles` field).  This
# avoids importing MonteCarloMeasurements into path-integral.jl.  MCM's
# `Particles` exposes `x.particles::Vector{Float64}` whose `maximum` is pmaximum.
function scalar_reduce(x)
    if hasfield(typeof(x), :particles)
        return Float64(maximum(getfield(x, :particles)))
    end
    return Float64(x)
end

# ----------------------------
# One exponential-midpoint step
# ----------------------------

"""
    sinhc()

sinhc(μ) = sinh(μ)/μ, with Taylor expansion around μ=0 to avoid numerical issues
(catastrophic cancellation).
"""
function sinhc(μ)
    # Elementwise arithmetic; lifts through MCM Particles entries automatically.
    # The small-μ branch uses the Taylor series; large-μ uses the analytic ratio.
    # Under MCM we cannot branch per-particle, so we reduce |μ| to a scalar first.
    if scalar_reduce(abs(μ)) < 1e-8
        μ2 = μ * μ
        return 1 + μ2 / 6 + μ2 * μ2 / 120 + μ2 * μ2 * μ2 / 5040
    end
    return sinh(μ) / μ
end

"""
    exp_jones_generator(A)

Closed-form exponential for a `2×2` Jones generator.

For an exactly traceless matrix `A`, Cayley-Hamilton gives `A^2 = -det(A) I`, so

`exp(A) = cosh(μ) I + sinh(μ)/μ * A`, where `μ^2 = -det(A)`.

This implementation also removes any tiny numerical trace drift by factoring out
`exp(tr(A)/2)` and applying the closed form to the centered traceless part.

Accepts any 2×2 AbstractMatrix with a Complex eltype; lifts elementwise through
MCM Particles entries.
"""
function exp_jones_generator(A::AbstractMatrix)
    @assert size(A) == (2, 2) "exp_jones_generator expects a 2x2 matrix"

    halftr = 0.5 * (A[1, 1] + A[2, 2])

    b11 = A[1, 1] - halftr
    b12 = A[1, 2]
    b21 = A[2, 1]
    b22 = A[2, 2] - halftr

    μ2 = -(b11 * b22 - b12 * b21)
    μ = sqrt(μ2)
    pref = exp(halftr)
    c = cosh(μ)
    s = sinhc(μ)

    T = promote_type(typeof(pref), typeof(c), typeof(s), typeof(b11))
    out = Matrix{T}(undef, 2, 2)
    out[1, 1] = pref * (c + s * b11)
    out[1, 2] = pref * (s * b12)
    out[2, 1] = pref * (s * b21)
    out[2, 2] = pref * (c + s * b22)
    return out
end

function exp_midpoint_step(K, s::Float64, h::Float64, J::AbstractMatrix)
    M = K(s + 0.5h)
    return exp_jones_generator(h * M) * J
end

"""
    exp_block_upper_triangular_2x2(hM, hMω)

Closed-form exponential for the block-upper-triangular 4×4 matrix

    A = [ M   Mω ]
        [ 0   M  ]

scaled by `h`.  The result is

    exp(h·A) = [ E   F ]
               [ 0   E ]

where `E = exp(hM)` and `F = h · φ(hM, Mω)` with `φ` the Frechet derivative
of `exp` along `Mω`:

    φ(A, B) = ∫₀¹ exp((1-τ)A) · B · exp(τA) dτ

For a 2×2 `M = c·I + M̃` with `M̃` traceless and `μ² = -det(M̃)`:

    exp(A) = exp(c) · [cosh(μ) I + sinhc(μ) M̃]

The Frechet derivative admits a closed form in the 2×2 case (see Al-Mohy &
Higham, "Computing the Fréchet derivative of the matrix exponential", 2009,
§5).  For `A = c·I + Ã` with `Ã` traceless and eigenvalues ±μ, one can write

    φ(A, B) = α₀·B + α₁·(ÃB + BÃ) + α₂·ÃBÃ

with scalar coefficients that depend only on `c` and `μ`.

We derive the coefficients directly: because the trace part commutes, it
factors as `exp(c)` on both `E` and `F`, leaving a traceless problem.  For a
traceless 2×2 `Ã` with eigenvalues ±μ, the identities `Ã² = μ² I` and the
representation `exp(Ã) = cosh(μ) I + sinhc(μ) Ã` let us evaluate the integral
in closed form using the basis {B, ÃB + BÃ, ÃBÃ}:

    exp(Ã) − I = (cosh(μ) − 1) I + sinhc(μ) Ã
    φ(Ã, B)   = β₀ B + β₁ (ÃB + BÃ) + β₂ ÃBÃ

Direct computation (expanding `exp((1-τ)Ã) = cosh((1-τ)μ)I + sinhc((1-τ)μ)(1-τ)Ã`
and the mirror factor, then integrating over τ ∈ [0,1]) gives:

    β₀ = sinhc(μ)                                    (coefficient of B)
    β₁ = (cosh(μ) − sinhc(μ)) / (2 μ²)                (coefficient of ÃB + BÃ)
    β₂ = (sinhc(μ) − cosh(μ) + μ² sinhc_int(μ)) / μ⁴  (coefficient of ÃBÃ)

where `sinhc_int(μ) = ∫₀¹ τ·sinhc((1-τ)μ)·sinhc(τμ) dτ`.  Rather than carry
that integral in closed form, we use the alternative decomposition based on
`Ã² = μ² I`: any polynomial in Ã on a 2×2 matrix reduces to `αI + βÃ`, so the
full `φ(Ã, B)` collapses to terms in `{B, ÃB, BÃ, ÃBÃ}` whose scalar weights
integrate to elementary functions of μ.

Because this derivation is fiddly, the MCM use case motivates it (LinearAlgebra's
generic `exp` relies on pivoting and bounds-style operations that do not lift
through `Particles`).  Instead of carrying the Frechet closed form symbolically,
we fall back to a numerically exact construction: compute `E = exp(hM)` via
`exp_jones_generator`, then compute `F` via Van Loan's identity

    exp([ hM   hMω ]) = [ exp(hM)   F       ]
        [ 0    hM  ]   [ 0         exp(hM) ]

and extract `F` by forming the 4×4 exponential through the 2×2 closed form
applied to `hM` and `hM'` plus an explicit integral representation.

Implementation: use Higham's scaling & squaring of the 2×2 closed form,
combined with the identity that for A block-upper-triangular with equal
diagonal blocks,

    exp(s·A)^2 = [ exp(2sM)    2·exp(sM)·F_s ]
                 [ 0           exp(2sM)      ]

where `F_s` is the off-diagonal of `exp(s·A)`.  We use s = 1 directly:

    F = ∑_{k=1}^∞ (hM)^(k-1)/k! · hMω · something...

which is equivalent to the series

    F = h · ∑_{k=0}^∞ (hM)^k · hMω · (hM)^j ... /(j+k+1)! ...

This is what Al-Mohy & Higham eq. (10.15) in Higham's "Functions of Matrices"
computes via scaling & squaring.  For 2×2 blocks the series closes in finite
form because of `M̃² = μ² I` for the traceless part.  Rather than re-derive,
we compute `F` by the Kronecker-product identity:

    vec(F) = (I ⊗ I) · L_A(hMω) · vec(I)
           ≈ h·sinhc(hM/2)  applied bilinearly

Concretely: factor `M = c I + M̃` with `M̃` traceless and `μ² = −det(M̃)`.
Then `exp(hM) = exp(hc) · [cosh(μh)·I + sinhc(μh)·M̃·h]` (note: we fold `h`
into the coefficient so `μh` is dimensionless).  The Frechet derivative along
`Mω` of the scalar-function part is then elementary, and we get:

    F = h · exp(hc) · (
            φ₀(hc, μh) · Mω
          + φ₁(hc, μh) · (M̃·Mω + Mω·M̃) · h
          + φ₂(hc, μh) · M̃·Mω·M̃ · h²
        )

with

    φ₀(a, b) = sinhc(b)
    φ₁(a, b) = (cosh(b) - sinhc(b)) / (2 b²)
    φ₂(a, b) = (b·sinh(b) - 3·cosh(b) + 3·sinhc(b)) / (2 b⁴)

evaluated at a = hc (trace part absorbed into `exp(hc)` prefactor) and
b = μh.  Small-b branches use Taylor expansions for numerical stability.
"""
# Compute the 2×2 product C = A · B given entries of A and B as tuples.
# Returns a 2×2 Matrix whose eltype is promoted from A and B entries so it
# lifts cleanly through MCM Particles.
@inline function _mul2x2(a11, a12, a21, a22, b11, b12, b21, b22)
    c11 = a11 * b11 + a12 * b21
    c12 = a11 * b12 + a12 * b22
    c21 = a21 * b11 + a22 * b21
    c22 = a21 * b12 + a22 * b22
    return c11, c12, c21, c22
end

function exp_block_upper_triangular_2x2(hM::AbstractMatrix, hMω::AbstractMatrix)
    @assert size(hM) == (2, 2) && size(hMω) == (2, 2)

    E = exp_jones_generator(hM)

    # Split hM = (half-trace)·I + hM̃, with hM̃ traceless.
    halftr = 0.5 * (hM[1, 1] + hM[2, 2])
    m11 = hM[1, 1] - halftr
    m12 = hM[1, 2]
    m21 = hM[2, 1]
    m22 = hM[2, 2] - halftr

    # μ² = −det(hM̃) ; since hM̃ already carries the h factor, μ is dimensionless.
    μ2 = -(m11 * m22 - m12 * m21)
    μ = sqrt(μ2)

    pref = exp(halftr)
    φ0, φ1, φ2 = _frechet_2x2_coeffs(μ, μ2)

    # F = pref · ( φ0·hMω + φ1·(M̃·hMω + hMω·M̃) + φ2·(M̃·hMω·M̃) )
    # where M̃ = hM - halftr·I (carries the h factor).
    w11, w12, w21, w22 = hMω[1,1], hMω[1,2], hMω[2,1], hMω[2,2]

    # M̃·Mω
    p11, p12, p21, p22 = _mul2x2(m11, m12, m21, m22, w11, w12, w21, w22)
    # Mω·M̃
    q11, q12, q21, q22 = _mul2x2(w11, w12, w21, w22, m11, m12, m21, m22)
    # M̃·Mω·M̃ = (M̃·Mω)·M̃
    r11, r12, r21, r22 = _mul2x2(p11, p12, p21, p22, m11, m12, m21, m22)

    T = promote_type(typeof(pref), typeof(φ0), typeof(w11), typeof(p11))
    F = Matrix{T}(undef, 2, 2)
    F[1, 1] = pref * (φ0 * w11 + φ1 * (p11 + q11) + φ2 * r11)
    F[1, 2] = pref * (φ0 * w12 + φ1 * (p12 + q12) + φ2 * r12)
    F[2, 1] = pref * (φ0 * w21 + φ1 * (p21 + q21) + φ2 * r21)
    F[2, 2] = pref * (φ0 * w22 + φ1 * (p22 + q22) + φ2 * r22)
    return E, F
end

# Fréchet-derivative scalar coefficients for exp along Ã on a 2×2 traceless matrix
# Ã with eigenvalues ±μ (so Ã² = μ²·I).  The derivative expands as
#     L_{c·I + Ã}(B) = exp(c) · [ φ0·B + φ1·(Ã·B + B·Ã) + φ2·Ã·B·Ã ]
# with, by direct integration of exp((1-τ)Ã)·B·exp(τÃ) over τ ∈ [0,1] using
# product-to-sum identities and Ã² = μ²·I:
#     φ0 = (cosh(μ) + sinhc(μ)) / 2                            [coeff of B]
#     φ1 = sinh(μ) / (2 μ)         = sinhc(μ) / 2              [coeff of Ã·B + B·Ã]
#     φ2 = (cosh(μ) − sinhc(μ)) / (2 μ²)                       [coeff of Ã·B·Ã]
# where sinhc(μ) = sinh(μ)/μ.  Small-|μ| branches from Taylor series in μ².
function _frechet_2x2_coeffs(μ, μ2)
    absμ = scalar_reduce(abs(μ))
    if absμ < 1e-4
        μ4 = μ2 * μ2
        μ6 = μ4 * μ2
        # φ0 = (cosh(μ) + sinhc(μ)) / 2
        #    = 1 + μ²/3 + μ⁴/30 + μ⁶/840 + ...   (after collecting)
        # Derivation:
        #   cosh(μ) = 1 + μ²/2 + μ⁴/24 + μ⁶/720 + ...
        #   sinhc(μ) = 1 + μ²/6 + μ⁴/120 + μ⁶/5040 + ...
        #   sum/2   = 1 + μ²·(1/2+1/6)/2 + μ⁴·(1/24+1/120)/2 + μ⁶·(1/720+1/5040)/2
        #           = 1 + μ²/3 + μ⁴/30 + μ⁶/840
        φ0 = 1 + μ2 / 3 + μ4 / 30 + μ6 / 840
        # φ1 = sinhc(μ) / 2
        φ1 = (1 + μ2 / 6 + μ4 / 120 + μ6 / 5040) / 2
        # φ2 = (cosh(μ) − sinhc(μ)) / (2 μ²)
        #    = (μ²/2 + μ⁴/24 + μ⁶/720 − μ²/6 − μ⁴/120 − μ⁶/5040) / (2 μ²)
        #    = (μ²·1/3 + μ⁴·(1/24−1/120) + μ⁶·(1/720−1/5040)) / (2 μ²)
        #    = 1/6 + μ²/30 + μ⁴/840 + ...
        φ2 = 1/6 + μ2 / 30 + μ4 / 840
    else
        s = sinh(μ)
        c = cosh(μ)
        shc = s / μ
        φ0 = (c + shc) / 2
        φ1 = shc / 2
        φ2 = (c - shc) / (2 * μ2)
    end
    return φ0, φ1, φ2
end

"""
    exp_sensitivity_midpoint_step(K, Kω, s, h, J, G)

One midpoint step for the coupled system `dJ/ds = K·J`, `dG/ds = Kω·J + K·G`.

Uses a closed-form exponential of the 4×4 block-upper-triangular generator

    A = [ K   Kω ]
        [ 0   K  ]

via `exp_block_upper_triangular_2x2`.  The motivation for avoiding
`LinearAlgebra.exp` on the 4×4 matrix is MCM compatibility: the generic Padé
approximation relies on operations (`typemax`, `isfinite` on matrix norms,
pivoting) that do not lift through `Particles`.  The closed form is also
faster and more accurate than the generic Padé even for pure Float64 inputs.
"""
function exp_sensitivity_midpoint_step(
    K,
    Kω,
    s::Float64,
    h::Float64,
    J::AbstractMatrix,
    G::AbstractMatrix
)
    M = K(s + 0.5h)
    Mω = Kω(s + 0.5h)

    E, F = exp_block_upper_triangular_2x2(h * M, h * Mω)

    # Y1 = exp(h·A) · Y0 with Y0 = [G; J] stacked as block column:
    #   exp(h·A) · [G; J] = [E·G + F·J; E·J]
    J_new = E * J
    G_new = E * G + F * J
    return J_new, G_new
end

# ----------------------------
# Matrix error metric
# Global phase insensitive
# ----------------------------

function phase_insensitive_error(A::AbstractMatrix, B::AbstractMatrix)
    # Remove best common phase from B relative to A.
    # Under MCM, opnorm of a Particles-entry matrix hits LinearAlgebra code paths
    # that don't lift (e.g. typemax); fall back to a Frobenius norm, which is
    # just arithmetic and lifts elementwise.
    α = tr(A' * B)
    absα = abs(α)
    if scalar_reduce(absα) == 0
        return _frobenius_norm(A - B)
    end
    ϕ = α / absα
    return _frobenius_norm(A - ϕ * B)
end

_frobenius_norm(M::AbstractMatrix) = sqrt(sum(abs2, M))

function align_global_phase(A::AbstractMatrix, B::AbstractMatrix)
    α = tr(A' * B)
    absα = abs(α)
    if scalar_reduce(absα) == 0
        return one(eltype(A)) + zero(eltype(A))
    end
    return α / absα
end

function sensitivity_phase_insensitive_error(
    J_ref::AbstractMatrix,
    G_ref::AbstractMatrix,
    J_cmp::AbstractMatrix,
    G_cmp::AbstractMatrix
)
    ϕ = align_global_phase(J_ref, J_cmp)
    errJ = _frobenius_norm(J_ref - ϕ * J_cmp)
    errG = _frobenius_norm(G_ref - ϕ * G_cmp)
    return max(scalar_reduce(errJ), scalar_reduce(errG))
end

# ----------------------------
# Adaptive step-doubling on one smooth interval
# ----------------------------

struct PropagatorStats
    accepted_steps::Int
    rejected_steps::Int
end

function propagate_interval!(
    K,
    s0::Float64,
    s1::Float64,
    J0::AbstractMatrix;
    rtol::Float64 = 1e-9,
    atol::Float64 = 1e-12,
    h_init::Float64 = (s1 - s0) / 100,
    h_min::Float64 = 1e-12,
    h_max::Float64 = s1 - s0,
    safety::Float64 = 0.9,
    growth_max::Float64 = 2.0,
    shrink_min::Float64 = 0.2
)
    s = s0
    J = copy(J0)
    h = min(h_init, s1 - s0)

    accepted = 0
    rejected = 0

    while s < s1
        h = min(h, s1 - s)
        @assert h >= h_min "Step size underflow near s=$s"

        # one full step
        J_full = exp_midpoint_step(K, s, h, J)

        # two half steps
        J_half = exp_midpoint_step(K, s, 0.5h, J)
        J_twohalf = exp_midpoint_step(K, s + 0.5h, 0.5h, J_half)

        # Collapse the error metric to a scalar here.  Under MCM, pmaximum is
        # conservative (worst-case particle picks the step); pmean would be a
        # cheaper compromise at the cost of occasionally under-refining.
        err_abs = scalar_reduce(phase_insensitive_error(J_full, J_twohalf))
        scale = max(
            scalar_reduce(_frobenius_norm(J_twohalf)),
            scalar_reduce(_frobenius_norm(J_full)),
        )
        tol = atol + rtol * scale

        if err_abs <= tol
            # accept better solution
            s += h
            J = J_twohalf
            accepted += 1

            # midpoint method is 2nd order, step-doubling estimate ~ O(h^3)
            if err_abs == 0.0
                fac = growth_max
            else
                fac = safety * (tol / err_abs)^(1/3)
                fac = clamp(fac, 1.0, growth_max)
            end
            h = min(h * fac, h_max)
        else
            rejected += 1
            fac = safety * (tol / err_abs)^(1/3)
            fac = clamp(fac, shrink_min, 0.5)
            h *= fac
            @assert h >= h_min "Step size fell below h_min near s=$s"
        end
    end

    return J, PropagatorStats(accepted, rejected)
end

function propagate_interval_sensitivity!(
    K,
    Kω,
    s0::Float64,
    s1::Float64,
    J0::AbstractMatrix;
    G0::AbstractMatrix = zeros(eltype(J0), 2, 2),
    rtol::Float64 = 1e-9,
    atol::Float64 = 1e-12,
    h_init::Float64 = (s1 - s0) / 100,
    h_min::Float64 = 1e-12,
    h_max::Float64 = s1 - s0,
    safety::Float64 = 0.9,
    growth_max::Float64 = 2.0,
    shrink_min::Float64 = 0.2
)
    s = s0
    J = copy(J0)
    G = copy(G0)
    h = min(h_init, s1 - s0)

    accepted = 0
    rejected = 0

    while s < s1
        h = min(h, s1 - s)
        @assert h >= h_min "Step size underflow near s=$s"

        J_full, G_full = exp_sensitivity_midpoint_step(K, Kω, s, h, J, G)

        J_half, G_half = exp_sensitivity_midpoint_step(K, Kω, s, 0.5h, J, G)
        J_twohalf, G_twohalf = exp_sensitivity_midpoint_step(K, Kω, s + 0.5h, 0.5h, J_half, G_half)

        # scalar_reduce: under MCM, pmaximum is conservative; pmean would be
        # a speed compromise — see note at top of this file.
        err_abs = sensitivity_phase_insensitive_error(J_full, G_full, J_twohalf, G_twohalf)
        scale = max(
            scalar_reduce(_frobenius_norm(J_full)),
            scalar_reduce(_frobenius_norm(J_twohalf)),
            scalar_reduce(_frobenius_norm(G_full)),
            scalar_reduce(_frobenius_norm(G_twohalf)),
        )
        tol = atol + rtol * scale

        if err_abs <= tol
            s += h
            J = J_twohalf
            G = G_twohalf
            accepted += 1

            if err_abs == 0.0
                fac = growth_max
            else
                fac = safety * (tol / err_abs)^(1/3)
                fac = clamp(fac, 1.0, growth_max)
            end
            h = min(h * fac, h_max)
        else
            rejected += 1
            fac = safety * (tol / err_abs)^(1/3)
            fac = clamp(fac, shrink_min, 0.5)
            h *= fac
            @assert h >= h_min "Step size fell below h_min near s=$s"
        end
    end

    return J, G, PropagatorStats(accepted, rejected)
end

# ----------------------------
# Domain decomposition with optional jump matrices
# ----------------------------

"""
    propagate_piecewise(K, breaks; jumps=Dict(), kwargs...)

Propagate dJ/ds = K(s) J from breaks[1] to breaks[end].

Arguments
---------
- K: callable generator
- breaks: sorted vector [s0, s1, ..., sN]
- jumps: optional Dict{Float64, Matrix{ComplexF64}} giving lumped jump matrices
         applied immediately AFTER reaching that breakpoint.

Returns
-------
- J_final
- stats_per_interval
"""
function propagate_piecewise(
    K,
    breaks::AbstractVector{<:Real};
    jumps::AbstractDict = Dict{Float64, Matrix{ComplexF64}}(),
    J0::AbstractMatrix = Matrix{ComplexF64}(I, 2, 2),
    verbose::Bool = true,
    kwargs...
)
    @assert issorted(breaks)
    J = copy(J0)
    stats = PropagatorStats[]
    n = length(breaks) - 1
    t_start = time()
    s_done = 0.0

    for i in 1:n
        sL = scalar_reduce(breaks[i])
        sR = scalar_reduce(breaks[i+1])

        J, st = propagate_interval!(K, sL, sR, J; kwargs...)
        push!(stats, st)
        s_done += sR - sL

        if verbose
            elapsed = time() - t_start
            rate = elapsed > 0 ? s_done / elapsed : 0.0
            @printf(stderr, "\r  propagate_piecewise: interval %d/%d (%.1f%%)  %.2f m/s  ", i, n, 100.0*i/n, rate)
        end

        if haskey(jumps, sR)
            J = jumps[sR] * J
        end
    end

    if verbose
        println(stderr)
    end

    return J, stats
end

function propagate_fiber(
    f::Fiber;
    λ_m::Real,
    jumps::AbstractDict = Dict{Float64, Matrix{ComplexF64}}(),
    kwargs...
)
    return propagate_piecewise(
        generator_K(f, λ_m),
        fiber_breakpoints(f);
        jumps = jumps,
        kwargs...,
    )
end

"""
    propagate_piecewise_sensitivity(K, Kω, breaks; jumps=Dict(), jump_omegas=Dict(), kwargs...)

Propagate the coupled system
`dJ/ds = K(s) J`,
`dG/ds = Kω(s) J + K(s) G`
from `breaks[1]` to `breaks[end]`.

Returns
-------
- `J_final`
- `G_final = ∂ωJ_final`
- `stats_per_interval`
"""
function propagate_piecewise_sensitivity(
    K,
    Kω,
    breaks::AbstractVector{<:Real};
    jumps::AbstractDict = Dict{Float64, Matrix{ComplexF64}}(),
    jump_omegas::AbstractDict = Dict{Float64, Matrix{ComplexF64}}(),
    J0::AbstractMatrix = Matrix{ComplexF64}(I, 2, 2),
    G0::AbstractMatrix = zeros(eltype(J0), 2, 2),
    verbose::Bool = true,
    kwargs...
)
    @assert issorted(breaks)
    J = copy(J0)
    G = copy(G0)
    stats = PropagatorStats[]
    n = length(breaks) - 1
    t_start = time()
    s_done = 0.0

    for i in 1:n
        sL = scalar_reduce(breaks[i])
        sR = scalar_reduce(breaks[i+1])

        J, G, st = propagate_interval_sensitivity!(K, Kω, sL, sR, J; G0 = G, kwargs...)
        push!(stats, st)
        s_done += sR - sL

        if verbose
            elapsed = time() - t_start
            rate = elapsed > 0 ? s_done / elapsed : 0.0
            @printf(stderr, "\r  propagate_piecewise_sensitivity: interval %d/%d (%.1f%%)  %.2f m/s  ", i, n, 100.0*i/n, rate)
        end

        if haskey(jumps, sR)
            J_pre = J
            G_pre = G
            J_jump = jumps[sR]
            G_jump = haskey(jump_omegas, sR) ? jump_omegas[sR] : zeros(eltype(J_pre), 2, 2)
            J = J_jump * J_pre
            G = G_jump * J_pre + J_jump * G_pre
        end
    end

    if verbose
        println(stderr)
    end

    return J, G, stats
end

function propagate_fiber_sensitivity(
    f::Fiber;
    λ_m::Real,
    jumps::AbstractDict = Dict{Float64, Matrix{ComplexF64}}(),
    jump_omegas::AbstractDict = Dict{Float64, Matrix{ComplexF64}}(),
    kwargs...
)
    return propagate_piecewise_sensitivity(
        generator_K(f, λ_m),
        generator_Kω(f, λ_m),
        fiber_breakpoints(f);
        jumps = jumps,
        jump_omegas = jump_omegas,
        kwargs...
    )
end

function pmd_generator(J::AbstractMatrix, G::AbstractMatrix; hermitianize::Bool = true)
    H = -1im * (J \ G)
    if hermitianize
        H = 0.5 * (H + H')
    end
    return H
end

"""
    output_dgd(J, G; hermitianize=true)

Differential group delay as the spread of real eigenvalues of the PMD generator
`H = -i·J⁻¹·G`.

Note: under MCM, `eigvals` on a `Matrix{<:Particles}` hits LinearAlgebra
routines that do not lift elementwise (typemax, pivoting).  For uncertain
propagations, compute the mean DGD per-sample via `bymap` at the call site,
or extract the Hermitian eigenvalues by hand for the 2×2 case (closed form:
eigenvalues of a 2×2 Hermitian matrix are `(tr/2) ± √((tr/2)² − det)`).
"""
function output_dgd(J::AbstractMatrix, G::AbstractMatrix; hermitianize::Bool = true)
    H = pmd_generator(J, G; hermitianize = hermitianize)
    λ = eigvals(H)
    return maximum(real.(λ)) - minimum(real.(λ))
end

"""
    output_dgd_2x2(J, G; hermitianize=true)

Closed-form DGD for a 2×2 PMD generator.  Equivalent to `output_dgd` for
the 2×2 case but avoids `eigvals`, so it lifts through MCM `Particles`.
"""
function output_dgd_2x2(J::AbstractMatrix, G::AbstractMatrix; hermitianize::Bool = true)
    @assert size(J) == (2, 2) && size(G) == (2, 2)
    H = pmd_generator(J, G; hermitianize = hermitianize)
    # Real eigenvalues of 2×2 Hermitian: (tr ± √(tr² − 4·det)) / 2.
    # DGD is the spread = √(tr² − 4·det).
    tr_H = real(H[1, 1] + H[2, 2])
    det_H = real(H[1, 1] * H[2, 2] - H[1, 2] * H[2, 1])
    disc = tr_H * tr_H - 4 * det_H
    return sqrt(max(disc, zero(disc)))
end
