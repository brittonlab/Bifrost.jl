using LinearAlgebra


# Code Overview:
# - fiber-path.jl defines the fiber model, source specification, and breakpoint assembly.
# - generator_K and generator_Kω assemble the local Jones dynamics (fiber-path.jl);
#   fiber-level wrappers propagate_fiber and propagate_fiber_Kω also live there.
# - the numerical propagation stack consists of
#     - exp_jones_generator
#     - exp_midpoint_step
#     - phase_insensitive_error
#     - propagate_interval!
#     - propagate_piecewise
# - the DGD Kω stack consists of
#     - exp_Kω_midpoint_step
#     - Kω_phase_insensitive_error
#     - propagate_interval_Kω!
#     - propagate_piecewise_Kω
#     - output_dgd
#
# fiber-path-plot.jl sits on top of that stack; it provides visual diagnostics.

include("fiber-path.jl")
include("fiber-path-plot.jl")

# ----------------------------
# One exponential-midpoint step
# ----------------------------

"""
    sinhc()

sinhc(μ) = sinh(μ)/μ, with Taylor expansion around μ=0 to avoid numerical issues
(catastrophic cancellation).
"""
function sinhc(μ::ComplexF64)
    if abs(μ) < 1e-8
        μ2 = μ * μ
        return 1 + μ2 / 6 + μ2 * μ2 / 120 + μ2 * μ2 * μ2 / 5040
    end
    return sinh(μ) / μ
end

"""
    exp_jones_generator(A::Matrix{ComplexF64})

Closed-form exponential for a `2×2` Jones generator.

For an exactly traceless matrix `A`, Cayley-Hamilton gives `A^2 = -det(A) I`, so

`exp(A) = cosh(μ) I + sinh(μ)/μ * A`, where `μ^2 = -det(A)`.

This implementation also removes any tiny numerical trace drift by factoring out
`exp(tr(A)/2)` and applying the closed form to the centered traceless part.
"""
function exp_jones_generator(A::Matrix{ComplexF64})
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

    return ComplexF64[
        pref * (c + s * b11)  pref * (s * b12)
        pref * (s * b21)      pref * (c + s * b22)
    ]
end

function exp_midpoint_step(K, s::Float64, h::Float64, J::Matrix{ComplexF64})
    M = K(s + 0.5h)
    return exp_jones_generator(h * M) * J
end

function exp_Kω_midpoint_step(
    K,
    Kω,
    s::Float64,
    h::Float64,
    J::Matrix{ComplexF64},
    G::Matrix{ComplexF64}
)
    M = K(s + 0.5h)
    Mω = Kω(s + 0.5h)

    A = zeros(ComplexF64, 4, 4)
    A[1:2, 1:2] = M
    A[1:2, 3:4] = Mω
    A[3:4, 3:4] = M

    Y0 = vcat(G, J)
    Y1 = exp(h * A) * Y0
    return Y1[3:4, :], Y1[1:2, :]
end

# ----------------------------
# Matrix error metric
# Global phase insensitive
# ----------------------------

function phase_insensitive_error(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    # Remove best common phase from B relative to A
    α = tr(A' * B)
    if abs(α) == 0
        return opnorm(A - B)
    end
    ϕ = α / abs(α)
    return opnorm(A - ϕ * B)
end

function align_global_phase(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    α = tr(A' * B)
    if abs(α) == 0
        return 1.0 + 0.0im
    end
    return α / abs(α)
end

function Kω_phase_insensitive_error(
    J_ref::Matrix{ComplexF64},
    G_ref::Matrix{ComplexF64},
    J_cmp::Matrix{ComplexF64},
    G_cmp::Matrix{ComplexF64}
)
    ϕ = align_global_phase(J_ref, J_cmp)
    errJ = opnorm(J_ref - ϕ * J_cmp)
    errG = opnorm(G_ref - ϕ * G_cmp)
    return max(errJ, errG)
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
    J0::Matrix{ComplexF64};
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

        err_abs = phase_insensitive_error(J_full, J_twohalf)
        scale = max(opnorm(J_twohalf), opnorm(J_full))
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

# ----------------------------
# Adaptive step-doubling on one smooth interval
# ----------------------------

function propagate_interval_Kω!(
    K,
    Kω,
    s0::Float64,
    s1::Float64,
    J0::Matrix{ComplexF64};
    G0::Matrix{ComplexF64} = zeros(ComplexF64, 2, 2),
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

        J_full, G_full = exp_Kω_midpoint_step(K, Kω, s, h, J, G)

        J_half, G_half = exp_Kω_midpoint_step(K, Kω, s, 0.5h, J, G)
        J_twohalf, G_twohalf = exp_Kω_midpoint_step(K, Kω, s + 0.5h, 0.5h, J_half, G_half)

        err_abs = Kω_phase_insensitive_error(J_full, G_full, J_twohalf, G_twohalf)
        scale = max(opnorm(J_full), opnorm(J_twohalf), opnorm(G_full), opnorm(G_twohalf))
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
    breaks::Vector{Float64};
    jumps::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}(),
    J0::Matrix{ComplexF64} = Matrix{ComplexF64}(I, 2, 2),
    kwargs...
)
    @assert issorted(breaks)
    J = copy(J0)
    stats = PropagatorStats[]

    for i in 1:length(breaks)-1
        sL = breaks[i]
        sR = breaks[i+1]

        J, st = propagate_interval!(K, sL, sR, J; kwargs...)
        push!(stats, st)

        if haskey(jumps, sR)
            J = jumps[sR] * J
        end
    end

    return J, stats
end

"""
    propagate_piecewise_Kω(K, Kω, breaks; jumps=Dict(), jump_omegas=Dict(), kwargs...)

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
function propagate_piecewise_Kω(
    K,
    Kω,
    breaks::Vector{Float64};
    jumps::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}(),
    jump_omegas::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}(),
    J0::Matrix{ComplexF64} = Matrix{ComplexF64}(I, 2, 2),
    G0::Matrix{ComplexF64} = zeros(ComplexF64, 2, 2),
    kwargs...
)
    @assert issorted(breaks)
    J = copy(J0)
    G = copy(G0)
    stats = PropagatorStats[]

    for i in 1:length(breaks)-1
        sL = breaks[i]
        sR = breaks[i+1]

        J, G, st = propagate_interval_Kω!(K, Kω, sL, sR, J; G0 = G, kwargs...)
        push!(stats, st)

        if haskey(jumps, sR)
            J_pre = J
            G_pre = G
            J_jump = jumps[sR]
            G_jump = haskey(jump_omegas, sR) ? jump_omegas[sR] : zeros(ComplexF64, 2, 2)
            J = J_jump * J_pre
            G = G_jump * J_pre + J_jump * G_pre
        end
    end

    return J, G, stats
end

function pmd_generator(J::Matrix{ComplexF64}, G::Matrix{ComplexF64}; hermitianize::Bool = true)
    H = -1im * (J \ G)
    if hermitianize
        H = 0.5 * (H + H')
    end
    return H
end

function output_dgd(J::Matrix{ComplexF64}, G::Matrix{ComplexF64}; hermitianize::Bool = true)
    H = pmd_generator(J, G; hermitianize = hermitianize)
    λ = eigvals(H)
    return maximum(real.(λ)) - minimum(real.(λ))
end
