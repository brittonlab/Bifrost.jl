"""
Example Use:

    julia julia-port/demo-path-integral.jl

Writes output/adaptive-step-doubling.html
"""

include("path-integral.jl")

using LinearAlgebra

"""
    demo_adaptive_step_doubling(; output, rtol, atol, title)

Illustrate the adaptive step-doubling integrator on a smooth, noncommuting generator

    K(s) = α·i·σ_x·cos(π·s) + β·i·σ_z·sin(2π·s),   s ∈ [0, 2]

The two Pauli components oscillate at different frequencies, so ‖K(s)‖ varies along the path
and the integrator must work harder near the fast oscillation peaks. The output plot shows:

- **Top panel**: accepted (green) and rejected (red) step sizes vs position, with the
  generator norm ‖K(s)‖ as a shaded overlay — small steps should cluster where K varies fastest.
- **Bottom panel**: err/tol ratio for every trial, with the acceptance threshold at 1.

Calls `collect_adaptive_steps` and `write_adaptive_steps_plot` from `fiber-path-plot.jl`;
the production solver in `path-integral.jl` is not modified.
"""
function demo_adaptive_step_doubling(;
    output::AbstractString = joinpath(@__DIR__, "..", "output", "adaptive-step-doubling.html"),
    rtol::Float64 = 1e-6,
    atol::Float64 = 1e-9,
    title::AbstractString = "Adaptive step-doubling: noncommuting K(s)"
)
    SX = ComplexF64[0 1; 1 0]
    SZ = ComplexF64[1 0; 0 -1]
    α = 1.2
    β = 0.9
    s0, s1 = 0.0, 2.0

    K = s -> α * im * cos(π * s) .* SX + β * im * sin(2π * s) .* SZ
    K_norm = s -> opnorm(K(s))

    J0 = Matrix{ComplexF64}(I, 2, 2)
    J_final, records = collect_adaptive_steps(K, s0, s1, J0; rtol = rtol, atol = atol)

    n_acc = count(r.accepted for r in records)
    n_rej = count(!r.accepted for r in records)
    println("Accepted steps: $n_acc,  rejected: $n_rej")
    println("Final Jones matrix:\n", J_final)

    plot_path = write_adaptive_steps_plot(
        records, K_norm, s0, s1;
        output = output,
        title  = title,
        rtol   = rtol,
        atol   = atol
    )
    println("Wrote adaptive step-doubling plot to:\n", plot_path)
    return (; J_final, records, plot_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_adaptive_step_doubling()
end
