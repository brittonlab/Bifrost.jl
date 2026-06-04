using Test
using LinearAlgebra
using Random
using Bifrost


# Verify the closed-form 4×4 block-upper-triangular exp against LinearAlgebra.exp.
# For A = [M Mω; 0 M] we should have exp(A) = [E F; 0 E] with E = exp(M).
# See docstring of exp_block_upper_triangular_2x2 for the Frechet-derivative closed form.

function reference_block_exp(M::Matrix{ComplexF64}, Mω::Matrix{ComplexF64})
    A = zeros(ComplexF64, 4, 4)
    A[1:2, 1:2] = M
    A[1:2, 3:4] = Mω
    A[3:4, 3:4] = M
    X = exp(A)
    return X[1:2, 1:2], X[1:2, 3:4]
end

@testset "Frechet 2×2 — closed form matches LinearAlgebra.exp" begin
    Random.seed!(42)
    tol = 1e-10
    fails = 0
    for trial in 1:500
        # Random complex 2×2 matrices with moderate norm so exp is well-conditioned
        M = randn(ComplexF64, 2, 2) * 0.3
        Mω = randn(ComplexF64, 2, 2) * 0.3
        E_ref, F_ref = reference_block_exp(M, Mω)
        E, F = exp_block_upper_triangular_2x2(M, Mω)
        if !(isapprox(E, E_ref; atol = tol, rtol = tol) &&
             isapprox(F, F_ref; atol = tol, rtol = tol))
            fails += 1
        end
    end
    @test fails == 0

    # Edge case: M traceless (halftr = 0)
    M = ComplexF64[0.1 0.2+1im; 0.3-0.5im -0.1]
    Mω = randn(ComplexF64, 2, 2) * 0.3
    E_ref, F_ref = reference_block_exp(M, Mω)
    E, F = exp_block_upper_triangular_2x2(M, Mω)
    @test isapprox(E, E_ref; atol = tol)
    @test isapprox(F, F_ref; atol = tol)

    # Edge case: very small μ (M ≈ c·I + tiny traceless part)
    M = ComplexF64[1e-10 + 0.5 1e-11; 1e-11 1e-10 + 0.5]
    Mω = ComplexF64[0.1 0.2; 0.3 0.4]
    E_ref, F_ref = reference_block_exp(M, Mω)
    E, F = exp_block_upper_triangular_2x2(M, Mω)
    @test isapprox(E, E_ref; atol = tol)
    @test isapprox(F, F_ref; atol = tol)

    # Edge case: M identically zero
    M = zeros(ComplexF64, 2, 2)
    Mω = ComplexF64[1 2; 3 4]
    E_ref, F_ref = reference_block_exp(M, Mω)
    E, F = exp_block_upper_triangular_2x2(M, Mω)
    @test isapprox(E, E_ref; atol = tol)
    @test isapprox(F, F_ref; atol = tol)
end
