using Test
using LinearAlgebra
using Bifrost


const SIGMA_X = ComplexF64[0 1; 1 0]
const SIGMA_Y = ComplexF64[0 -im; im 0]
const SIGMA_Z = ComplexF64[1 0; 0 -1]
const DGD_ATOL = 1e-9
const MATRIX_ATOL = 1e-9

pauli_axis(axis::AbstractVector{<:Real}) = axis[1] .* SIGMA_X + axis[2] .* SIGMA_Y + axis[3] .* SIGMA_Z

function piecewise_function(breaks::Vector{Float64}, pieces)
    @assert issorted(breaks)
    @assert length(pieces) == length(breaks) - 1
    return function (s::Real)
        idx = min(searchsortedlast(breaks, Float64(s)), length(pieces))
        return pieces[idx](s)
    end
end

function piecewise_constant_matrix(breaks::Vector{Float64}, values::Vector{Matrix{ComplexF64}})
    @assert issorted(breaks)
    @assert length(values) == length(breaks) - 1
    return function (s::Real)
        idx = min(searchsortedlast(breaks, Float64(s)), length(values))
        return values[idx]
    end
end

function propagate_dgd_case(K, Komega, breaks; rtol = 1e-11, atol = 1e-13, h_init = 0.05)
    J, G, stats = propagate_piecewise_sensitivity(K, Komega, breaks; rtol = rtol, atol = atol, h_init = h_init)
    return J, G, stats
end

accepted_steps(stats) = sum(st.accepted_steps for st in stats)

@testset "DGD exact cases" begin
    @testset "Zero Komega gives zero DGD" begin
        K = s -> 0.3im .* SIGMA_Z
        Komega = s -> zeros(ComplexF64, 2, 2)
        J, G, stats = propagate_dgd_case(K, Komega, [0.0, 2.0])

        @test accepted_steps(stats) > 0
        @test opnorm(J - Matrix{ComplexF64}(I, 2, 2)) > 1e-3
        @test opnorm(G) <= MATRIX_ATOL
        @test output_dgd(J, G) <= DGD_ATOL
    end

    @testset "Closed form along sigma_z" begin
        length_scale = 2.0
        beta = 0.4
        beta_prime = 0.7
        K = s -> beta * im .* SIGMA_Z
        Komega = s -> beta_prime * im .* SIGMA_Z
        J, G, _ = propagate_dgd_case(K, Komega, [0.0, length_scale])

        expected_pmd = length_scale * beta_prime .* SIGMA_Z
        @test pmd_generator(J, G) ≈ expected_pmd atol = MATRIX_ATOL
        @test output_dgd(J, G) ≈ 2 * abs(length_scale * beta_prime) atol = DGD_ATOL
    end

    @testset "Closed form along arbitrary axis" begin
        length_scale = 2.0
        beta = 0.4
        beta_prime = 0.7
        axis = [1.0, 1.0, 1.0]
        axis ./= norm(axis)
        axis_matrix = pauli_axis(axis)
        K = s -> beta * im .* axis_matrix
        Komega = s -> beta_prime * im .* axis_matrix
        J, G, _ = propagate_dgd_case(K, Komega, [0.0, length_scale])

        expected_pmd = length_scale * beta_prime .* axis_matrix
        @test pmd_generator(J, G) ≈ expected_pmd atol = MATRIX_ATOL
        @test output_dgd(J, G) ≈ 2 * abs(length_scale * beta_prime) atol = DGD_ATOL
    end

    @testset "Piecewise cancellation on one axis" begin
        breaks = [0.0, 1.0, 2.0]
        K = s -> 0.5im .* SIGMA_Z
        Komega = piecewise_constant_matrix(
            breaks,
            [
                0.8im .* SIGMA_Z,
                -0.8im .* SIGMA_Z
            ]
        )
        J, G, stats = propagate_dgd_case(K, Komega, breaks)

        @test length(stats) == 2
        @test pmd_generator(J, G) ≈ zeros(ComplexF64, 2, 2) atol = MATRIX_ATOL
        @test output_dgd(J, G) <= DGD_ATOL
    end

    @testset "Piecewise commuting accumulation" begin
        breaks = [0.0, 0.5, 1.25, 2.0]
        beta_vals = [0.3, 0.4, -0.2]
        beta_prime_vals = [0.25, -0.5, 1.25]
        K = piecewise_constant_matrix(breaks, [b * im .* SIGMA_Z for b in beta_vals])
        Komega = piecewise_constant_matrix(breaks, [bp * im .* SIGMA_Z for bp in beta_prime_vals])
        J, G, _ = propagate_dgd_case(K, Komega, breaks)

        expected_coeff =
            (breaks[2] - breaks[1]) * beta_prime_vals[1] +
            (breaks[3] - breaks[2]) * beta_prime_vals[2] +
            (breaks[4] - breaks[3]) * beta_prime_vals[3]

        @test pmd_generator(J, G) ≈ expected_coeff .* SIGMA_Z atol = MATRIX_ATOL
        @test output_dgd(J, G) ≈ 2 * abs(expected_coeff) atol = DGD_ATOL
    end
end

@testset "DGD invariants and regression" begin
    @testset "Constant basis rotation preserves DGD" begin
        breaks = [0.0, 0.5, 1.25, 2.0]
        K = piecewise_constant_matrix(
            breaks,
            [
                0.3im .* SIGMA_Z,
                0.4im .* SIGMA_Z,
                -0.2im .* SIGMA_Z
            ]
        )
        Komega = piecewise_constant_matrix(
            breaks,
            [
                0.25im .* SIGMA_Z,
                -0.5im .* SIGMA_Z,
                1.25im .* SIGMA_Z
            ]
        )
        J_ref, G_ref, _ = propagate_dgd_case(K, Komega, breaks)
        H_ref = pmd_generator(J_ref, G_ref)

        rotation = exp_jones_generator(0.37im .* SIGMA_Y)
        K_rot = s -> rotation * K(s) * rotation'
        Komega_rot = s -> rotation * Komega(s) * rotation'
        J_rot, G_rot, _ = propagate_dgd_case(K_rot, Komega_rot, breaks)
        H_rot = pmd_generator(J_rot, G_rot)

        @test H_rot ≈ rotation * H_ref * rotation' atol = MATRIX_ATOL
        @test output_dgd(J_rot, G_rot) ≈ output_dgd(J_ref, G_ref) atol = DGD_ATOL
    end

    @testset "Phase scaling leaves PMD and DGD unchanged" begin
        J, G, _ = propagate_dgd_case(s -> 0.4im .* SIGMA_Z, s -> 0.7im .* SIGMA_Z, [0.0, 2.0])
        phase = cis(0.73)

        @test pmd_generator(phase .* J, phase .* G) ≈ pmd_generator(J, G) atol = MATRIX_ATOL
        @test output_dgd(phase .* J, phase .* G) ≈ output_dgd(J, G) atol = DGD_ATOL
    end

    @testset "Noncommuting K with zero Komega still gives zero DGD" begin
        breaks = [0.0, 0.7, 1.4, 2.0]
        mixed_axis = (SIGMA_X .+ SIGMA_Z) ./ sqrt(2)
        K = piecewise_constant_matrix(
            breaks,
            [
                0.6im .* SIGMA_X,
                -0.2im .* SIGMA_Z,
                0.4im .* mixed_axis
            ]
        )
        Komega = s -> zeros(ComplexF64, 2, 2)
        J, G, _ = propagate_dgd_case(K, Komega, breaks)

        @test opnorm(J - Matrix{ComplexF64}(I, 2, 2)) > 1e-3
        @test opnorm(G) <= MATRIX_ATOL
        @test output_dgd(J, G) <= DGD_ATOL
    end

    @testset "Direct sensitivity agrees with finite difference" begin
        length_scale = 1.75
        breaks = [0.0, length_scale]
        a0, a1 = 0.31, 0.07
        b0, b1 = -0.22, 0.05
        c0, c1 = 0.18, -0.03

        function K_param(s, omega)
            return im * ((a0 + a1 * omega) * sin(0.7 * s)) .* SIGMA_X +
                   im * ((b0 + b1 * omega) * cos(0.4 * s)) .* SIGMA_Y +
                   im * (c0 + c1 * omega) .* SIGMA_Z
        end

        K_of_omega(omega) = s -> K_param(s, omega)
        Komega = s -> im * (a1 * sin(0.7 * s)) .* SIGMA_X +
                      im * (b1 * cos(0.4 * s)) .* SIGMA_Y +
                      im * c1 .* SIGMA_Z

        J_direct, G_direct, _ = propagate_dgd_case(
            K_of_omega(0.0),
            Komega,
            breaks;
            rtol = 1e-10,
            atol = 1e-12,
            h_init = 0.02
        )

        delta_omega = 1e-5
        J_plus, _ = propagate_piecewise(K_of_omega(delta_omega), breaks; rtol = 1e-10, atol = 1e-12, h_init = 0.02)
        J_minus, _ = propagate_piecewise(K_of_omega(-delta_omega), breaks; rtol = 1e-10, atol = 1e-12, h_init = 0.02)
        G_fd = (J_plus - J_minus) ./ (2 * delta_omega)

        @test opnorm(G_direct - G_fd) <= 1e-8
        @test opnorm(pmd_generator(J_direct, G_direct) - pmd_generator(J_direct, G_fd)) <= 1e-8
        @test abs(output_dgd(J_direct, G_direct) - output_dgd(J_direct, G_fd)) <= 1e-8
    end
end

# TODO: spin refactor — DGD fiber integration test depends on a spin overlay
# adding a breakpoint at s=2.0; pending per-segment-meta spin subsystem.
@testset "DGD fiber integration" begin
    @test_skip true
end
