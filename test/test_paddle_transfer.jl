using Test
using LinearAlgebra
using Bifrost

const PADDLE_BEND_SCALE = 7.5
const STATE_TOL = 1e-6

const H_STATE = ComplexF64[1.0 + 0.0im, 0.0 + 0.0im]
const V_STATE = ComplexF64[0.0 + 0.0im, 1.0 + 0.0im]
const D_STATE = ComplexF64[inv(sqrt(2)) + 0.0im, inv(sqrt(2)) + 0.0im]
const A_STATE = ComplexF64[inv(sqrt(2)) + 0.0im, -inv(sqrt(2)) + 0.0im]
const R_STATE = ComplexF64[inv(sqrt(2)) + 0.0im, -im * inv(sqrt(2))]
const L_STATE = ComplexF64[inv(sqrt(2)) + 0.0im, im * inv(sqrt(2))]

function state_phase_error(actual::Vector{ComplexF64}, expected::Vector{ComplexF64})
    actual_norm = actual / norm(actual)
    expected_norm = expected / norm(expected)
    α = dot(expected_norm, actual_norm)
    ϕ = abs(α) == 0 ? (1.0 + 0.0im) : α / abs(α)
    return norm(actual_norm - ϕ * expected_norm)
end

function paddle_breaks_and_values(paddles)
    @assert !isempty(paddles) "Need at least one paddle"

    breaks = Float64[0.0]
    values = Matrix{ComplexF64}[]

    for (radius_mm, turns, angle_deg) in paddles
        # The proposed unit-test vectors use the opposite sign convention from
        # the internal bend-axis angle used by the fiber generator assembly.
        θ = -deg2rad(angle_deg)
        k2 = inv(radius_mm)^2
        Δβ = PADDLE_BEND_SCALE * k2
        c2φ = cos(2 * θ)
        s2φ = sin(2 * θ)
        push!(values, ComplexF64[
             0.5im * Δβ * c2φ   0.5im * Δβ * s2φ
             0.5im * Δβ * s2φ  -0.5im * Δβ * c2φ
        ])
        push!(breaks, breaks[end] + 2π * radius_mm * turns)
    end

    return breaks, values
end

function propagate_test_state(paddles, input_state::Vector{ComplexF64})
    breaks, values = paddle_breaks_and_values(paddles)
    K = function (s::Real)
        idx = min(searchsortedlast(breaks, Float64(s)), length(values))
        return values[idx]
    end
    J, stats = propagate_piecewise(K, breaks; rtol = 1e-11, atol = 1e-13, h_init = 0.1)
    return J * input_state, J, stats
end

const PADDLE_TEST_CASES = [
    (
        name = "1P-1",
        input_state = H_STATE,
        paddles = [(30.0, 2.0, 45.0)],
        expected_state = V_STATE
    ),
    (
        name = "1P-2",
        input_state = H_STATE,
        paddles = [(30.0, 1.0, 45.0)],
        expected_state = R_STATE
    ),
    (
        name = "1P-3",
        input_state = D_STATE,
        paddles = [(30.0, 2.0, 0.0)],
        expected_state = A_STATE
    ),
    (
        name = "1P-4",
        input_state = H_STATE,
        paddles = [(60.0, 1.0, 45.0)],
        expected_state = ComplexF64[0.923880 + 0.0im, 0.0 - 0.382683im]
    ),
    (
        name = "2P-1",
        input_state = H_STATE,
        paddles = [(30.0, 1.0, 45.0), (30.0, 1.0, 45.0)],
        expected_state = V_STATE
    ),
    (
        name = "2P-2",
        input_state = H_STATE,
        paddles = [(30.0, 1.0, 45.0), (30.0, 2.0, 22.5)],
        expected_state = L_STATE
    ),
    (
        name = "2P-3",
        input_state = D_STATE,
        paddles = [(60.0, 1.0, 0.0), (30.0, 1.0, 45.0)],
        expected_state = ComplexF64[0.382683 + 0.0im, 0.923880 + 0.0im]
    ),
    (
        name = "2P-4",
        input_state = R_STATE,
        paddles = [(30.0, 2.0, 0.0), (30.0, 2.0, 22.5)],
        expected_state = R_STATE
    ),
    (
        name = "3P-1",
        input_state = H_STATE,
        paddles = [(30.0, 1.0, 45.0), (30.0, 2.0, 22.5), (30.0, 1.0, 45.0)],
        expected_state = H_STATE
    ),
    (
        name = "3P-2",
        input_state = D_STATE,
        paddles = [(30.0, 1.0, 0.0), (30.0, 2.0, 45.0), (30.0, 1.0, 0.0)],
        expected_state = D_STATE
    ),
    (
        name = "3P-3",
        input_state = H_STATE,
        paddles = [(60.0, 1.0, 30.0), (30.0, 2.0, 10.0), (20.0, 1.0, 70.0)],
        expected_state = ComplexF64[0.754499 + 0.0im, 0.159053 - 0.636736im]
    ),
    (
        name = "3P-4",
        input_state = R_STATE,
        paddles = [(30.0, 1.0, 15.0), (60.0, 1.0, 70.0), (30.0, 1.0, -20.0)],
        expected_state = ComplexF64[0.466695 + 0.0im, -0.845613 + 0.259104im]
    ),
    (
        name = "5P-1",
        input_state = H_STATE,
        paddles = [
            (30.0, 1.0, 45.0),
            (30.0, 1.0, 45.0),
            (30.0, 2.0, 45.0),
            (30.0, 1.0, 45.0),
            (30.0, 1.0, 45.0)
        ],
        expected_state = V_STATE
    ),
    (
        name = "5P-2",
        input_state = H_STATE,
        paddles = [
            (30.0, 2.0, 0.0),
            (30.0, 2.0, 0.0),
            (30.0, 1.0, 45.0),
            (30.0, 1.0, 45.0),
            (30.0, 2.0, 45.0)
        ],
        expected_state = H_STATE
    ),
    (
        name = "5P-3",
        input_state = H_STATE,
        paddles = [
            (60.0, 1.0, 10.0),
            (30.0, 2.0, 35.0),
            (20.0, 1.0, 70.0),
            (30.0, 1.0, -15.0),
            (60.0, 1.0, 40.0)
        ],
        expected_state = ComplexF64[0.603608 + 0.0im, -0.495796 + 0.624375im]
    ),
    (
        name = "5P-4",
        input_state = D_STATE,
        paddles = [
            (30.0, 1.0, 0.0),
            (60.0, 1.0, 25.0),
            (30.0, 2.0, 50.0),
            (20.0, 1.0, -10.0),
            (30.0, 1.0, 40.0)
        ],
        expected_state = ComplexF64[0.639963 + 0.0im, 0.767355 + 0.040164im]
    ),
    (
        name = "5P-5",
        input_state = R_STATE,
        paddles = [
            (30.0, 2.0, 45.0),
            (60.0, 1.0, 0.0),
            (30.0, 1.0, 30.0),
            (30.0, 2.0, 15.0),
            (20.0, 1.0, 75.0)
        ],
        expected_state = ComplexF64[0.354458 + 0.0im, -0.664321 - 0.658055im]
    )
]

@testset "Paddle transfer cases" begin
    for case in PADDLE_TEST_CASES
        @testset "$(case.name)" begin
            output_state, J, stats = propagate_test_state(case.paddles, case.input_state)

            @test size(J) == (2, 2)
            @test !isempty(stats)
            @test length(stats) == length(case.paddles)
            @test all(st -> st.accepted_steps > 0, stats)
            @test isapprox(norm(output_state), norm(case.input_state); atol = 1e-10)
            @test state_phase_error(output_state, case.expected_state) <= STATE_TOL
        end
    end
end
