using Test
using LinearAlgebra

if !isdefined(Main, :sample_fiber_centerline)
    # path-integral.jl already includes fiber-path.jl and fiber-path-plot.jl
    include("../path-integral.jl")
end

const CENTERLINE_ATOL = 5e-8
const TANGENT_ATOL = 5e-13
const SAMPLE_COUNT = 8001

function canonical_bend_segment(seg)
    if length(seg) == 2
        R, θ = seg
        α = 0.0
    else
        R, α, θ = seg
    end

    @assert R > 0.0 "Bend radius must be positive"
    αeff = θ >= 0 ? α : α + π
    return (R = Float64(R), α = Float64(αeff), θ = abs(Float64(θ)))
end

function skew_matrix(u::AbstractVector{<:Real})
    return Float64[
         0.0   -u[3]   u[2]
         u[3]   0.0   -u[1]
        -u[2]   u[1]   0.0
    ]
end

function rotate_vector_about_axis(v::Vector{Float64}, u::Vector{Float64}, θ::Float64)
    c = cos(θ)
    s = sin(θ)
    return c .* v + s .* cross(u, v) + (1 - c) .* dot(u, v) .* u
end

function exact_endpoint_and_tangent(segments)
    r = zeros(3)
    T = [0.0, 0.0, 1.0]

    for raw_seg in segments
        seg = canonical_bend_segment(raw_seg)
        if seg.θ == 0.0
            continue
        end

        u = [-sin(seg.α), cos(seg.α), 0.0]
        U = skew_matrix(u)
        A =
            seg.R * seg.θ * Matrix{Float64}(I, 3, 3) +
            seg.R * (1 - cos(seg.θ)) * U +
            seg.R * (seg.θ - sin(seg.θ)) * (U * U)

        r .+= A * T
        T = rotate_vector_about_axis(T, u, seg.θ)
        T ./= norm(T)
    end

    return r, T
end

function build_centerline_test_fiber(segments)
    @assert !isempty(segments) "Use build_straight_test_fiber for the straight case"

    breaks = Float64[0.0]
    pieces_Rb = Function[]
    pieces_theta = Function[]

    for raw_seg in segments
        seg = canonical_bend_segment(raw_seg)
        push!(pieces_Rb, _ -> seg.R)
        push!(pieces_theta, _ -> seg.α)
        push!(breaks, breaks[end] + seg.R * seg.θ)
    end

    bend = BendSource(
        PiecewiseProfile(copy(breaks), copy(pieces_Rb)),
        PiecewiseProfile(copy(breaks), copy(pieces_theta)),
        k2 -> 0.0;
        breakpoints = copy(breaks)
    )
    twist = TwistSource(
        _ -> 0.0,
        _ -> 0.0;
        coverage = [(first(breaks), last(breaks))],
        breakpoints = [first(breaks), last(breaks)]
    )

    return Fiber(first(breaks), last(breaks), AbstractBirefringenceSource[bend, twist])
end

function build_straight_test_fiber(length::Real)
    L = Float64(length)
    bend = BendSource(
        _ -> Inf,
        _ -> 0.0,
        k2 -> 0.0;
        coverage = [(0.0, L)],
        breakpoints = [0.0, L]
    )
    twist = TwistSource(
        _ -> 0.0,
        _ -> 0.0;
        coverage = [(0.0, L)],
        breakpoints = [0.0, L]
    )
    return Fiber(0.0, L, AbstractBirefringenceSource[bend, twist])
end

function sample_final_centerline_state(fiber::Fiber; n::Int = SAMPLE_COUNT)
    path = sample_fiber_centerline(fiber, fiber.s_start, fiber.s_end; n = n)
    endpoint = [path.x[end], path.y[end], path.zc[end]]
    tangent = [path.tx[end], path.ty[end], path.tz[end]]
    return path, endpoint, tangent
end

const EXPLICIT_PATH_CASES = [
    (
        name = "single_quarter_circle_xz_plane",
        segments = [(1.0, 0.0, π / 2)],
        xyz = [1.0, 0.0, 1.0],
        tangent = [1.0, 0.0, 0.0],
        planar = true
    ),
    (
        name = "single_half_circle_xz_plane",
        segments = [(2.0, 0.0, π)],
        xyz = [4.0, 0.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = true
    ),
    (
        name = "two_quarters_make_half_circle",
        segments = [(1.0, 0.0, π / 2), (1.0, 0.0, π / 2)],
        xyz = [2.0, 0.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = true
    ),
    (
        name = "three_sixty_degree_bends_make_half_circle",
        segments = [(1.0, 0.0, π / 3), (1.0, 0.0, π / 3), (1.0, 0.0, π / 3)],
        xyz = [2.0, 0.0, 0.0],
        tangent = [0.0, 0.0, -1.0],
        planar = true
    ),
    (
        name = "single_quarter_circle_yz_plane",
        segments = [(1.0, π / 2, π / 2)],
        xyz = [0.0, 1.0, 1.0],
        tangent = [0.0, 1.0, 0.0],
        planar = false
    ),
    (
        name = "single_oblique_quarter_circle",
        segments = [(2.0, π / 4, π / 2)],
        xyz = [sqrt(2), sqrt(2), 2.0],
        tangent = [inv(sqrt(2)), inv(sqrt(2)), 0.0],
        planar = false
    ),
    (
        name = "orthogonal_quarters_lab_frame",
        segments = [(1.0, 0.0, π / 2), (1.0, π / 2, π / 2)],
        xyz = [1.0 + π / 2, 0.0, 1.0],
        tangent = [1.0, 0.0, 0.0],
        planar = false
    ),
    (
        name = "orthogonal_quarters_swapped_order",
        segments = [(1.0, π / 2, π / 2), (1.0, 0.0, π / 2)],
        xyz = [0.0, 1.0 + π / 2, 1.0],
        tangent = [0.0, 1.0, 0.0],
        planar = false
    )
]

@testset "fiber-path centerline geometry" begin
    @testset "Straight section" begin
        fiber = build_straight_test_fiber(2.5)
        path, endpoint, tangent = sample_final_centerline_state(fiber)

        @test endpoint ≈ [0.0, 0.0, 2.5] atol = CENTERLINE_ATOL
        @test tangent ≈ [0.0, 0.0, 1.0] atol = TANGENT_ATOL
        @test maximum(abs.(path.x)) <= CENTERLINE_ATOL
        @test maximum(abs.(path.y)) <= CENTERLINE_ATOL
    end

    @testset "Explicit endpoint cases" begin
        for case in EXPLICIT_PATH_CASES
            @testset "$(case.name)" begin
                exact_xyz, exact_tangent = exact_endpoint_and_tangent(case.segments)
                fiber = build_centerline_test_fiber(case.segments)
                path, sampled_xyz, sampled_tangent = sample_final_centerline_state(fiber)

                @test exact_xyz ≈ case.xyz atol = 1e-12
                @test exact_tangent ≈ case.tangent atol = 1e-12
                @test sampled_xyz ≈ exact_xyz atol = CENTERLINE_ATOL
                @test sampled_tangent ≈ exact_tangent atol = TANGENT_ATOL
                @test norm(sampled_tangent) ≈ 1.0 atol = TANGENT_ATOL
                @test path.s[end] ≈ fiber.s_end atol = eps(Float64)

                if case.planar
                    @test maximum(abs.(path.y)) <= CENTERLINE_ATOL
                end
            end
        end
    end

    @testset "Segmentation invariance" begin
        reference = [(1.0, 0.0, π / 2)]
        split_two = [(1.0, 0.0, π / 4), (1.0, 0.0, π / 4)]
        split_four = [(1.0, 0.0, π / 8), (1.0, 0.0, π / 8), (1.0, 0.0, π / 8), (1.0, 0.0, π / 8)]

        ref_xyz_exact, ref_t_exact = exact_endpoint_and_tangent(reference)
        ref_path, ref_xyz, ref_t = sample_final_centerline_state(build_centerline_test_fiber(reference))
        two_path, two_xyz, two_t = sample_final_centerline_state(build_centerline_test_fiber(split_two))
        four_path, four_xyz, four_t = sample_final_centerline_state(build_centerline_test_fiber(split_four))

        @test ref_xyz_exact ≈ [1.0, 0.0, 1.0] atol = 1e-12
        @test ref_t_exact ≈ [1.0, 0.0, 0.0] atol = 1e-12

        @test ref_xyz ≈ ref_xyz_exact atol = CENTERLINE_ATOL
        @test two_xyz ≈ ref_xyz_exact atol = CENTERLINE_ATOL
        @test four_xyz ≈ ref_xyz_exact atol = CENTERLINE_ATOL
        @test ref_t ≈ ref_t_exact atol = TANGENT_ATOL
        @test two_t ≈ ref_t_exact atol = TANGENT_ATOL
        @test four_t ≈ ref_t_exact atol = TANGENT_ATOL

        @test ref_xyz ≈ two_xyz atol = CENTERLINE_ATOL
        @test ref_xyz ≈ four_xyz atol = CENTERLINE_ATOL
        @test two_xyz ≈ four_xyz atol = CENTERLINE_ATOL
        @test maximum(abs.(ref_path.y)) <= CENTERLINE_ATOL
        @test maximum(abs.(two_path.y)) <= CENTERLINE_ATOL
        @test maximum(abs.(four_path.y)) <= CENTERLINE_ATOL
    end

    @testset "Hard exact oracle case" begin
        segments = [(1.0, π / 4, π / 3), (1.0, -π / 4, π / 3)]
        exact_xyz, exact_tangent = exact_endpoint_and_tangent(segments)
        fiber = build_centerline_test_fiber(segments)
        _, sampled_xyz, sampled_tangent = sample_final_centerline_state(fiber)

        @test sampled_xyz ≈ exact_xyz atol = CENTERLINE_ATOL
        @test sampled_tangent ≈ exact_tangent atol = TANGENT_ATOL
    end

    @testset "Order matters in 3D" begin
        case_a = [(1.0, 0.0, π / 2), (1.0, π / 2, π / 2)]
        case_b = [(1.0, π / 2, π / 2), (1.0, 0.0, π / 2)]

        _, xyz_a, _ = sample_final_centerline_state(build_centerline_test_fiber(case_a))
        _, xyz_b, _ = sample_final_centerline_state(build_centerline_test_fiber(case_b))

        @test norm(xyz_a - xyz_b) > 1.0
    end
end

@testset "FiberSpec temperature profile" begin
    xs = FiberCrossSection(
        GermaniaSilicaGlass(0.036),
        GermaniaSilicaGlass(0.0),
        8.2e-6,
        125e-6
    )

    function fiber_sources(spec)
        return build(spec).sources
    end

    @testset "default T_K is 297.15" begin
        spec = FiberSpec(0.0, 10.0; cross_section = xs)
        @test spec.T_K(0.0) ≈ 297.15
        @test spec.T_K(5.0) ≈ 297.15
        twist!(spec, 0.0, 10.0; rate = 0.1)
        src = only(fiber_sources(spec))
        @test src.T_K(3.0) ≈ 297.15
    end

    @testset "scalar global T_K override" begin
        spec = FiberSpec(0.0, 10.0; cross_section = xs, T_K = 310.0)
        @test spec.T_K(0.0) ≈ 310.0
        twist!(spec, 0.0, 10.0; rate = 0.1)
        src = only(fiber_sources(spec))
        @test src.T_K(5.0) ≈ 310.0
    end

    @testset "function global T_K override" begin
        spec = FiberSpec(0.0, 10.0; cross_section = xs, T_K = s -> 300.0 + s)
        @test spec.T_K(0.0) ≈ 300.0
        @test spec.T_K(4.0) ≈ 304.0
        bend!(spec, 0.0, 10.0; angle = π / 4, axis = 0.0)
        src = only(fiber_sources(spec))
        @test src.T_K(4.0) ≈ 304.0
    end

    @testset "per-segment T_K overrides spec T_K" begin
        spec = FiberSpec(0.0, 20.0; cross_section = xs, T_K = 310.0)
        twist!(spec, 0.0, 10.0; T_K = 280.0, rate = 0.1)   # explicit override
        twist!(spec, 10.0, 20.0; rate = 0.1)                # inherits spec T_K
        srcs = [s for s in fiber_sources(spec) if s isa TwistSource]
        @test length(srcs) == 2
        # Identify sources by their T_K value at an interior point
        vals = sort([srcs[1].T_K(5.0), srcs[2].T_K(15.0)])
        @test vals ≈ [280.0, 310.0]
    end
end
