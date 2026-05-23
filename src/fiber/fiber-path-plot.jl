using LinearAlgebra


function rotate_about_axis(v::Vector{Float64}, u::Vector{Float64}, α::Float64)
    return v * cos(α) + cross(u, v) * sin(α) + u * dot(u, v) * (1 - cos(α))
end

function stokes_from_jones(ψ::AbstractVector{ComplexF64})
    s0 = real(abs2(ψ[1]) + abs2(ψ[2]))
    if s0 == 0.0
        return (0.0, 0.0, 0.0)
    end
    s1 = real(abs2(ψ[1]) - abs2(ψ[2])) / s0
    coh = ψ[1] * conj(ψ[2])
    s2 = 2 * real(coh) / s0
    s3 = -2 * imag(coh) / s0
    return (s1, s2, s3)
end

""" Explain sampled_path_* methods here.

The sampled path methods provide utilities for working with a discretely sampled fiber path, 
  where the geometry and polarization state are defined at specific arc-length coordinates. These 
  methods allow for interpolation of scalar and vector fields along the path, as well as computation 
  of the Frenet-Serret frame, which describes the local geometric properties of the path. 
  
  - The `render_pol_circle` method generates a visual representation of the polarization state at a given 
  point along the path.
  - The `render_poincare_sphere` method creates a visualization of the polarization state on the Poincare sphere. 
  
sample_path_bracket, sasmpled_path_scalar, and sampled_path_vector are helper functions for interpolating values along 
the path based on the arc-length parameter `s`. The frenet_serret_frame function computes the tangent, normal, and 
binormal vectors at a given point on the path, which are essential for understanding the local geometry of the fiber.
"""

"""
    sampled_path_bracket(path, s)

Given a sampled path with arc-length parameter `s`, find the indices of the two sample points that 
bracket `s` and the interpolation parameter `t` for linear interpolation.
"""
function sampled_path_bracket(path, s::Real)
    ss = path.s
    @assert first(ss) <= s <= last(ss) "s out of bounds"

    sf = Float64(s)
    idx = searchsortedlast(ss, sf)
    if idx <= 0
        return (left = 1, right = 1, t = 0.0)
    elseif idx >= length(ss)
        return (left = length(ss), right = length(ss), t = 0.0)
    end

    sL = ss[idx]
    sR = ss[idx + 1]
    t = sR == sL ? 0.0 : (sf - sL) / (sR - sL)
    return (left = idx, right = idx + 1, t = t)
end

function sampled_path_scalar(path, field::Symbol, s::Real)
    br = sampled_path_bracket(path, s)
    vals = getproperty(path, field)
    return (1 - br.t) * vals[br.left] + br.t * vals[br.right]
end

function sampled_path_vector(path, fields::NTuple{3,Symbol}, s::Real)
    return [
        sampled_path_scalar(path, fields[1], s),
        sampled_path_scalar(path, fields[2], s),
        sampled_path_scalar(path, fields[3], s)
    ]
end

"""
    frenet_serret_frame(path, s; curvature_tol = 1e-8)

Return the Frenet-Serret frame of a sampled path at arc-length coordinate `s`.

The `path` input must provide `s`, `x`, `y`, and `zc` arrays, where `s` is the arc-length
parameter. On locally straight sections, where the classical Frenet normal is undefined, the
function falls back to the path's transported frame `e1/e2` when available.
"""
function frenet_serret_frame(path, s::Real; curvature_tol::Float64 = 1e-8)
    ss = path.s
    @assert first(ss) <= s <= last(ss) "s out of bounds"

    idx = searchsortedfirst(ss, Float64(s))
    i = clamp(idx, 2, length(ss) - 1)

    rL = sampled_path_vector(path, (:x, :y, :zc), ss[i - 1])
    r0 = sampled_path_vector(path, (:x, :y, :zc), ss[i])
    rR = sampled_path_vector(path, (:x, :y, :zc), ss[i + 1])
    dsL = ss[i] - ss[i - 1]
    dsR = ss[i + 1] - ss[i]

    t_raw = (rR - rL) / (ss[i + 1] - ss[i - 1])
    tangent = t_raw / norm(t_raw)
    d2r = 2 * (((rR - r0) / dsR) - ((r0 - rL) / dsL)) / (dsL + dsR)
    normal_raw = d2r - dot(d2r, tangent) * tangent
    curvature = norm(normal_raw)

    if curvature > curvature_tol
        normal = normal_raw / curvature
        binormal = cross(tangent, normal)
        binormal /= norm(binormal)
    elseif all(hasproperty(path, p) for p in (:e1x, :e1y, :e1z, :e2x, :e2y, :e2z))
        normal = [path.e1x[i], path.e1y[i], path.e1z[i]]
        binormal = [path.e2x[i], path.e2y[i], path.e2z[i]]
    else
        error("Frenet normal undefined on a straight section and no fallback frame is available")
    end

    return (; tangent, normal, binormal, curvature)
end

"""
    poincare_vector_representation(ψ)

Convert a Jones state vector `ψ` into a Poincare-sphere vector and its equatorial projection.

Returns a named tuple with:
- `sphere`: Cartesian `(x, y, z)` coordinates on the unit Poincare sphere
- `linear`: `(angle_rad, radius)` for the equatorial-plane projection, where `0` is
  horizontal, `π/2` is vertical, and `radius` is the degree of linear polarization
"""
function poincare_vector_representation(ψ::AbstractVector{ComplexF64})
    s1, s2, s3 = stokes_from_jones(ψ)
    angle_rad = mod(0.5 * atan(s2, s1), π)
    radius = hypot(s1, s2)
    return (
        sphere = (x = s1, y = s2, z = s3),
        linear = (angle_rad = angle_rad, radius = radius)
    )
end

"""
    render_pol_circle(path, s, rep; npts = 181)

Build Plotly-ready traces for a unit polarization circle centered on the fiber centerline at `s`.

The circle lies in the plane spanned by the Frenet-Serret normal and binormal vectors.
The arrow starts at the circle center and points along the equatorial projection of the
Poincare vector, with angle `rep.linear.angle_rad` and length `rep.linear.radius`.
"""
function render_pol_circle(path, s::Real, rep::NamedTuple; npts::Int = 181)
    center = sampled_path_vector(path, (:x, :y, :zc), s)
    if all(hasproperty(path, p) for p in (:nx, :ny, :nz, :bx, :by, :bz))
        n̂ = sampled_path_vector(path, (:nx, :ny, :nz), s)
        b̂ = sampled_path_vector(path, (:bx, :by, :bz), s)
    else
        fs = frenet_serret_frame(path, s)
        n̂ = fs.normal
        b̂ = fs.binormal
    end

    ts = collect(range(0.0, 2π, length = npts))
    circle_x = center[1] .+ cos.(ts) .* n̂[1] .+ sin.(ts) .* b̂[1]
    circle_y = center[2] .+ cos.(ts) .* n̂[2] .+ sin.(ts) .* b̂[2]
    circle_z = center[3] .+ cos.(ts) .* n̂[3] .+ sin.(ts) .* b̂[3]
    tick_angles = (0.0, π / 2, π, 3π / 2)
    tick_inner = 0.82
    tick_x = Float64[]
    tick_y = Float64[]
    tick_z = Float64[]
    for ϕ in tick_angles
        rim = cos(ϕ) .* n̂ .+ sin(ϕ) .* b̂
        inner = tick_inner .* rim
        append!(tick_x, [center[1] + inner[1], center[1] + rim[1], NaN])
        append!(tick_y, [center[2] + inner[2], center[2] + rim[2], NaN])
        append!(tick_z, [center[3] + inner[3], center[3] + rim[3], NaN])
    end

    ψ = rep.linear.angle_rad
    ρ = rep.linear.radius
    arrow_dir = cos(ψ) .* n̂ .+ sin(ψ) .* b̂
    tip = (
        x = center[1] + ρ * arrow_dir[1],
        y = center[2] + ρ * arrow_dir[2],
        z = center[3] + ρ * arrow_dir[3]
    )

    return (
        circle_trace = (
            type = "scatter3d",
            mode = "lines",
            x = circle_x,
            y = circle_y,
            z = circle_z,
            hoverinfo = "skip",
            line = (width = 4, color = "rgba(17, 17, 17, 0.45)"),
            showlegend = false
        ),
        tick_trace = (
            type = "scatter3d",
            mode = "lines",
            x = tick_x,
            y = tick_y,
            z = tick_z,
            hoverinfo = "skip",
            line = (width = 8, color = "#111111"),
            showlegend = false
        ),
        arrow_trace = (
            type = "scatter3d",
            mode = "lines",
            x = [center[1], tip.x],
            y = [center[2], tip.y],
            z = [center[3], tip.z],
            hoverinfo = "skip",
            line = (width = 8, color = "#111111"),
            showlegend = false
        ),
        head_trace = (
            type = "cone",
            x = [tip.x],
            y = [tip.y],
            z = [tip.z],
            u = [ρ * arrow_dir[1]],
            v = [ρ * arrow_dir[2]],
            w = [ρ * arrow_dir[3]],
            anchor = "tip",
            sizemode = "absolute",
            sizeref = 0.28,
            colorscale = (((0.0, "#111111"), (1.0, "#111111"))),
            showscale = false,
            hoverinfo = "skip",
            showlegend = false
        ),
        linear = rep.linear
    )
end

"""
    render_poincare_sphere(rep; trail = nothing, title = "Poincare Sphere")

Build a renderable specification for a Poincare sphere with a graphical state vector.

`rep` should be the output of [`poincare_vector_representation`](fiber-path-plot.jl).
The returned named tuple contains the sphere surface, guide circles, the vector, the point,
and any optional trail data, ready to be consumed by the HTML/Plotly renderer.
"""
function render_poincare_sphere(
    rep::NamedTuple;
    trail::Union{Nothing, NamedTuple} = nothing
)
    sphere_u = [Float64[] for _ in 1:31]
    sphere_v = [Float64[] for _ in 1:31]
    sphere_w = [Float64[] for _ in 1:31]
    for i in 0:30
        θ = π * i / 30
        for j in 0:60
            ϕ = 2π * j / 60
            push!(sphere_u[i + 1], sin(θ) * cos(ϕ))
            push!(sphere_v[i + 1], sin(θ) * sin(ϕ))
            push!(sphere_w[i + 1], cos(θ))
        end
    end

    circle_t = collect(range(0.0, 2π, length = 181))
    eq_xy = (x = cos.(circle_t), y = sin.(circle_t), z = zeros(length(circle_t)))
    eq_xz = (x = cos.(circle_t), y = zeros(length(circle_t)), z = sin.(circle_t))
    eq_yz = (x = zeros(length(circle_t)), y = cos.(circle_t), z = sin.(circle_t))

    point = rep.sphere
    vector = (
        x = [0.0, point.x],
        y = [0.0, point.y],
        z = [0.0, point.z]
    )
    labels = (
        x = [1.12, -1.12, 0.0, 0.0],
        y = [0.0, 0.0, 0.0, 0.0],
        z = [0.0, 0.0, 1.12, -1.37],
        text = ["H", "V", "LCP", "RCP"]
    )
    tick_inner = 1.00
    tick_outer = 1.10
    tick_x = Float64[
        tick_inner, tick_outer, NaN,
        -tick_inner, -tick_outer, NaN,
        0.0, 0.0, NaN,
        0.0, 0.0, NaN
    ]
    tick_y = Float64[
        0.0, 0.0, NaN,
        0.0, 0.0, NaN,
        0.0, 0.0, NaN,
        0.0, 0.0, NaN
    ]
    tick_z = Float64[
        0.0, 0.0, NaN,
        0.0, 0.0, NaN,
        tick_inner, tick_outer, NaN,
        -tick_inner, -tick_outer, NaN
    ]
    return (
        title = "Poincare Sphere",
        surface = (x = sphere_u, y = sphere_v, z = sphere_w),
        circles = (xy = eq_xy, xz = eq_xz, yz = eq_yz),
        ticks = (x = tick_x, y = tick_y, z = tick_z),
        vector = vector,
        point = point,
        trail = trail,
        labels = labels
    )
end

# ----------------------------
# visualize the trajectory of a Fiber by sampling the centerline
# geometry and polarization evolution
# ----------------------------
module PlotRuntime
    using LinearAlgebra
    # Bring in runtime API from sibling Bifrost submodules. PlotRuntime lives
    # at depth Bifrost.Plots.PlotRuntime, so `...FiberPath` (= Bifrost.FiberPath)
    # is two parents up. In-file helpers (rotate_about_axis, frenet_serret_frame,
    # poincare_vector_representation, render_pol_circle, render_poincare_sphere)
    # live in `Bifrost.Plots`, reachable as `..Plots`.
    using ...FiberPath: Fiber, bend_geometry, twist_rate, fiber_path
    using ...PathIntegral: generator_K, generator_Kω,
                           exp_sensitivity_midpoint_step, output_dgd
    using ...PathGeometry: frame
    using ..Plots: rotate_about_axis, frenet_serret_frame,
                   poincare_vector_representation,
                   render_pol_circle, render_poincare_sphere

    function js_real(x::Real)
        xf = Float64(x)
        if isnan(xf)
            return "NaN"
        elseif xf == Inf
            return "Infinity"
        elseif xf == -Inf
            return "-Infinity"
        else
            return string(xf)
        end
    end

    js_array(xs::AbstractVector{<:Real}) = "[" * join(js_real.(xs), ", ") * "]"
    js_array2d(xss::AbstractVector{<:AbstractVector{<:Real}}) = "[" * join(js_array.(xss), ", ") * "]"

    function js_string_array(xs::AbstractVector{<:AbstractString})
        escaped = replace.(xs, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
        return "[" * join(["\"" * x * "\"" for x in escaped], ", ") * "]"
    end

    function propagate_state_segment(K, jumps, state, a::Float64, b::Float64)
        if b <= a
            return state
        end

        z = a
        ψ = state
        jump_points = sort([Float64(ζ) for ζ in keys(jumps) if a < ζ <= b])

        for ζ in jump_points
            h1 = ζ - z
            if h1 > 0
                ψ = exp(h1 * K(z + 0.5h1)) * ψ
            end
            ψ = jumps[ζ] * ψ
            z = ζ
        end

        h2 = b - z
        if h2 > 0
            ψ = exp(h2 * K(z + 0.5h2)) * ψ
        end

        return ψ
    end

    function propagate_sensitivity_segment(
        K,
        Kω,
        jumps,
        jump_omegas,
        J::Matrix{ComplexF64},
        G::Matrix{ComplexF64},
        a::Float64,
        b::Float64
    )
        if b <= a
            return J, G
        end

        s = a
        Jseg = J
        Gseg = G
        jump_points = sort([Float64(ζ) for ζ in keys(jumps) if a < ζ <= b])

        for ζ in jump_points
            h1 = ζ - s
            if h1 > 0
                Jseg, Gseg = exp_sensitivity_midpoint_step(K, Kω, s, h1, Jseg, Gseg)
            end
            J_pre = Jseg
            G_pre = Gseg
            J_jump = jumps[ζ]
            G_jump = haskey(jump_omegas, ζ) ? jump_omegas[ζ] : zeros(ComplexF64, 2, 2)
            Jseg = J_jump * J_pre
            Gseg = G_jump * J_pre + J_jump * G_pre
            s = ζ
        end

        h2 = b - s
        if h2 > 0
            Jseg, Gseg = exp_sensitivity_midpoint_step(K, Kω, s, h2, Jseg, Gseg)
        end

        return Jseg, Gseg
    end

    """
        sample_fiber_centerline(f::Fiber, s1, s2; n = 1001)

    Reconstruct a fiber centerline from a [`Fiber`](path-integral.jl)
    on the interval `[s1, s2]`, treating the input coordinate as arc length.

    The bend field defines the local curvature vector

    `κ(s) = (cos(θ_b(s))/Rb(s), sin(θ_b(s))/Rb(s), 0)`

    in the lab frame, and the centerline is recovered by integrating the tangent vector.
    `dtwist` is included in the returned metadata in rad/m but does not affect the centerline
    geometry.

    Returns a named tuple containing `s`, `x`, `y`, `zc`, `tx`, `ty`, `tz`, Frenet-Serret
    `nx`, `ny`, `nz`, `bx`, `by`, `bz`, the transported frame `e1x`, `e1y`, `e1z`, `e2x`,
    `e2y`, `e2z`, `Rb`, `theta_b`, `dtwist`, `kx`, `ky`, and `k2`.
    """
    function sample_fiber_centerline(f::Fiber, s1::Real, s2::Real; n::Int = 1001)
        @assert s2 > s1 "Require s2 > s1"
        @assert n >= 2 "Require at least two sample points"

        path_geometry = fiber_path(f)
        if !isnothing(path_geometry)
            ss = collect(range(Float64(s1), Float64(s2), length = n))
            frames = [frame(path_geometry, s) for s in ss]

            x = [Float64(fr.position[1]) for fr in frames]
            y = [Float64(fr.position[2]) for fr in frames]
            zc = [Float64(fr.position[3]) for fr in frames]
            tx = [Float64(fr.tangent[1]) for fr in frames]
            ty = [Float64(fr.tangent[2]) for fr in frames]
            tz = [Float64(fr.tangent[3]) for fr in frames]
            nx = [Float64(fr.normal[1]) for fr in frames]
            ny = [Float64(fr.normal[2]) for fr in frames]
            nz = [Float64(fr.normal[3]) for fr in frames]
            bx = [Float64(fr.binormal[1]) for fr in frames]
            by = [Float64(fr.binormal[2]) for fr in frames]
            bz = [Float64(fr.binormal[3]) for fr in frames]
            e1x = copy(nx)
            e1y = copy(ny)
            e1z = copy(nz)
            e2x = copy(bx)
            e2y = copy(by)
            e2z = copy(bz)
            kx = [Float64(fr.curvature) for fr in frames]
            ky = zeros(Float64, n)
            k2 = [κ^2 for κ in kx]
            Rb = [iszero(fr.curvature) ? Inf : Float64(inv(fr.curvature)) for fr in frames]
            theta_b = zeros(Float64, n)
            dtwist = [Float64(fr.geometric_torsion + fr.material_twist) for fr in frames]

            return (; s = ss, x, y, zc, tx, ty, tz, e1x, e1y, e1z, e2x, e2y, e2z,
                      Rb, theta_b, dtwist, kx, ky, k2, nx, ny, nz, bx, by, bz)
        end

        ss = collect(range(Float64(s1), Float64(s2), length = n))
        x = zeros(Float64, n)
        y = zeros(Float64, n)
        zc = zeros(Float64, n)
        Rb = Vector{Float64}(undef, n)
        theta_b = Vector{Float64}(undef, n)
        dtwist = Vector{Float64}(undef, n)
        kx = zeros(Float64, n)
        ky = zeros(Float64, n)
        k2 = zeros(Float64, n)
        tx = zeros(Float64, n)
        ty = zeros(Float64, n)
        tz = ones(Float64, n)
        e1x = ones(Float64, n)
        e1y = zeros(Float64, n)
        e1z = zeros(Float64, n)
        e2x = zeros(Float64, n)
        e2y = ones(Float64, n)
        e2z = zeros(Float64, n)
        nx = zeros(Float64, n)
        ny = zeros(Float64, n)
        nz = zeros(Float64, n)
        bx = zeros(Float64, n)
        by = zeros(Float64, n)
        bz = zeros(Float64, n)

        T = [0.0, 0.0, 1.0]
        r = [0.0, 0.0, 0.0]
        E1 = [1.0, 0.0, 0.0]
        E2 = [0.0, 1.0, 0.0]

        tx[1], ty[1], tz[1] = T
        e1x[1], e1y[1], e1z[1] = E1
        e2x[1], e2y[1], e2z[1] = E2

        for i in eachindex(ss)
            si = ss[i]
            bend = bend_geometry(f, si)
            R = Float64(bend.Rb)
            θ = Float64(bend.theta_b)
            τ = Float64(twist_rate(f, si))

            Rb[i] = R
            theta_b[i] = θ
            dtwist[i] = τ

            kx[i] = bend.kx
            ky[i] = bend.ky
            k2[i] = bend.k2

            if i < length(ss)
                ds = ss[i + 1] - ss[i]
                Ω = [-ky[i], kx[i], 0.0]
                ω = norm(Ω)

                if ω == 0.0
                    Tnext = T
                    E1next = E1
                    E2next = E2
                else
                    u = Ω / ω
                    α = ω * ds
                    Tnext = rotate_about_axis(T, u, α)
                    E1next = rotate_about_axis(E1, u, α)
                    E2next = rotate_about_axis(E2, u, α)
                end

                E1next .-= dot(E1next, Tnext) .* Tnext
                E1next ./= norm(E1next)
                E2next = cross(Tnext, E1next)
                E2next ./= norm(E2next)

                Tmid = 0.5 .* (T .+ Tnext)
                Tmid ./= norm(Tmid)
                r .+= ds .* Tmid
                T = Tnext / norm(Tnext)
                E1 = E1next
                E2 = E2next

                x[i + 1], y[i + 1], zc[i + 1] = r
                tx[i + 1], ty[i + 1], tz[i + 1] = T
                e1x[i + 1], e1y[i + 1], e1z[i + 1] = E1
                e2x[i + 1], e2y[i + 1], e2z[i + 1] = E2
            end
        end

        path = (; s = ss, x, y, zc, tx, ty, tz, e1x, e1y, e1z, e2x, e2y, e2z, Rb, theta_b, dtwist, kx, ky, k2)
        for i in eachindex(ss)
            fs = frenet_serret_frame(path, ss[i])
            nx[i], ny[i], nz[i] = fs.normal
            bx[i], by[i], bz[i] = fs.binormal
        end

        return (; path..., nx, ny, nz, bx, by, bz)
    end

    sample_fiber_input = sample_fiber_centerline

    function sample_polarization_trajectory(
        f::Fiber,
        s1::Real,
        s2::Real;
        λ_m::Real,
        n::Int = 1001,
        input_state::AbstractVector{ComplexF64} = ComplexF64[1.0 + 0.0im, 0.0 + 0.0im],
        jumps::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}()
    )
        @assert length(input_state) == 2 "input_state must have length 2"

        ss = collect(range(Float64(s1), Float64(s2), length = n))
        ψs = Vector{Vector{ComplexF64}}(undef, n)
        s1 = zeros(Float64, n)
        s2 = zeros(Float64, n)
        s3 = zeros(Float64, n)
        linear_angle_rad = zeros(Float64, n)
        linear_radius = zeros(Float64, n)

        K = generator_K(f, λ_m)
        ψ = ComplexF64[input_state[1], input_state[2]]
        ψ ./= norm(ψ)

        ψs[1] = copy(ψ)
        rep = poincare_vector_representation(ψ)
        s1[1], s2[1], s3[1] = rep.sphere
        linear_angle_rad[1] = rep.linear.angle_rad
        linear_radius[1] = rep.linear.radius

        for i in 1:n-1
            ψ = propagate_state_segment(K, jumps, ψ, ss[i], ss[i + 1])
            ψ ./= norm(ψ)
            ψs[i + 1] = copy(ψ)
            rep = poincare_vector_representation(ψ)
            s1[i + 1], s2[i + 1], s3[i + 1] = rep.sphere
            linear_angle_rad[i + 1] = rep.linear.angle_rad
            linear_radius[i + 1] = rep.linear.radius
        end

        return (; s = ss, state = ψs, s1, s2, s3, linear_angle_rad, linear_radius)
    end

    function sample_dgd_trajectory(
        f::Fiber,
        s1::Real,
        s2::Real;
        λ_m::Real,
        n::Int = 1001,
        jumps::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}(),
        jump_omegas::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}()
    )
        ss = collect(range(Float64(s1), Float64(s2), length = n))
        dgds = zeros(Float64, n)

        K = generator_K(f, λ_m)
        Kω = generator_Kω(f, λ_m)
        J = Matrix{ComplexF64}(I, 2, 2)
        G = zeros(ComplexF64, 2, 2)
        dgds[1] = output_dgd(J, G)

        for i in 1:n-1
            J, G = propagate_sensitivity_segment(K, Kω, jumps, jump_omegas, J, G, ss[i], ss[i + 1])
            dgds[i + 1] = output_dgd(J, G)
        end

        return (; s = ss, dgd = dgds)
    end

    """
        write_fiber_input_plot3d(
            f::Fiber,
            s1,
            s2;
            n = 1001,
            output = "fiberinput_3d.html",
            title = "Fiber 3D centerline",
            input_state = ComplexF64[1.0 + 0.0im, 0.0 + 0.0im],
            jumps = Dict{Float64, Matrix{ComplexF64}}()
        )

    Write an interactive Plotly HTML file showing the sampled `Fiber` over `s ∈ [s1, s2]`.

    The 3D curve is plotted as the reconstructed fiber centerline `(x(s), y(s), zc(s))`.
    Plotly provides mouse-based rotation, zooming, and panning in the generated HTML.

    For interactive geometry of an analytic [`PathGeometry.Path`](path-geometry-plot.jl) built
    from `path-geometry.jl` (no fiber or polarization layers), see [`write_path_geometry_plot3d`](path-geometry-plot.jl).

    Returns the output path.
    """
    function write_fiber_input_plot3d(
        f::Fiber,
        s1::Real,
        s2::Real;
        λ_m::Real,
        n::Int = 1001,
        output::AbstractString = "fiberinput_3d.html",
        title::AbstractString = "Fiber 3D centerline",
        input_state::AbstractVector{ComplexF64} = ComplexF64[1.0 + 0.0im, 0.0 + 0.0im],
        jumps::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}(),
        jump_omegas::Dict{Float64, Matrix{ComplexF64}} = Dict{Float64, Matrix{ComplexF64}}()
    )
        samples = sample_fiber_centerline(f, s1, s2; n = n)
        pol = sample_polarization_trajectory(
            f,
            s1,
            s2;
            λ_m = λ_m,
            n = n,
            input_state = ComplexF64[input_state[1], input_state[2]],
            jumps = jumps
        )
        dgd = sample_dgd_trajectory(
            f,
            s1,
            s2;
            λ_m = λ_m,
            n = n,
            jumps = jumps,
            jump_omegas = jump_omegas
        )
        pol_circle = render_pol_circle(samples, samples.s[1], poincare_vector_representation(pol.state[1]))
        sphere_spec = render_poincare_sphere(
            poincare_vector_representation(pol.state[1]);
            trail = (x = [pol.s1[1]], y = [pol.s2[1]], z = [pol.s3[1]], color = [samples.s[1]])
        )

        hover_text = [
            "s=$(samples.s[i]) m<br>" *
            "x=$(samples.x[i]) m<br>" *
            "y=$(samples.y[i]) m<br>" *
            "z=$(samples.zc[i]) m<br>" *
            "Rb=$(samples.Rb[i]) m<br>" *
            "theta_b=$(samples.theta_b[i]) rad<br>" *
            "dtwist=$(samples.dtwist[i]) rad/m<br>" *
            "DGD=$(dgd.dgd[i]) s<br>" *
            "S1=$(pol.s1[i])<br>" *
            "S2=$(pol.s2[i])<br>" *
            "S3=$(pol.s3[i])"
            for i in eachindex(samples.s)
        ]

        html = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
          <meta charset="utf-8" />
          <meta name="viewport" content="width=device-width, initial-scale=1" />
          <title>$title</title>
          <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
          <style>
            html, body {
              width: 100%;
              height: 100%;
              margin: 0;
              font-family: sans-serif;
            }
            #viewer {
              position: relative;
              width: 100%;
              height: 100%;
              overflow: hidden;
            }
            #plot {
              width: 100%;
              height: 100%;
            }
            #inset {
              position: absolute;
              top: 16px;
              right: 16px;
              width: min(28vw, 340px);
              height: min(28vw, 340px);
              min-width: 220px;
              min-height: 220px;
              border: 1px solid rgba(0, 0, 0, 0.15);
              background: rgba(255, 255, 255, 0.88);
              box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
            }
            #status {
              position: absolute;
              top: calc(16px + min(28vw, 340px) + 10px);
              right: 16px;
              width: min(28vw, 340px);
              min-width: 220px;
              padding: 10px 12px;
              box-sizing: border-box;
              border: 1px solid rgba(0, 0, 0, 0.15);
              background: rgba(255, 255, 255, 0.9);
              box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
              font-size: 13px;
              line-height: 1.4;
              white-space: pre-line;
            }
            #scrub-help {
              position: absolute;
              left: 16px;
              bottom: 18px;
              padding: 8px 10px;
              font-size: 13px;
              background: rgba(255, 255, 255, 0.86);
              border: 1px solid rgba(0, 0, 0, 0.12);
            }
          </style>
        </head>
        <body>
          <div id="viewer">
            <div id="plot"></div>
            <div id="inset"></div>
            <div id="status"></div>
            <div id="scrub-help">Move the mouse left-to-right over the plot to scrub the fiber cursor and polarization state.</div>
          </div>
          <script>
            const xs = $(js_array(samples.x));
            const ys = $(js_array(samples.y));
            const zs = $(js_array(samples.zc));
            const ss = $(js_array(samples.s));
            const s1 = $(js_array(pol.s1));
            const s2 = $(js_array(pol.s2));
            const s3 = $(js_array(pol.s3));
            const nx = $(js_array(samples.nx));
            const ny = $(js_array(samples.ny));
            const nz = $(js_array(samples.nz));
            const bx = $(js_array(samples.bx));
            const by = $(js_array(samples.by));
            const bz = $(js_array(samples.bz));
            const rb = $(js_array(samples.Rb));
            const theta = $(js_array(samples.theta_b));
            const dtwist = $(js_array(samples.dtwist));
            const dgd = $(js_array(dgd.dgd));
            const hoverText = $(js_string_array(hover_text));
            const linearAngle = $(js_array(pol.linear_angle_rad));
            const arrowX = new Array(xs.length).fill(0);
            const arrowY = new Array(xs.length).fill(0);
            const arrowZ = new Array(xs.length).fill(0);
            const arrowLen = new Array(xs.length).fill(0);
            const circleRadX = new Array(xs.length).fill(0);
            const circleRadY = new Array(xs.length).fill(0);
            const circleRadZ = new Array(xs.length).fill(0);

            function rawArrowVector(index) {
              const linMag = Math.min(1.0, Math.hypot(s1[index], s2[index]));
              if (linMag < 1e-12) {
                return { x: 0, y: 0, z: 0, length: 0 };
              }
              const psi = linearAngle[index];
              return {
                x: linMag * (Math.cos(psi) * nx[index] + Math.sin(psi) * bx[index]),
                y: linMag * (Math.cos(psi) * ny[index] + Math.sin(psi) * by[index]),
                z: linMag * (Math.cos(psi) * nz[index] + Math.sin(psi) * bz[index]),
                length: linMag
              };
            }

            for (let i = 0; i < xs.length; i += 1) {
              const a = rawArrowVector(i);
              arrowX[i] = a.x;
              arrowY[i] = a.y;
              arrowZ[i] = a.z;
              arrowLen[i] = a.length;
              circleRadX[i] = Math.hypot(nx[i], bx[i]);
              circleRadY[i] = Math.hypot(ny[i], by[i]);
              circleRadZ[i] = Math.hypot(nz[i], bz[i]);
            }

            for (let i = 1; i < xs.length; i += 1) {
              const dot =
                arrowX[i] * arrowX[i - 1] +
                arrowY[i] * arrowY[i - 1] +
                arrowZ[i] * arrowZ[i - 1];
              if (dot < 0) {
                arrowX[i] = -arrowX[i];
                arrowY[i] = -arrowY[i];
                arrowZ[i] = -arrowZ[i];
              }
            }

            function arrowVector(index) {
              return {
                x: arrowX[index],
                y: arrowY[index],
                z: arrowZ[index],
                length: arrowLen[index],
                angleRad: linearAngle[index]
              };
            }

            function polCircle(index, npts = 121) {
              const cx = [];
              const cy = [];
              const cz = [];
              const tickX = [];
              const tickY = [];
              const tickZ = [];
              for (let i = 0; i < npts; i += 1) {
                const t = 2 * Math.PI * i / (npts - 1);
                cx.push(xs[index] + Math.cos(t) * nx[index] + Math.sin(t) * bx[index]);
                cy.push(ys[index] + Math.cos(t) * ny[index] + Math.sin(t) * by[index]);
                cz.push(zs[index] + Math.cos(t) * nz[index] + Math.sin(t) * bz[index]);
              }
              const tickAngles = [0, 0.5 * Math.PI, Math.PI, 1.5 * Math.PI];
              const tickInner = 0.82;
              for (const phi of tickAngles) {
                const rimX = Math.cos(phi) * nx[index] + Math.sin(phi) * bx[index];
                const rimY = Math.cos(phi) * ny[index] + Math.sin(phi) * by[index];
                const rimZ = Math.cos(phi) * nz[index] + Math.sin(phi) * bz[index];
                tickX.push(xs[index] + tickInner * rimX, xs[index] + rimX, NaN);
                tickY.push(ys[index] + tickInner * rimY, ys[index] + rimY, NaN);
                tickZ.push(zs[index] + tickInner * rimZ, zs[index] + rimZ, NaN);
              }
              return { x: cx, y: cy, z: cz, tickX, tickY, tickZ };
            }

            function computeBounds(values, offsets) {
              let vmin = Infinity;
              let vmax = -Infinity;
              for (let i = 0; i < values.length; i += 1) {
                const a = values[i];
                const b = values[i] + offsets[i];
                if (a < vmin) vmin = a;
                if (a > vmax) vmax = a;
                if (b < vmin) vmin = b;
                if (b > vmax) vmax = b;
              }
              return [vmin, vmax];
            }

            const headMargin = 1.01;
            const xRangeRaw = computeBounds(xs, circleRadX.map((r, i) => Math.max(r, Math.abs(arrowX[i]))));
            const yRangeRaw = computeBounds(ys, circleRadY.map((r, i) => Math.max(r, Math.abs(arrowY[i]))));
            const zRangeRaw = computeBounds(zs, circleRadZ.map((r, i) => Math.max(r, Math.abs(arrowZ[i]))));
            const xPad = Math.max(0.08 * (xRangeRaw[1] - xRangeRaw[0]), headMargin);
            const yPad = Math.max(0.08 * (yRangeRaw[1] - yRangeRaw[0]), headMargin);
            const zPad = Math.max(0.08 * (zRangeRaw[1] - zRangeRaw[0]), headMargin);
            const xRange = [xRangeRaw[0] - xPad, xRangeRaw[1] + xPad];
            const yRange = [yRangeRaw[0] - yPad, yRangeRaw[1] + yPad];
            const zRange = [zRangeRaw[0] - zPad, zRangeRaw[1] + zPad];

            const circle0 = polCircle(0);

            const fiberTrace = {
              type: "scatter3d",
              mode: "lines+markers",
              x: xs,
              y: ys,
              z: zs,
              hoverinfo: "skip",
              line: {
                width: 4,
                color: "rgba(30, 30, 30, 0.22)"
              },
              marker: {
                size: 3.5,
                color: ss,
                colorscale: "Turbo",
                cmin: ss[0],
                cmax: ss[ss.length - 1],
                showscale: false
              }
            };

            const startTrace = {
              type: "scatter3d",
              mode: "markers",
              x: [$(samples.x[1])],
              y: [$(samples.y[1])],
              z: [$(samples.zc[1])],
              hovertemplate: "Start<br>x=%{x} m<br>y=%{y} m<br>z=%{z} m<extra></extra>",
              marker: {
                size: 7,
                color: "#2ca02c",
                symbol: "circle"
              },
              name: "start"
            };

            const endTrace = {
              type: "scatter3d",
              mode: "markers",
              x: [$(samples.x[end])],
              y: [$(samples.y[end])],
              z: [$(samples.zc[end])],
              hovertemplate: "End<br>x=%{x} m<br>y=%{y} m<br>z=%{z} m<extra></extra>",
              marker: {
                size: 7,
                color: "#d62728",
                symbol: "circle"
              },
              name: "end"
            };

            const cursorTrace = {
              type: "scatter3d",
              mode: "markers",
              x: [xs[0]],
              y: [ys[0]],
              z: [zs[0]],
              hoverinfo: "skip",
              marker: {
                size: 6,
                color: "#111111",
                symbol: "circle"
              },
              name: "cursor"
            };

            const polCircleTrace = {
              type: $(repr(pol_circle.circle_trace.type)),
              mode: $(repr(pol_circle.circle_trace.mode)),
              x: $(js_array(pol_circle.circle_trace.x)),
              y: $(js_array(pol_circle.circle_trace.y)),
              z: $(js_array(pol_circle.circle_trace.z)),
              hoverinfo: $(repr(pol_circle.circle_trace.hoverinfo)),
              line: {
                width: $(pol_circle.circle_trace.line.width),
                color: $(repr(pol_circle.circle_trace.line.color))
              },
              showlegend: $(pol_circle.circle_trace.showlegend)
            };

            const polTickTrace = {
              type: $(repr(pol_circle.tick_trace.type)),
              mode: $(repr(pol_circle.tick_trace.mode)),
              x: $(js_array(pol_circle.tick_trace.x)),
              y: $(js_array(pol_circle.tick_trace.y)),
              z: $(js_array(pol_circle.tick_trace.z)),
              hoverinfo: $(repr(pol_circle.tick_trace.hoverinfo)),
              line: {
                width: $(pol_circle.tick_trace.line.width),
                color: $(repr(pol_circle.tick_trace.line.color))
              },
              showlegend: $(pol_circle.tick_trace.showlegend)
            };

            const polArrowTrace = {
              type: $(repr(pol_circle.arrow_trace.type)),
              mode: $(repr(pol_circle.arrow_trace.mode)),
              x: $(js_array(pol_circle.arrow_trace.x)),
              y: $(js_array(pol_circle.arrow_trace.y)),
              z: $(js_array(pol_circle.arrow_trace.z)),
              hoverinfo: $(repr(pol_circle.arrow_trace.hoverinfo)),
              line: {
                width: $(pol_circle.arrow_trace.line.width),
                color: $(repr(pol_circle.arrow_trace.line.color))
              },
              showlegend: $(pol_circle.arrow_trace.showlegend)
            };

            const polArrowHead = {
              type: $(repr(pol_circle.head_trace.type)),
              x: $(js_array(pol_circle.head_trace.x)),
              y: $(js_array(pol_circle.head_trace.y)),
              z: $(js_array(pol_circle.head_trace.z)),
              u: $(js_array(pol_circle.head_trace.u)),
              v: $(js_array(pol_circle.head_trace.v)),
              w: $(js_array(pol_circle.head_trace.w)),
              anchor: $(repr(pol_circle.head_trace.anchor)),
              sizemode: $(repr(pol_circle.head_trace.sizemode)),
              sizeref: $(pol_circle.head_trace.sizeref),
              colorscale: [[0, "#111111"], [1, "#111111"]],
              showscale: $(pol_circle.head_trace.showscale),
              hoverinfo: $(repr(pol_circle.head_trace.hoverinfo)),
              showlegend: $(pol_circle.head_trace.showlegend)
            };

            const layout = {
              title: $((repr(title))),
              showlegend: true,
              legend: { x: 0.02, y: 0.98 },
              margin: { l: 0, r: 0, b: 0, t: 48 },
              scene: {
                xaxis: { title: "x (m)", range: xRange, autorange: false, showspikes: false },
                yaxis: { title: "y (m)", range: yRange, autorange: false, showspikes: false },
                zaxis: { title: "z (m)", range: zRange, autorange: false, showspikes: false },
                camera: {
                  eye: { x: 1.6, y: 1.6, z: 0.9 }
                }
              }
            };

            Plotly.newPlot("plot", [fiberTrace, startTrace, endTrace, cursorTrace, polCircleTrace, polTickTrace, polArrowTrace, polArrowHead], layout, {
              responsive: true,
              scrollZoom: true,
              displaylogo: false
            });

            const sphereSurface = {
              type: "surface",
              x: $(js_array2d(sphere_spec.surface.x)),
              y: $(js_array2d(sphere_spec.surface.y)),
              z: $(js_array2d(sphere_spec.surface.z)),
              opacity: 0.18,
              showscale: false,
              colorscale: [[0, "#cfd8dc"], [1, "#cfd8dc"]],
              hoverinfo: "skip"
            };

            const sphereLines = [
              { x: $(js_array(sphere_spec.circles.xy.x)), y: $(js_array(sphere_spec.circles.xy.y)), z: $(js_array(sphere_spec.circles.xy.z)) },
              { x: $(js_array(sphere_spec.circles.xz.x)), y: $(js_array(sphere_spec.circles.xz.y)), z: $(js_array(sphere_spec.circles.xz.z)) },
              { x: $(js_array(sphere_spec.circles.yz.x)), y: $(js_array(sphere_spec.circles.yz.y)), z: $(js_array(sphere_spec.circles.yz.z)) }
            ].map(curve => ({
              type: "scatter3d",
              mode: "lines",
              x: curve.x,
              y: curve.y,
              z: curve.z,
              line: { color: "#90a4ae", width: 4 },
              hoverinfo: "skip",
              showlegend: false
            }));

            const sphereTicks = {
              type: "scatter3d",
              mode: "lines",
              x: $(js_array(sphere_spec.ticks.x)),
              y: $(js_array(sphere_spec.ticks.y)),
              z: $(js_array(sphere_spec.ticks.z)),
              line: { color: "#111111", width: 6 },
              hoverinfo: "skip",
              showlegend: false
            };

            const stateTrail = {
              type: "scatter3d",
              mode: "markers",
              x: $(js_array(sphere_spec.trail.x)),
              y: $(js_array(sphere_spec.trail.y)),
              z: $(js_array(sphere_spec.trail.z)),
              marker: {
                size: 3.5,
                color: $(js_array(sphere_spec.trail.color)),
                colorscale: "Turbo",
                cmin: ss[0],
                cmax: ss[ss.length - 1],
                showscale: false
              },
              hoverinfo: "skip",
              showlegend: false
            };

            const stateVector = {
              type: "scatter3d",
              mode: "lines",
              x: $(js_array(sphere_spec.vector.x)),
              y: $(js_array(sphere_spec.vector.y)),
              z: $(js_array(sphere_spec.vector.z)),
              line: { color: "#111111", width: 6 },
              hoverinfo: "skip",
              showlegend: false
            };

            const statePoint = {
              type: "scatter3d",
              mode: "markers",
              x: [$(sphere_spec.point.x)],
              y: [$(sphere_spec.point.y)],
              z: [$(sphere_spec.point.z)],
              hovertemplate: "Poincare state<br>S1=%{x:.4f}<br>S2=%{y:.4f}<br>S3=%{z:.4f}<extra></extra>",
              marker: {
                size: 6,
                color: [ss[0]],
                colorscale: "Turbo",
                cmin: ss[0],
                cmax: ss[ss.length - 1],
                showscale: false
              },
              name: "state"
            };

            const sphereLabels = {
              type: "scatter3d",
              mode: "text",
              x: $(js_array(sphere_spec.labels.x)),
              y: $(js_array(sphere_spec.labels.y)),
              z: $(js_array(sphere_spec.labels.z)),
              text: $(js_string_array(sphere_spec.labels.text)),
              textfont: { size: 14, color: "#111111", family: "sans-serif" },
              hoverinfo: "skip",
              showlegend: false
            };

            const insetLayout = {
              margin: { l: 18, r: 18, b: 18, t: 34 },
              title: { text: "Poincare Sphere", x: 0.5, xanchor: "center", font: { size: 14 } },
              showlegend: false,
              scene: {
                xaxis: { title: { text: "S1", font: { size: 12 } }, range: [-1.18, 1.18] },
                yaxis: { title: { text: "S2", font: { size: 12 } }, range: [-1.18, 1.18] },
                zaxis: { title: { text: "S3", font: { size: 12 } }, range: [-1.18, 1.18] },
                aspectmode: "cube",
                camera: { eye: { x: 1.35, y: 1.35, z: 1.0 } }
              }
            };

            Plotly.newPlot(
              "inset",
              [sphereSurface, ...sphereLines, sphereTicks, stateTrail, stateVector, statePoint, sphereLabels],
              insetLayout,
              { responsive: true, displaylogo: false, scrollZoom: true }
            );

            const statusBox = document.getElementById("status");
            let activeIndex = 0;

            function formatStatus(index) {
              const arrow = arrowVector(index);
              return [
                "Status at fiber path point s: \\\\n",
                "s = " + ss[index].toFixed(4) + " m",
                "x = " + xs[index].toFixed(4) + " m",
                "y = " + ys[index].toFixed(4) + " m",
                "z = " + zs[index].toFixed(4) + " m",
                "Rb = " + (Number.isFinite(rb[index]) ? rb[index].toFixed(4) : "Inf") + " m",
                "theta_b = " + theta[index].toFixed(4) + " rad",
                "dtwist = " + dtwist[index].toFixed(4) + " rad/m",
                "DGD = " + dgd[index].toExponential(6) + " s",
                "S1 = " + s1[index].toFixed(4),
                "S2 = " + s2[index].toFixed(4),
                "S3 = " + s3[index].toFixed(4),
                "linear angle = " + arrow.angleRad.toFixed(4) + " rad",
                "arrow length = " + arrow.length.toFixed(4)
              ].join("\\n");
            }

            function updateCursor(index) {
              activeIndex = Math.max(0, Math.min(xs.length - 1, index));
              const circle = polCircle(activeIndex);
              const arrow = arrowVector(activeIndex);
              Plotly.restyle("plot", {
                x: [[xs[activeIndex]]],
                y: [[ys[activeIndex]]],
                z: [[zs[activeIndex]]]
              }, [3]);
              Plotly.restyle("plot", {
                x: [circle.x],
                y: [circle.y],
                z: [circle.z]
              }, [4]);
              Plotly.restyle("plot", {
                x: [circle.tickX],
                y: [circle.tickY],
                z: [circle.tickZ]
              }, [5]);
              Plotly.restyle("plot", {
                x: [[xs[activeIndex], xs[activeIndex] + arrow.x]],
                y: [[ys[activeIndex], ys[activeIndex] + arrow.y]],
                z: [[zs[activeIndex], zs[activeIndex] + arrow.z]]
              }, [6]);
              Plotly.restyle("plot", {
                x: [[xs[activeIndex] + arrow.x]],
                y: [[ys[activeIndex] + arrow.y]],
                z: [[zs[activeIndex] + arrow.z]],
                u: [[arrow.x]],
                v: [[arrow.y]],
                w: [[arrow.z]]
              }, [7]);
              Plotly.restyle("inset", {
                x: [s1.slice(0, activeIndex + 1)],
                y: [s2.slice(0, activeIndex + 1)],
                z: [s3.slice(0, activeIndex + 1)],
                "marker.color": [ss.slice(0, activeIndex + 1)]
              }, [4]);
              Plotly.restyle("inset", {
                x: [[0, s1[activeIndex]]],
                y: [[0, s2[activeIndex]]],
                z: [[0, s3[activeIndex]]]
              }, [5]);
              Plotly.restyle("inset", {
                x: [[s1[activeIndex]]],
                y: [[s2[activeIndex]]],
                z: [[s3[activeIndex]]],
                "marker.color": [[ss[activeIndex]]]
              }, [6]);
              statusBox.textContent = formatStatus(activeIndex);
            }

            const viewer = document.getElementById("viewer");
            viewer.addEventListener("mousemove", event => {
              if (event.buttons !== 0) {
                return;
              }
              const rect = viewer.getBoundingClientRect();
              const t = Math.max(0, Math.min(1, (event.clientX - rect.left) / rect.width));
              const index = Math.round(t * (xs.length - 1));
              if (index !== activeIndex) {
                updateCursor(index);
              }
            });

            viewer.addEventListener("touchmove", event => {
              if (event.touches.length === 0) {
                return;
              }
              const rect = viewer.getBoundingClientRect();
              const t = Math.max(0, Math.min(1, (event.touches[0].clientX - rect.left) / rect.width));
              const index = Math.round(t * (xs.length - 1));
              if (index !== activeIndex) {
                updateCursor(index);
              }
            }, { passive: true });

            updateCursor(0);
          </script>
        </body>
        </html>
        """

        open(output, "w") do io
            write(io, html)
        end

        return output
    end

end

const sample_fiber_centerline = PlotRuntime.sample_fiber_centerline
const sample_fiber_input = PlotRuntime.sample_fiber_input
const write_fiber_input_plot3d = PlotRuntime.write_fiber_input_plot3d

# ----------------------------------------------------------------------
# Adaptive step-doubling diagnostic
# ----------------------------------------------------------------------

struct StepRecord
    s_start::Float64   # position at the start of this trial
    h::Float64         # step size tried
    err_abs::Float64   # phase-insensitive error between full and two-half steps
    tol::Float64       # tolerance threshold at this step
    accepted::Bool
end

"""
    collect_adaptive_steps(K, s0, s1, J0; kwargs...)

Diagnostic wrapper that mirrors `propagate_interval!` exactly but records every trial step
(accepted and rejected) as a `StepRecord`. The production solver is not modified.

Returns `(J_final, records::Vector{StepRecord})`.
"""
function collect_adaptive_steps(
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
    records = StepRecord[]

    while s < s1
        h = min(h, s1 - s)
        h >= h_min || error("Step size underflow near s=$s")

        J_full     = exp_midpoint_step(K, s, h, J)
        J_half     = exp_midpoint_step(K, s, 0.5h, J)
        J_twohalf  = exp_midpoint_step(K, s + 0.5h, 0.5h, J_half)

        err_abs = scalar_reduce(phase_insensitive_error(J_full, J_twohalf))
        scale = max(
            scalar_reduce(_frobenius_norm(J_twohalf)),
            scalar_reduce(_frobenius_norm(J_full)),
        )
        tol = atol + rtol * scale

        if err_abs <= tol
            push!(records, StepRecord(s, h, err_abs, tol, true))
            s += h
            J  = J_twohalf
            fac = err_abs == 0.0 ? growth_max : clamp(safety * (tol / err_abs)^(1/3), 1.0, growth_max)
            h = min(h * fac, h_max)
        else
            push!(records, StepRecord(s, h, err_abs, tol, false))
            fac = clamp(safety * (tol / err_abs)^(1/3), shrink_min, 0.5)
            h *= fac
            h >= h_min || error("Step size fell below h_min near s=$s")
        end
    end

    return J, records
end

"""
    write_adaptive_steps_plot(records, K_norm_fn, s0, s1; output, title, rtol, atol)

Write an HTML file with two stacked Plotly panels:

- **Top**: accepted step sizes (green) and rejected step sizes (red crosses) vs position.
- **Bottom**: `err_abs / tol` ratio for every trial step, with the acceptance threshold
  drawn as a horizontal dashed line at 1.

`K_norm_fn(s)` should return a scalar representing the local "stiffness" of the generator
(e.g. `opnorm(K(s))`). It is plotted as a shaded background overlay on the top panel so the
viewer can see where small steps coincide with rapid generator variation.

Returns the output path.
"""
function write_adaptive_steps_plot(
    records::Vector{StepRecord},
    K_norm_fn,
    s0::Float64,
    s1::Float64;
    output::AbstractString = "adaptive-steps.html",
    title::AbstractString  = "Adaptive step-doubling diagnostics",
    rtol::Float64 = 1e-9,
    atol::Float64 = 1e-12,
    components = nothing,
)
    acc = filter(r -> r.accepted,  records)
    rej = filter(r -> !r.accepted, records)

    s_acc   = [r.s_start       for r in acc]
    h_acc   = [r.h             for r in acc]
    s_rej   = [r.s_start       for r in rej]
    h_rej   = [r.h             for r in rej]

    s_all   = [r.s_start       for r in records]
    ratio   = [r.err_abs / max(r.tol, 1e-300) for r in records]
    accepted_flag = [r.accepted for r in records]

    bar_x  = Float64[]
    bar_y  = Float64[]
    bar2_x = Float64[]
    bar2_y = Float64[]
    ratios = [r.err_abs / max(r.tol, 1e-300) for r in records]
    i = 1
    while i <= length(records)
        if !records[i].accepted
            s_pos   = records[i].s_start
            h_top   = records[i].h
            r_top   = ratios[i]
            j = i + 1
            while j <= length(records) && !records[j].accepted
                j += 1
            end
            if j <= length(records) && records[j].accepted
                append!(bar_x,  [s_pos, s_pos, NaN])
                append!(bar_y,  [h_top, records[j].h, NaN])
                append!(bar2_x, [s_pos, s_pos, NaN])
                append!(bar2_y, [r_top, ratios[j], NaN])
            end
            i = j
        else
            i += 1
        end
    end

    ns = 400
    s_bg = collect(range(s0, s1, length = ns))
    knorm = [Float64(opnorm(K_norm_fn(si))) for si in s_bg]

    js_vec(v) = "[" * join(string.(Float64.(v)), ", ") * "]"
    js_bool_vec(v) = "[" * join(ifelse.(v, "true", "false"), ", ") * "]"

    has_components = components !== nothing && !isempty(components)
    comp_traces_js = ""
    if has_components
        buf = IOBuffer()
        for (idx, c) in enumerate(components)
            cname, cfn, ccolor = c[1], c[2], c[3]
            cy = [Float64(cfn(si)) for si in s_bg]
            println(buf, "const comp$(idx)_y = $(js_vec(cy));")
            println(buf, "const compTrace$(idx) = {")
            println(buf, "  type: \"scatter\", mode: \"lines\", name: $(repr(cname)),")
            println(buf, "  x: s_bg, y: comp$(idx)_y,")
            println(buf, "  line: { color: $(repr(ccolor)), width: 1.8 },")
            println(buf, "  yaxis: \"y3\"")
            println(buf, "};")
        end
        comp_traces_js = String(take!(buf))
    end
    comp_trace_names = has_components ?
        join(["compTrace$(i)" for i in 1:length(components)], ", ") : ""

    if has_components
        layout_grid       = "{ rows: 3, columns: 1, pattern: \"independent\", roworder: \"top to bottom\" }"
        yaxis_domain      = "[0.55, 1.0]"
        yaxis2_domain     = "[0.26, 0.50]"
        yaxis3_domain     = "[0.0, 0.20]"
        ann_step_y        = "1.0"
        ann_err_y         = "0.49"
        ann_comp          = ",\n    { text: \"Generator components\", xref: \"paper\", yref: \"paper\", x: 0.01, y: 0.20, xanchor: \"left\", yanchor: \"top\", showarrow: false, font: { size: 13 } }"
        xaxis3_block      = ",\n  xaxis3: { title: \"s (arc length)\", domain: [0, 1] },\n  yaxis3: { title: \"coefficient\", domain: $(yaxis3_domain), zeroline: true, zerolinecolor: \"#999\" }"
        threshold_zero_trace = ",\n  xaxis2: { title: \"s (arc length)\", domain: [0, 1], showticklabels: false }"
        xaxis_block       = "{ title: \"s (arc length)\", domain: [0, 1], showticklabels: false }"
        plot_traces_extra = ", " * comp_trace_names
    else
        layout_grid       = "{ rows: 2, columns: 1, pattern: \"independent\", roworder: \"top to bottom\" }"
        yaxis_domain      = "[0.42, 1.0]"
        yaxis2_domain     = "[0.0, 0.38]"
        ann_step_y        = "1.0"
        ann_err_y         = "0.37"
        ann_comp          = ""
        xaxis3_block      = ""
        threshold_zero_trace = ",\n  xaxis2: { title: \"s (arc length)\", domain: [0, 1] }"
        xaxis_block       = "{ title: \"s (arc length)\", domain: [0, 1] }"
        plot_traces_extra = ""
    end

    html = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>$title</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    html, body { margin: 0; width: 100%; height: 100%; font-family: sans-serif; background: #fafafa; }
    #container { width: 100%; height: 100%; }
  </style>
</head>
<body>
<div id="container"></div>
<script>
const s_acc   = $(js_vec(s_acc));
const h_acc   = $(js_vec(h_acc));
const s_rej   = $(js_vec(s_rej));
const h_rej   = $(js_vec(h_rej));
const s_all   = $(js_vec(s_all));
const ratio   = $(js_vec(ratio));
const acc_flag = $(js_bool_vec(accepted_flag));
const s_bg    = $(js_vec(s_bg));
const knorm   = $(js_vec(knorm));
const bar_x   = $(js_vec(bar_x));
const bar_y   = $(js_vec(bar_y));
const bar2_x  = $(js_vec(bar2_x));
const bar2_y  = $(js_vec(bar2_y));

$(comp_traces_js)
const hAll   = h_acc.concat(h_rej);
const hMax   = Math.max(...hAll) * 1.1;
const hMin   = Math.min(...hAll) * 0.85;
const knMax  = Math.max(...knorm) || 1;
const knNorm = knorm.map(k => hMin + (hMax - hMin) * 0.25 * (k / knMax));

const stiffTrace = {
  type: "scatter", mode: "lines", name: "‖K(s)‖ (scaled)",
  x: s_bg, y: knNorm,
  fill: "tozeroy", fillcolor: "rgba(200,220,255,0.35)",
  line: { color: "rgba(100,140,220,0.55)", width: 1.5 },
  yaxis: "y", hoverinfo: "skip", showlegend: true
};

const hAccTrace = {
  type: "scatter", mode: "markers", name: "accepted step",
  x: s_acc, y: h_acc,
  marker: { color: "#2ca02c", size: 6, symbol: "circle" },
  yaxis: "y"
};

const hRejTrace = {
  type: "scatter", mode: "markers", name: "rejected step",
  x: s_rej, y: h_rej,
  marker: { color: "#d62728", size: 8, symbol: "x", line: { width: 2 } },
  yaxis: "y"
};

const barTrace = {
  type: "scatter", mode: "lines", name: "shrink chain",
  x: bar_x, y: bar_y,
  line: { color: "rgba(140,140,140,0.6)", width: 1 },
  yaxis: "y", hoverinfo: "skip", showlegend: false
};

const ratioColors = acc_flag.map(a => a ? "#2ca02c" : "#d62728");
const ratioTrace = {
  type: "scatter", mode: "markers", name: "err / tol",
  x: s_all, y: ratio,
  marker: { color: ratioColors, size: 5, symbol: "circle" },
  yaxis: "y2"
};

const bar2Trace = {
  type: "scatter", mode: "lines", name: "shrink chain",
  x: bar2_x, y: bar2_y,
  line: { color: "rgba(140,140,140,0.6)", width: 1 },
  yaxis: "y2", hoverinfo: "skip", showlegend: false
};

const threshTrace = {
  type: "scatter", mode: "lines", name: "threshold (err = tol)",
  x: [$(s0), $(s1)], y: [1.0, 1.0],
  line: { color: "#888888", width: 1.5, dash: "dash" },
  yaxis: "y2", showlegend: true
};

const layout = {
  title: { text: $(repr(title)), font: { size: 15 } },
  grid: $(layout_grid),
  xaxis:  $(xaxis_block),
  yaxis:  { title: "step size h", type: "log", domain: $(yaxis_domain) }$(threshold_zero_trace),
  yaxis2: { title: "err / tol", type: "log", domain: $(yaxis2_domain) }$(xaxis3_block),
  legend: { x: 1.02, y: 1.0, xanchor: "left" },
  margin: { l: 70, r: 160, t: 60, b: 50 },
  annotations: [
    { text: "Step size", xref: "paper", yref: "paper", x: 0.01, y: $(ann_step_y),
      xanchor: "left", yanchor: "top", showarrow: false, font: { size: 13 } },
    { text: "Error / tolerance", xref: "paper", yref: "paper", x: 0.01, y: $(ann_err_y),
      xanchor: "left", yanchor: "top", showarrow: false, font: { size: 13 } }$(ann_comp)
  ]
};

Plotly.newPlot("container", [stiffTrace, barTrace, hAccTrace, hRejTrace, bar2Trace, ratioTrace, threshTrace$(plot_traces_extra)],
               layout, { responsive: true, displaylogo: false });
</script>
</body>
</html>
"""
    open(output, "w") do io
        write(io, html)
    end
    return output
end
