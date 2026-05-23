if !isdefined(@__MODULE__, :FiberCrossSection)
    include(joinpath(@__DIR__, "fiber-cross-section.jl"))
end
if !isdefined(@__MODULE__, :PathSpecCached)
    include(joinpath(@__DIR__, "..", "geometry", "path-geometry.jl"))
end

if !isdefined(@__MODULE__, :Nickname)
    include(joinpath(@__DIR__, "fiber-path-meta.jl"))
end

"""
Fiber assembly on top of `path-geometry.jl`.

High-level authoring happens in `path-geometry.jl`:
- build piecewise geometry with `PathSpecBuilder`
- compile it to `PathSpecCached` with `build(path_spec)`
- bind geometry to a cross section with `Fiber(path; cross_section, T_ref_K)`

`Fiber` is the compiled query object consumed downstream by `path-integral.jl`.
It owns:
- the immutable built `PathSpecCached`
- the `FiberCrossSection`
- a single reference temperature `T_ref_K` that is the reference for both the
  path geometry (segment lengths/radii valid at `T_ref_K`) and the cross-section
  dimensions (core/cladding diameters valid at `T_ref_K`).
- the fiber domain `[s_start, s_end]`

Operating wavelength `λ_m` is NOT stored on `Fiber`; it is an argument to
`generator_K` / `generator_Kω` (and to `propagate_fiber` in `path-integral.jl`),
so the same `Fiber` can be queried at multiple wavelengths. Temperature is
fixed at `T_ref_K` for all queries; per-segment temperature overrides will be
introduced later via segment metadata (MetaList / MCMadd), not via a T(s) closure.

# ----------------------------
# Example Use
# ----------------------------

xs = FiberCrossSection(
    GermaniaSilicaGlass(0.036),
    GermaniaSilicaGlass(0.0),
    8.2e-6,
    125e-6;
    manufacturer = "Corning",
    model_number = "SMF-like"
)

path_spec = PathSpecBuilder()
straight!(path_spec; length = 5.0)
bend!(path_spec;
    radius = 4.458, angle = π / 2, axis_angle = 0.0,
    meta = [
        Nickname("90° bend"),
        MCMadd(:T_K, Normal(0.0, 2.0)),   # +ΔT_K ~ N(0, 2 K) on this segment
    ],
)
straight!(path_spec; length = 8.0)
# twist!(path_spec; s_start = 0.0, length = 13.0, rate = 0.15)  # TODO: twist refactor

path = build(path_spec)
fiber = Fiber(path; cross_section = xs, T_ref_K = 297.15)

# Operating wavelength is supplied per query; temperature is f.T_ref_K.
K  = generator_K(fiber, 1550e-9)
Kω = generator_Kω(fiber, 1550e-9)
"""

if !isdefined(@__MODULE__, :DEFAULT_T_REF_K)
    const DEFAULT_T_REF_K = 297.15
end

function bend_components(path::PathSpecCached, s::Real)
    κ = curvature(path, s)
    if κ == zero(κ)
        z = zero(κ)
        return (kx = z, ky = z, k2 = z)
    end

    # Path-backed fibers use the path's local normal/binormal frame. The bend
    # axis is the curvature normal, so the local transverse bend components are
    # (κ, 0) in that frame. Frame rotation enters through the path twist rate.
    z = zero(κ)
    return (kx = κ, ky = z, k2 = κ * κ)
end

struct Fiber{P,T,S}
    path::P
    cross_section::FiberCrossSection
    T_ref_K::T
    s_start::S
    s_end::S
end

function Fiber(
    path::PathSpecCached;
    cross_section::FiberCrossSection,
    T_ref_K = DEFAULT_T_REF_K,
)
    s_start, s_end = promote(path.spec.s_start, path.s_end)
    @assert s_end > s_start "Fiber requires s_end > s_start"
    return Fiber{typeof(path),typeof(T_ref_K),typeof(s_start)}(
        path,
        cross_section,
        T_ref_K,
        s_start,
        s_end,
    )
end

# TODO: twist refactor — material_twist is currently a stub; restore once the
# per-segment-meta twist subsystem lands.
path_twist_rate(path::PathSpecCached, s::Real) = geometric_torsion(path, s) + material_twist(path, s)

fiber_path(f::Fiber) = f.path

# ----------------------------
# Generator K(s) and Curvature Kω(s)
# ----------------------------

zero_generator() = zeros(ComplexF64, 2, 2)

function bend_generator_K(f::Fiber, s::Real, λ_m::Real)
    curv = bend_components(f.path, s)
    if curv.k2 == zero(curv.k2)
        return zero_generator()
    end

    T = f.T_ref_K
    R = inv(sqrt(curv.k2))
    Δβb = bending_birefringence(f.cross_section, λ_m, T; bend_radius_m = R)
    c2φ = (curv.kx * curv.kx - curv.ky * curv.ky) / curv.k2
    s2φ = (2 * curv.kx * curv.ky) / curv.k2
    return [
         0.5im * Δβb * c2φ             0.5im * Δβb * s2φ
         0.5im * Δβb * s2φ            -0.5im * Δβb * c2φ
    ]
end

function bend_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    curv = bend_components(f.path, s)
    if curv.k2 == zero(curv.k2)
        return zero_generator()
    end

    T = f.T_ref_K
    R = inv(sqrt(curv.k2))
    Δβbω = bending_birefringence(
        WithDerivative(),
        f.cross_section,
        λ_m,
        T;
        bend_radius_m = R
    ).dω
    c2φ = (curv.kx * curv.kx - curv.ky * curv.ky) / curv.k2
    s2φ = (2 * curv.kx * curv.ky) / curv.k2
    return [
         0.5im * Δβbω * c2φ             0.5im * Δβbω * s2φ
         0.5im * Δβbω * s2φ            -0.5im * Δβbω * c2φ
    ]
end

function twist_generator_K(f::Fiber, s::Real, λ_m::Real)
    tau = path_twist_rate(f.path, s)
    T = f.T_ref_K
    Δβt = twisting_birefringence(f.cross_section, λ_m, T; twist_rate_rad_per_m = tau)
    return [
         0.0           -0.5 * Δβt
         0.5 * Δβt      0.0
    ]
end

function twist_generator_Kω(f::Fiber, s::Real, λ_m::Real)
    tau = path_twist_rate(f.path, s)
    T = f.T_ref_K
    Δβtω = twisting_birefringence(
        WithDerivative(),
        f.cross_section,
        λ_m,
        T;
        twist_rate_rad_per_m = tau
    ).dω
    return [
         0.0           -0.5 * Δβtω
         0.5 * Δβtω     0.0
    ]
end

fiber_breakpoints(f::Fiber) = breakpoints(f.path)

"""
    generator_K(fiber, λ_m) -> (s -> 2×2 ComplexF64)

Return a closure that evaluates the local Jones generator `K(s)` at the given
operating wavelength `λ_m` (metres). Temperature is `fiber.T_ref_K`.
"""
function generator_K(f::Fiber, λ_m::Real)
    return function (s::Real)
        return bend_generator_K(f, s, λ_m) +
               twist_generator_K(f, s, λ_m)
    end
end

"""
    generator_Kω(fiber, λ_m) -> (s -> 2×2 ComplexF64)

Frequency-derivative counterpart of `generator_K`.
"""
function generator_Kω(f::Fiber, λ_m::Real)
    return function (s::Real)
        return bend_generator_Kω(f, s, λ_m) +
               twist_generator_Kω(f, s, λ_m)
    end
end

# ----------------------------
# Fiber diagnostics for plotting
# ----------------------------

function bend_geometry(f::Fiber, s::Real)
    curv = bend_components(f.path, s)
    kx = curv.kx
    ky = curv.ky
    k2 = kx * kx + ky * ky
    if k2 == 0.0
        return (Rb = Inf, theta_b = 0.0, kx = 0.0, ky = 0.0, k2 = 0.0)
    end

    return (Rb = inv(sqrt(k2)), theta_b = atan(ky, kx), kx = kx, ky = ky, k2 = k2)
end

function twist_rate(f::Fiber, s::Real)
    return path_twist_rate(f.path, s)
end
