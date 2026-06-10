#################################################
#
# Abstract structures
#
#################################################

"""
    FiberCrossSection

Supertype for transverse fiber cross sections.

Concrete cross sections model only quantities that are meaningful for a single
transverse slice of fiber of infinitesimal length — guided-mode quantities and
local birefringence response magnitudes — and depend only on the material layer.
Anything that depends on fiber length, path through space, accumulated phase, or
concatenation of segments belongs to the fiber assembly layer.
"""
abstract type FiberCrossSection end

"""
    BirefringenceResponse(Δβ, dω)

Pair a local birefringence magnitude `Δβ` (rad/m) with its angular-frequency
derivative `dω`.
"""
struct BirefringenceResponse{T}
    Δβ::T
    dω::T
end

#################################################
#
# Validation utilities
#
#################################################

"""
    validate_positive_length(value, name) -> Float64

Return `value` as a `Float64`, throwing an `ArgumentError` naming `name` unless
it is a finite positive length (m).
"""
function validate_positive_length(value::Real, name::AbstractString)
    x = float(value)
    if !(isfinite(x) && x > zero(x))
        throw(ArgumentError("$(name) must be a finite positive value in meters"))
    end
    return x
end

"""
    validate_bend_radius(bend_radius_m)

Return `bend_radius_m`, throwing an `ArgumentError` unless it is a finite
positive radius (m) or `Inf` (a straight fiber).
"""
function validate_bend_radius(bend_radius_m)
    if isinf(bend_radius_m) && bend_radius_m > zero(bend_radius_m)
        return bend_radius_m
    end
    if !(isfinite(bend_radius_m) && bend_radius_m > zero(bend_radius_m))
        throw(ArgumentError(
            "bend_radius_m must be a finite positive value in meters or Inf"
        ))
    end
    return bend_radius_m
end

"""
    validate_axis_ratio(axis_ratio)

Return `axis_ratio`, throwing an `ArgumentError` unless it is a finite
major/minor ellipse axis ratio `≥ 1`.
"""
function validate_axis_ratio(axis_ratio)
    if !(isfinite(axis_ratio) && axis_ratio >= one(axis_ratio))
        throw(ArgumentError("axis_ratio must be a finite value >= 1 (major/minor)"))
    end
    return axis_ratio
end