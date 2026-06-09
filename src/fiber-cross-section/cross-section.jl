"""
Shared base for the concrete cross sections in the `fiber-cross-section`
directory. Provides the `FiberCrossSection` abstract type, the
`BirefringenceResponse` structure, and the common validation methods used by
those files.

Files in this directory intentionally model only quantities that are meaningful
for a single transverse slice of fiber of infinitesimal length. They exclude
any property that depends on fiber length, path through space, accumulated
phase, or concatenation of segments, and depend only on `material-properties.jl`.
"""

#################################################
#
# Abstract structures
#
#################################################

abstract type FiberCrossSection end

struct BirefringenceResponse{T}
    Δβ::T
    dω::T
end

#################################################
#
# Validation utilities
#
#################################################

function validate_positive_length(value::Real, name::AbstractString)
    x = float(value)
    if !(isfinite(x) && x > zero(x))
        throw(ArgumentError("$(name) must be a finite positive value in meters"))
    end
    return x
end

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

function validate_axis_ratio(axis_ratio)
    if !(isfinite(axis_ratio) && axis_ratio >= one(axis_ratio))
        throw(ArgumentError("axis_ratio must be a finite value >= 1 (major/minor)"))
    end
    return axis_ratio
end