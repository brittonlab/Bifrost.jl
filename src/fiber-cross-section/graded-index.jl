"""

Local optical properties of a graded-index fiber cross section.

This fiber cross-section type has not yet been cleanly implemented. The file serves only
as a placeholder to show that other cross-section types can be included.

{This file is a stub. You can help Wikipedia by expanding it.}

"""

struct GradedIndexCrossSection{T<:Real} <: FiberCrossSection
    manufacturer::Union{Nothing, String}
    model_number::Union{Nothing, String}
    core_material::AbstractMaterial
    cladding_material::AbstractMaterial
    core_diameter_m::T
    cladding_diameter_m::T

    function GradedIndexCrossSection(
        core_material::AbstractMaterial,
        cladding_material::AbstractMaterial,
        core_diameter_m::Real,
        cladding_diameter_m::Real;
        manufacturer::Union{Nothing, AbstractString} = nothing,
        model_number::Union{Nothing, AbstractString} = nothing
    )
        core_diameter = validate_positive_length(core_diameter_m, "core diameter")
        cladding_diameter = validate_positive_length(cladding_diameter_m, "cladding diameter")
        core_diameter < cladding_diameter || throw(ArgumentError(
            "core diameter must be smaller than cladding diameter"
        ))

        T = promote_type(typeof(core_diameter), typeof(cladding_diameter))
        new{T}(
            isnothing(manufacturer) ? nothing : String(manufacturer),
            isnothing(model_number) ? nothing : String(model_number),
            core_material,
            cladding_material,
            core_diameter,
            cladding_diameter
        )
    end
end