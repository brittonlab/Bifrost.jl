"""
    GradedIndexCrossSection(core_material, cladding_material, core_diameter_m,
                            cladding_diameter_m; manufacturer = nothing,
                            model_number = nothing)

Placeholder graded-index fiber cross section.

Stores geometry and materials only; no guided-mode or birefringence model is
implemented yet, and the fiber-level generators error for this type. It exists
to show that cross-section types beyond step-index can be added.
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