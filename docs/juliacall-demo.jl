module Zoo

export Animal, describe, daily_food_fraction, scale_intake, weekly_intake, herd_total_mass

struct Animal
    name::String
    species::String
    gender::String      # "F", "M", or "U"
    weight_kg::Float64
    food_kg_per_day::Float64
end

# Multiple dispatch on description
describe(x::Number) = "number: $x"
describe(a::Animal) =
    "$(a.name) the $(a.species) ($(a.gender)): " *
    "$(a.weight_kg) kg, eats $(a.food_kg_per_day) kg/day"

# Daily food as a fraction of body weight
daily_food_fraction(a::Animal) = a.food_kg_per_day / a.weight_kg

# Scale a vector of daily intakes (kg/day) by a factor
scale_intake(intake::AbstractMatrix{<:Real}, α::Real) = α .* intake

# Weekly intake per animal: 7 * food_kg_per_day, returned as a 1D array
weekly_intake(animals::Vector{Animal}) =
    [7.0 * a.food_kg_per_day for a in animals]

# Total herd mass (kg)
herd_total_mass(animals::Vector{Animal}) = sum(a.weight_kg for a in animals)

end
