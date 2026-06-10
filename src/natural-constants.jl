"""
natural-constants.jl

Single source of truth for fundamental physical constants.

Values are taken from `PhysicalConstants.CODATA2022` and stripped to plain
`Float64` (SI units) at definition time. Stripping is deliberate: `CODATA2022`
returns Unitful `Quantity` values, which are not `<:Real` and would break the
`MonteCarloMeasurements.Particles` code paths in the material, cross-section,
and propagation layers. Plain `Float64` fields keep those paths `::Real`-clean.

Constants are bundled in the immutable [`NaturalConstants`](@ref) struct and
exposed as the `const` singleton [`u`](@ref), accessed by field: `u.c`, `u.ħ`,
`u.kB`. Because `u` is a `const` of an immutable struct with concrete fields,
field access constant-folds and is type-stable with no runtime overhead.
"""

import PhysicalConstants.CODATA2022 as CODATA
using Unitful: ustrip, @u_str

"""
    NaturalConstants

Bundle of fundamental physical constants in SI units, stored as `Float64`.

# Fields
- `c`: speed of light in vacuum, m/s
- `ħ`: reduced Planck constant, J·s
- `kB`: Boltzmann constant, J/K
"""
struct NaturalConstants
    c::Float64
    ħ::Float64
    kB::Float64
end

"""
    u

Singleton bundle of fundamental physical constants (CODATA 2022, SI units).
Access fields directly:

```jldoctest
julia> using Bifrost

julia> u.c
2.99792458e8
```

See [`NaturalConstants`](@ref) for the available fields (`c`, `ħ`, `kB`).
"""
const u = let
    c  = ustrip(u"m/s", float(CODATA.SpeedOfLightInVacuum))
    ħ  = ustrip(u"J*s", float(CODATA.ReducedPlanckConstant))
    kB = ustrip(u"J/K", float(CODATA.BoltzmannConstant))
    NaturalConstants(c, ħ, kB)
end
