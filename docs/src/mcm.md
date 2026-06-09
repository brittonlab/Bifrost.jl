# Monte Carlo Measurements Compatibility

Several files are written to lift through
`MonteCarloMeasurements.Particles`:

- `src/material/material-properties.jl` and `materials/` files
- `src/fiber/fiber-cross-section.jl` and `fiber-cross-sections/` files
- `src/geometry/path-geometry.jl`
- `src/fiber/fiber-path.jl`
- `src/fiber/fiber-path-modify.jl`
- `src/path-integral.jl`

Important conventions:

- avoid `::Real` annotations on uncertain-input slots,
- avoid `Float64(...)` coercions on paths that may carry `Particles`,
- avoid conditionals that would need to branch independently per particle,
- use `MonteCarloMeasurements.unsafe_comparisons(true)` in MCM tests,
- reduce ensemble-wide adaptive decisions through `scalar_reduce`.

`scalar_reduce` currently uses the maximum particle value when reducing an
MCM-valued scalar error metric. That makes the adaptive controller conservative:
the whole ensemble takes one step size selected by the worst particle. A
`pmean`-style reduction is the documented performance compromise if worst-case
step counts become too high.

Per-segment uncertainty and annotations live in the `meta` vector:

- `Nickname(label)` labels a segment for visual diagnostics,
- `MCMadd(symbol, distribution)` applies additive perturbations,
- `MCMmul(symbol, distribution)` applies multiplicative perturbations.

`src/fiber/fiber-path-modify.jl` interprets those annotations. For example,
`:T_K` metadata is converted into thermal length scaling using the cladding
material CTE at the fiber reference temperature.
