# Developing

This page is the source of authority for contributor workflow, extension recipes, and
the Monte Carlo Measurements (MCM) compatibility contract. Repository-root files such as
`AGENTS.md` and `ARCHITECTURE.md` summarize the hard invariants for agent tooling and
cross-reference this page.

## Development setup

Install Julia and instantiate the project environment as described in the repository
`README.md` (Installation section).

The repository owns reusable agent skills under `skills/`. The Julia documentation skill
lives at `skills/julia-docstrings/` and should be used when creating, revising, or
auditing inline documentation for Julia code. Register repo-owned skills with your local
agent tools by running:

```bash
skills/install-agent-skills.sh
```

The installer symlinks the repo skill into local Codex and Claude skill directories, so
updates from `git pull` are picked up without copying files. Use
`skills/install-agent-skills.sh --dry-run` to preview the target paths.

## Running the tests

Run the Julia test suite from the repository root:

```bash
julia --project=. test/runtests.jl
```

Tests are classified by the taxonomy in `AGENTS.md` (`T-PHYSICS`, `T-VALIDATION`,
`T-SIM-REGRESSION`, `T-GUARDRAIL`, and visual demos); place new tests in one of those
categories and follow its standard. Any change that touches an uncertain-input code path
must add or extend an MCM test — see [Testing MCM changes](@ref) below.

## Building the documentation

Documentation is built with [Documenter.jl](https://documenter.juliadocs.org). Source
markdown lives under `docs/src/`; `docs/make.jl` defines the build and navigation. A
GitHub Actions workflow (`.github/workflows/Documenter.yml`) rebuilds and deploys the
site to GitHub Pages.

Build the docs locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The rendered site is written to `docs/build/` (open `docs/build/index.html`).

## Extending Bifrost

The package is layered with strictly hierarchical dependencies - material ≤ cross section ≤ fiber, geometry ≤ fiber — with the
propagation layer consuming callable generators. `ARCHITECTURE.md` (Layered Design)
maps the layers to files. Each extension point below names the layer it belongs to.
Whatever you add, write code that satisfies the
[MCM compatibility](@ref mcm-compatibility) contract from the start: retrofitting
`Particles` support is much harder than designing for it.

### Adding a birefringence source

Fiber assembly layer (`src/fiber/fiber-path.jl`).

1. Define a struct as a subtype of `AbstractBirefringenceSource`.
2. Implement `generator_K_contribution(source, s)` — returns the local 2×2 generator
   contribution.
3. Implement `generator_Kω_contribution(source, s)` — returns ∂K/∂ω. May return a zero
   matrix if not yet modeled.
4. Declare `coverage_intervals(source)` and `breakpoints(source)`. Every source must
   cover the full fiber domain `[s_start, s_end]`; gaps are a hard error, not silent
   zero. The `Fiber` merges breakpoints globally and the propagator never steps across
   one.
5. Extend the fiber-level generator assembly so the new source contributes to `K(s)`
   and `Kω(s)`.
6. Add guardrail tests *before* physics tests.

Changing the source interface contracts themselves (`generator_K_contribution`,
`generator_Kω_contribution`) requires user authorization — see `AGENTS.md`.

### Adding a material

Material layer (`src/material/`).

Implement `refractive_index(::ValueOnly, material, λ, T_K)` and
`refractive_index(::WithDerivative, material, λ, T_K)`. For Sellmeier-form materials the
helpers in `src/material/material-properties.jl` (`sellmeier_index_from_coefficients`,
`sellmeier_index_from_coefficients_dω`, and the `_evaluate_sellmeier_*` utilities for
parameter-dependent coefficients) do the heavy lifting; see `silica.jl` and
`germania.jl` for worked examples. Routing all spectral evaluation through the Sellmeier
helpers also gives `Particles` compatibility for free — the helpers are branch-free and
coercion-free.

### Adding a cross section

Cross-section layer (`src/fiber-cross-section/`).

The base interface lives in `src/fiber-cross-section/cross-section.jl`; use
`step-index.jl` as the model implementation. A cross section converts material
properties into guided-index, dispersion, nonlinearity, and local birefringence response
coefficients, and must keep operating wavelength a per-query argument rather than
stored state.

## [MCM compatibility](@id mcm-compatibility)

BIFROST propagates measurement uncertainty with
[MonteCarloMeasurements.jl](https://github.com/baggepinnen/MonteCarloMeasurements.jl)
(MCM): an uncertain input is a `Particles` ensemble that flows through the entire
computation, so one simulation yields a distribution of outputs. This only works if
every function on the path of an uncertain value lifts through `Particles`. This section
is the contract for keeping it that way.

### What must lift

`Particles` may arrive on these inputs:

- temperature: the user sets it through the thermal `:T_K` segment annotation (and
  the baseline `T_ref_K`); `temperature(f, s)` resolves these to the internal `T_K`
  slot that the cross-section and material layers consume, so that slot must lift too,

- bend, twist, tension, and axis-ratio properties,
- segment shrinkage and field-level `MCMadd`/`MCMmul` perturbations,
- the per-entry element type of the Jones matrix `J` and sensitivity `G`.

The files on those paths — all of `src/material/`, `src/fiber-cross-section/`,
`src/geometry/path-geometry*.jl`, `src/fiber/fiber-path.jl`, and
`src/path-integral.jl` — must obey the rules below on every uncertain-input slot.

### The rules

1. **No annotations or coercions that exclude `Particles`.** Leave uncertain-input
   slots unannotated where possible. This is because `Particles <: Real` (but not `<: AbstractFloat`). Use of  `::Float64` or `::AbstractFloat` silently
   excludes the ensemble and degrade performance. Never coerce with `Float64(·)` on a path that may carry  `Particles`: promotion does the right thing, coercion destroys the ensemble.
2. **No per-particle branching.** A conditional on a `Particles` value is an error
   (which branch would the ensemble take?). Either write branch-free code (`sign`,
   `flipsign`, `clamp`, elementwise arithmetic), or reduce to a representative scalar
   first and branch on that: `scalar_reduce` in `src/path-integral.jl` (worst-case
   `pmaximum` reduction) and `_qc_nominalize` in
   `src/geometry/path-geometry-connector.jl` (`pmean` nominalization) are the two
   established patterns. Loop bounds and quadrature limits must be deterministic
   `Float64`.
3. **No generic dense linear algebra on `Particles` matrices.** `LinearAlgebra.exp`
   (Padé with pivoting), `opnorm`, and `eigvals` hit code paths that do not lift. 
   Bifrost.jl implements private variants that work with MCM.
   For example `exp_jones_generator` and
   `exp_block_upper_triangular_2x2` (a closed-form Fréchet derivative of the 2×2
   exponential) for propagation, and `output_dgd_2x2` for DGD extraction
   (`output_dgd` uses `eigvals` and is Float64-only).
4. **Nominalize deliberately, not defensively.** Reducing a `Particles` value to a
   scalar is correct for *control flow* (step acceptance, table bisection, Newton
   iteration counts) but wrong for *results* — an accumulated phase or coefficient
   must keep the full ensemble. When you nominalize, say why in a short comment.

### Ensemble-wide adaptive decisions

The adaptive step controller makes one decision per step for the whole ensemble — it
cannot take half a step for some particles and a full step for others. `scalar_reduce`
collapses the error metric via `pmaximum`, so the worst-case particle drives the step
size. That is conservative; a `pmean`-based reduction is the documented performance
compromise if step counts become too large under tight tolerances. `scalar_reduce` is
the single switch point for that change.

### Testing MCM changes

`test/test_mcm_compatability.jl` is the primary MCM suite; thermal, tension, and
perturbation coverage also lives in `test/test_fiber_thermal.jl`,
`test/test_fiber_tension.jl`, and `test/test_path_geometry_perturb.jl`. When you touch
an uncertain-input code path, extend one of these (or add a new MCM test) so the change
is exercised with `Particles` inputs, not just `Float64`.

Conventions:

- wrap MCM test blocks in `MonteCarloMeasurements.unsafe_comparisons(true)`;
- under unsafe comparisons, structural invariants ("breakpoints are sorted and
  deduplicated") reduce via `pmean` rather than failing;
- compare `pmean` of a `Particles` result against the nominal `Float64` baseline.

