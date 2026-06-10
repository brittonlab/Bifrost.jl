# Project Structure

This is a high-level schematic. Do not update it to reflect every file.

```text
.
├── .cursor                          [1]
├── AGENTS.md                        [2]
├── ARCHITECTURE.md
├── README.md                        [3]
├── Project.toml
├── Manifest.toml
├── src                              [8]
│   ├── path-integral.jl             [11]
│   ├── material                     [9]
|   ├── material/nonlinear           [9b]
│   ├── geometry                     [10]
│   ├── fiber                        [12, 13, 14]
│   └── nonlinear                    [16]
├── test                             [19]
│   ├── human                        [15]
│   └── legacy-python                [4L]
├── docs                             [7]
└── output                           [20]
```

- [2] Notes for tooling workflows.
- [3] Primary project overview and scientific context.
- [4L] Legacy Python implementation for birefringence simulation.
- [7] Documenter.jl documentation, research references, and source material.
- [8] Active Julia source tree and solver architecture.
- [9] Standalone material models and refractive-index behavior.
- [9b] Nonlinear material properties. 
- [10] Standalone path construction and differential geometry.
- [11] Generic adaptive propagation for callable Jones generators.
- [12] Cross-sectional fiber optics and local birefringence responses.
- [13] Path-backed fiber assembly and generator construction.
- [14] 3D geometry and propagation visualization pipeline.
- [15] Runnable Julia demos; visual demos write HTML files to `output/`.
- [16] Reserved Julia nonlinear namespace; legacy Raman and Brillouin Python
  scripts remain under `test/legacy-python/`.
- [19] Julia tests.
- [20] Output of demo methods and generated visual artifacts.

The folder marked `[L]` contains legacy files for the old Python implementation. Do not
read them unless a specific workflow requires it. They are authoritative for
legacy behavior and must not be modified without explicit user authorization.

## Architectural Intent

- Dependencies between modules is strictly limited. Here, A <= B means B depends on A. 
  - material <= fiber-cross-section <= fiber 
  - geometry <= fiber 
- Keep the core propagation API usable with any callable `K(s)` and `Kω(s)`.
- Support continuous/function-valued geometry and spin rather than only fixed
  pre-sliced segment grids.
- Keep lossless Jones propagation isolated from any future gain/loss model.
- Preserve MCM compatibility on uncertainty-carrying code paths.

## Standalone Building Blocks

These files are intentionally useful on their own:

| File | Standalone role |
| --- | --- |
| `material-properties.jl` | Material constants and spectra; no path or fiber geometry. |
| `path-geometry.jl` | Three-dimensional path construction and geometric queries; no optics. |
| `path-integral.jl` | Adaptive propagation for callable `K(s)` and `Kω(s)` generators. |

The fiber-specific layers combine those pieces:

| File | How it extends the standalone pieces |
| --- | --- |
| `fiber-cross-section/` | Adds step-index fiber optics and birefringence responses (`cross-section.jl` base + concrete cross sections). |
| `fiber-path.jl` | Binds path geometry to a cross section and assembles bend/spin `K` and `Kω`. |

## Layered Design

0. **Geometry layer** (`geometry/path-geometry.jl`,
   `geometry/path-geometry-connector.jl`, `geometry/path-geometry-meta.jl`,
   `geometry/path-geometry-perturb.jl`, `geometry/path-geometry-plot.jl`)

   - Builds and queries three-dimensional paths.
   - Authoring lifecycle on a `SubpathBuilder`: `start!` → segment-adding calls
     (`straight!`, `bend!`, `helix!`, `catenary!`, `jumpby!`) → seal
     (`jumpto!` to a global target, or `seal!` to end at the natural exit) →
     `build()`.
   - `build(Subpath(builder)) → SubpathBuilt` compiles to an immutable form;
     `build(::Vector{Subpath}) → PathBuilt` concatenates multiple independent
     subpaths under a shared global arc length.
   - `start!(builder, :inherit)` (or per-field `point=:inherit`, etc.) defers a
     subpath's start state to the previous subpath's endpoint; `:inherit`
     is rejected on the first
     subpath and on any standalone build (no predecessor). Hand-loaded start
     coordinates remain valid and are still validated by the conformity check.
   - `build(...; perturb=true)` applies the field-level `MCMadd`/`MCMmul` that
     name a segment's own fields (the mechanism lives in
     `geometry/path-geometry-perturb.jl`); `perturb=false` (default) is nominal.
   - **Invariant:** the geometry layer carries any meta it cannot interpret
     blindly and never errors on it — in particular it never references `:T_K`.
     Interpretation of foreign meta is a consuming layer's job (the fiber).
   - Resolves the per-Subpath `spin_rate` into path-coordinate material spin.
   - Resolves `JumpBy` and the terminal `jumpto!` connector into G2 quintic
     connectors at build time.
   - The `AbstractMeta` vocabulary (`Nickname`, `MCMadd`, `MCMmul`)
     lives in `geometry/path-geometry-meta.jl`. It makes no reference to fiber.
     (Material spin is a `start!(; spin_rate=…)` keyword, not meta.)

1. **Material layer** (`material-properties.jl`)

   - Encodes intrinsic optical material properties.
   - Provides spectral responses and derivatives needed by DGD calculations.

2. **Cross-section layer** (`fiber-cross-section/cross-section.jl` and concrete
   cross sections in `fiber-cross-section/`)

   - Encodes transverse step-index fiber geometry.
   - Converts material properties into guided-index, dispersion, nonlinearity,
     and local birefringence response coefficients.

3. **Fiber assembly layer** (`fiber/fiber-path.jl`)

   - `Fiber(spec; cross_section, T_ref_K)` accepts authored geometry (a
     `SubpathBuilder`, `Subpath`, or `Vector{Subpath}`) and builds it once; it
     also binds an already-built `SubpathBuilt`/`PathBuilt` as-is.
   - Sole interpreter of the thermal `:T_K` meta: computes
     `α_lin = cte(cladding_material, T_ref_K)` (lazily, only when a `:T_K`
     segment is present), bakes the isotropic length scaling `1 + α_lin·ΔT` into
     the affected segments, strips `:T_K`, and lets the geometry build apply any
     field-level `MCMadd`/`MCMmul`. (`modify` has been removed.)
   - A `jumpto!` seal may itself carry `:T_K`: the terminal connector then
     thermally expands — its arc length scales by τ while still landing at the
     fixed `jumpto_point` — by passing `build(...; jumpto_target_length=τ·L0)`.
     `min_bend_radius` is still honored (validated post-hoc when a
     target length is set).
   - Keeps operating wavelength as a per-query argument rather than `Fiber`
     state.
   - Assembles fiber-level bend and spin generators `K(s)` and `Kω(s)`.

4. **Propagation layer** (`path-integral.jl`)

   - Solves `dJ/ds = K(s)J` with adaptive step-doubling exponential midpoint
     integration.
   - Solves the coupled sensitivity system for `G = ∂ωJ`.
   - Derives DGD from `J` and `G`.
   - Uses breakpoint-aware interval decomposition to avoid integrating across
     discontinuities.
   - Uses phase-insensitive error metrics for Jones propagation.

5. **Presentation layer** (`fiber-path-plot.jl`, `demo*.jl`)

   - Generates visual diagnostics and runnable examples.
   - Keeps visual demos as human-inspected outputs rather than reusable library
     code.

## Runtime Flow

0. See `test/human/demo-smallest.jl` for the smallest runnable example.
1. Build a `SubpathBuilder` with path primitives and optional metadata
   (`start!` → segment calls → `jumpto!` or `seal!`).
2. Bind it into `Fiber(builder; cross_section, T_ref_K)` — the constructor
   builds the geometry once (interpreting `:T_K` thermal meta and applying
   field-level MCM) and accepts a `SubpathBuilder`, `Subpath`, or
   `Vector{Subpath}`. (To inspect nominal geometry directly, use
   `build(Subpath(builder))`, or `build(...; perturb=true)` for field-MCM.)
3. Propagate with `propagate_fiber(fiber; λ_m=...)` for Jones output.
4. Use `propagate_fiber_sensitivity(fiber; λ_m=...)` when DGD is needed.
5. Post-process outputs for diagnostics, plots, demos, and regression checks.

## Contracts and Invariants

- Path breakpoints are normalized and globally merged before piecewise
  propagation.
- The propagator must not step across path segment or spin-run boundaries.
- Numerical tolerances (`rtol`, `atol`, step controls) are explicit API inputs,
  not hidden globals.
- Global phase-insensitive error metrics are used in adaptive acceptance checks.
- `path-integral.jl` assumes lossless Jones propagation.
- MCM-compatible paths must avoid scalar coercions and particle-dependent
  conditionals on uncertainty-carrying values.

## Testing Strategy

- Test entrypoint: `julia --project=. test/runtests.jl`.
- Emphasis areas:
  - path geometry and connector invariants,
  - material and cross-section calculations,
  - fiber assembly and breakpoint behavior,
  - propagation behavior and DGD computation,
  - MCM compatibility,
  - paddle-like path construction.
- Tests should follow the taxonomy in `AGENTS.md`: `T-PHYSICS`,
  `T-VALIDATION`, `T-SIM-REGRESSION`, `T-GUARDRAIL`, and visual demos.
- Validation against published fiber data and direct legacy `fiber.py`
  comparisons is still limited.

## Extension Guidance

- Add new path shapes by implementing the `AbstractPathSegment` local geometry
  interface in `path-geometry.jl`, and declare its length-dimensioned fields via
  `_length_fields` in `path-geometry-perturb.jl` so isotropic scaling is defined
  (the fallback errors loudly if omitted).
- Add new per-segment annotations by extending the `AbstractMeta` vocabulary and
  keeping interpretation in the consuming layer — the geometry layer must carry
  meta it cannot interpret blindly (never naming it, never erroring on it).
- Add new fiber-level birefringence mechanisms by extending generator assembly
  in `fiber-path.jl` and adding guardrail tests first.
- Keep solver changes in `path-integral.jl` deliberate; step controller, error
  metric, and exponential formulas are core numerical contracts.
- Preserve separation between lossless Jones propagation and any future
  attenuation, gain, or polarization-dependent loss model.

## Operational Best Practices

- Keep architecture docs schematic; avoid drifting into exhaustive file
  inventory.
- When a feature is still moving, delay broad docs/interface updates until the
  design has settled.
- When a feature settles, update docs and tests together.
- Favor deterministic demos and tests so numerical regressions are quickly
  detectable.
