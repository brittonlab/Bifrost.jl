# Guidance for Automated Agents
AGENTS.md is read by Codex by default. 
Read @ARCHITECTURE.md. 

## Running Tests
All test files live under `test/`.

The test orchestrator is `test/runtests.jl`.

    ```bash
    julia --project=. test/runtests.jl
    ```

## Test Taxonomy

When writing or evaluating tests, place them in one of four categories. Each serves a
different purpose and has different standards.

### 1. Physics-Motivated Tests

Derive expected results from first principles or analytic limits. These are the
highest-value tests. Mark them with the comment "T-PHYSICS".

Examples:

- A straight fiber with no twist produces an identity Jones matrix.
- Pure twist at rate τ over length L produces a rotation matrix with angle ∝ τL.
- Bending at small radius produces birefringence that scales as 1/R².
- DGD of a constant-birefringence fiber equals the analytic beat-length formula.

Standard: The expected value should be computed from the physics equation in a comment,
not inferred from the code output.

### 2. Validation Tests

Compare against published external data (datasheets, papers, known fibers). Mark them
with the comment "T-VALIDATION".

Examples:

- Corning SMF-28 chromatic dispersion vs. datasheet values.
- Refractive index vs. Sellmeier coefficients from the literature.

Standard: Cite the source. Use tolerances that reflect the precision of the published
data, not the precision of the solver. Known acceptable mismatches should be documented
as such.

### 3. Simulation Regression Tests

Run a full propagation over a reference fiber and check that outputs are numerically
stable across code changes. Mark them with the comment "T-SIM-REGRESSION".

Examples:

- Propagate a standard bend-twist fiber at 1550 nm and verify the final Jones matrix and
  DGD against stored expected values.

Standard: Tolerance should be set to something physically meaningful (e.g., DGD to 1 fs,
Jones matrix elements to 1e-8). Update the reference values intentionally, not
silently, when the physics model changes.

### 4. Guardrail Tests

Prevent agents from introducing silent regressions in invariants that the code enforces.
Mark them with the comment "T-GUARDRAIL".

Examples:

- A `Fiber` with a source that does not cover the full domain raises an error.
- Breakpoints are sorted and deduplicated.
- `propagate_fiber` and `propagate_piecewise` produce identical results for the same
  fiber.
- `output_dgd` returns a non-negative scalar.

Standard: These should be cheap, structural, and not depend on tuned numerical
thresholds.

### 5. Visual Tests

These are tests that generate a plot or other visual output for debugging by humans.
They live in the notebook `test/human/bifrost-demos.ipynb`; its demo narratives are
recorded in `test/human/demo-intent.md`, and all plotting machinery lives in
`test/human/demo-helper.jl`.
- Demo cells are not reusable library code; the Bifrost-API portion of each cell is
  written to be read, and duplication across cells is deliberate.
  - **Inline `Nickname` meta directly on every segment call** rather than
    deriving names at render time from type information.  The name lives with the
    segment definition, not in a helper function.

    ```julia
    # good
    PG.straight!(spec; length = 0.1, meta = [PG.Nickname("Straight")])
    PG.bend!(spec;   radius = 0.05, angle = π/2, meta = [PG.Nickname("Bend")])

    # bad — abstracts what should be readable at a glance
    function _segment_type_label(seg) ... end
    ```
  - For segments whose role is obvious from context (lead-in, sag, etc.) use a
  descriptive nickname instead of just the type name:
  `PG.Nickname("lead-in")`, `PG.Nickname("90° bend")`.
- Abstraction related to visual presentation is fine and preferred — that is what
  `demo-helper.jl` is for; helpers there may be as abstract as convenient.
- One atomic demo, one cell, one figure. Markdown cells carry each demo's
  description immediately above its implementation; section order mirrors the
  class index at the top of the notebook.

## Key Invariants

Do not break these without explicit user discussion:

- **No Manifest.toml commits**: Never commit `Manifest.toml`; it is gitignored and only
  snapshotted on tagged releases (issue #77). CI rejects PRs that track it.
- **Source coverage**: Every `AbstractBirefringenceSource` must cover the full fiber
  domain `[s_start, s_end]`. Gaps are a hard error, not silent zero.
- **Breakpoints**: Each source declares its own breakpoints. The `Fiber` merges them
  globally. The propagator never steps across a breakpoint.
- **Lossless assumption**: `path-integral.jl` assumes Jones matrices are in SU(2). Do
  not introduce gain/loss there; that requires a separate module.
- **Phase-insensitive error**: Adaptive step-doubling uses `phase_insensitive_error`, not
  raw matrix difference. This is intentional.
- **Function-valued inputs**: Physical profiles (bend radius, twist rate, temperature,
  axis angle) must be callable at arbitrary `s`, not just on a fixed grid.
- **MCM compatibility**: `material-properties.jl`, the `fiber-cross-section/`
  files, `path-geometry.jl`, and `path-integral.jl` accept
  `MonteCarloMeasurements.Particles` on the uncertain inputs (`T_K`,
  bend/twist/tension/axis-ratio properties, segment shrinkage, and the
  per-entry eltype of the Jones matrices `J` and sensitivity `G`). Keep
  these files `::Real`-free on uncertain-input slots and avoid `Float64(·)`
  coercions on those paths. Test files using MCM must wrap blocks in
  `MonteCarloMeasurements.unsafe_comparisons(true)`; under unsafe comparisons,
  invariants like "breakpoints are sorted and deduplicated" reduce via `pmean`
  rather than failing. In `path-integral.jl` specifically:
  - The adaptive step controller reduces the error metric through
    `scalar_reduce` (pmaximum under MCM) so the ensemble takes one step at a
    time. A `pmean`-based reduction is a reasonable performance compromise if
    step counts become too large under tight tolerances.
  - The 4×4 sensitivity exp is a closed-form Fréchet derivative of the 2×2
    exp (`exp_block_upper_triangular_2x2`), not `LinearAlgebra.exp` — this is
    required for MCM compatibility and is also faster/exacter for Float64.
  - `output_dgd_2x2` is the MCM-friendly DGD extractor (closed form, no
    `eigvals`); `output_dgd` still uses `eigvals` and is Float64-only.
  - MCM prohibits conditionals anywhere that needs to propagate Particles. 
    Here this includes most methods in material-properties.jl, the
    fiber-cross-section/ files, path-geometry.jl, fiber-path.jl, path-integral.jl.
- **python wrapper**: A Python wrapper for `Bifrost.jl` is supplied in `src/wrapper.py`. The spirit 
    for this wrapper is that Python is the workspace and Julia modules are guest libraries. No Julia
    methods should be used in th python environment except those explicitly exposed via
    `wrapper.wrap()`. It is never acceptable to expose general julia methods
    in the python workspace.   

## Adding a New Birefringence Source

1. Define a struct as a subtype of `AbstractBirefringenceSource`.
2. Implement `generator_K_contribution(source, s)` — returns the local 2×2 generator
   contribution.
3. Implement `generator_Kω_contribution(source, s)` — returns ∂K/∂ω. May return zero
   matrix if not yet modeled.
4. Declare `coverage_intervals(source)` and `breakpoints(source)`.
5. Add guardrail tests before adding physics tests.

## What Requires User Authorization

- Starting any task listed in TODO.md.
- Changing the solver algorithm in `path-integral.jl` (step controller, error metric,
  exponential formula).
- Changing the source interface contracts (`generator_K_contribution`,
  `generator_Kω_contribution`).
- Updating interfaces that affect demo.jl or the existing test suite.
- Any modification to legacy Python files.

## Physics Scope

The model must be run only be run in the following limits. 

- Single-mode fiber: V < 2.405
- Weakly guiding fiber: Δn ≪ 1
- Nearly circular core: ellipticity e² ≪ 1
- Bend radii R ≫ r_cladding
- Temperature: ~200–300 K
- Wavelength: 1–2 μm

Mechanisms currently **not modeled**: cladding noncircularity, non-concentric core,
external stress, electric/magnetic fields, polarization-dependent loss, nonlinear
scattering.

## Best practices

- When citing literature only use sources that you can verify in a library catalogue or
  database.
- Markdown, comments and code must line wrap at 92 characters.
- Modules, methods, and structs should have inline documentation following Julia norms.
  Write for package users and current behavior; avoid development history, migration
  notes, and references to retired APIs or design choices.
- Use the `julia-docstrings` skill (registered at the account level) for Julia docstring
  creation, revision, and audits.
- Obey the rule JuliaIndexFromLength. Use eachindex() not length() for iterating over
  vectors/tensors in loops.
- Citations of literature should include doi in only this format doi:xxxxx