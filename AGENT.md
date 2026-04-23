# AGENT.md — Guidance for Automated Agents

## Orientation

Start here, then read ARCHITECTURE.md and README.md.

- **Active code**: `julia-port/` — Julia refactor, all new work goes here.

## Running Tests

From the shell:

```bash
julia julia-port/test/runtests.jl
```

The test orchestrator is `julia-port/test/runtests.jl`. All test files live under
`julia-port/test/`. The test for the path-integral solver (`test_path_integral.jl`) is a
known work in progress.

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

## Key Invariants

Do not break these without explicit user discussion:

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

The model is valid for:

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
- Line wrap markdown and comments at 100 characters.
- Send the output of all demo.jl methods to the folder in the repository root 
named output. 
- If a file has an example usage section at the top, always keep it 
  up to date.  