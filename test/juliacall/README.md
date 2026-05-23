# test/juliacall

Coverage of the Python → Julia bridge (juliacall / PythonCall) for BIFROST.
These tests exist because BIFROST is a Julia library that we expect to drive
from Python in real workflows, and the juliacall runtime has its own failure
modes (Julia version selection, multithreading, signal handling, project
environments) that the native Julia test suite cannot reach.

These probes run as part of the standard Julia test suite via
[`../test_juliacall.jl`](../test_juliacall.jl), which is included from
[`../runtests.jl`](../runtests.jl):

```bash
julia --project=. test/runtests.jl
```

Individual probes can also be run directly (`uv run python <probe>.py`) for
faster iteration. All scripts assume `julia` is on `PATH` — we pin the
juliacall-managed Julia to the system install via `PYTHON_JULIAPKG_EXE` so it
matches the version used by `julia --project=.`.

## Files

### `../../docs/juliacall-demo.py` + `../../docs/juliacall-demo.jl`

Smoke test for the juliacall bridge. Defines a small `Zoo` Julia module
(animals with weight and per-day food intake) and drives it from Python:

- multiple dispatch (`describe(::Number)` vs `describe(::Animal)`),
- Julia structs constructed from Python,
- numpy ↔ Julia matrix round-trip via `scale_intake`,
- vector-of-struct operations (`weekly_intake`, `herd_total_mass`).

If this script runs cleanly, the basics of the bridge are healthy: Julia
selection, package resolution in `.venv/julia_env/`, signal handling, and
data marshalling all work.

Run with
```bash
uv run python docs/juliacall-demo.py
```

### `juliacall-mcm.py` + `juliacall-mcm.jl` + `juliacall-mcm-native.jl`

Targeted probe for `MonteCarloMeasurements.jl` (MCM) under juliacall with
Julia started multi-threaded. BIFROST relies on MCM throughout
[`src/material-properties.jl`](../../src/material-properties.jl),
[`src/fiber/fiber-cross-section.jl`](../../src/fiber/fiber-cross-section.jl),
and the propagator — so any MCM thread-safety or juliacall numerical drift
matters.

The Julia module `MCMDemo` (in `juliacall-mcm.jl`) builds `M` independent
`T ~ Normal(T_nom, T_sigma)` `Particles` ensembles (default
`T_nom = 293 K`, `T_sigma = 10 K`, `N = 2000` samples each), calls
`refractive_index(PURE_SILICA, λ, T)` on each, and reduces the result to
`pmean` / `pstd` scalars. The Python driver passes `T_sigma` explicitly to
both the in-process call and the native subprocess, so both runners
provably use the same value. The same module is driven two ways:

1. From Python (`juliacall-mcm.py`) under juliacall, with
   `PYTHON_JULIACALL_THREADS=auto` so `Threads.nthreads() > 1`. Runs once
   serially and once with `Threads.@threads`; asserts both produce the
   same scalars; writes `output/juliacall-mcm.python.csv`.
2. From a native Julia process (`juliacall-mcm-native.jl`) shelled out by
   the Python driver, using the same module, same RNG seed, same iteration
   order. Writes `output/juliacall-mcm.julia.csv`.

The Python driver then `filecmp`s the two CSVs. They are written with
`@sprintf("%.17g", ...)` so every `Float64` round-trips exactly, and the
RNG is seeded deterministically — so a true match is byte-identical. On
mismatch the driver falls back to a numerical diff and prints the worst
delta before failing.

Outputs:

- `output/juliacall-mcm.python.csv` — written from juliacall.
- `output/juliacall-mcm.julia.csv` — written from native Julia.

Failure modes the probe surfaces:

- juliacall starts Julia single-threaded → `AssertionError`.
- MCM not resolvable in the repo project → `Pkg.instantiate()` raises.
- Serial vs threaded mismatch inside one Julia runtime → real MCM thread
  safety issue; the rows are printed.
- juliacall-driven CSV diverges from native-Julia CSV → juliacall is
  altering numerical behavior (signal handling, RNG state, library
  versions); worst |Δ| is logged.

## Adding a new juliacall test

Follow the file-naming pattern (`<purpose>.py` plus an optional sibling
`<purpose>.jl`). Boot juliacall with:

```python
from bifrost import start
jl = start()                   # or start(threads="auto") for multi-threaded
```

`start()` (defined in [`../../src/bifrost.py`](../../src/bifrost.py)) handles
the env-var ordering juliacall requires, locates the system `julia`, and
activates + instantiates the repo's Julia project. Never hard-code paths to
a specific machine's Julia install. If a probe also needs to shell out to
native Julia, import `julia_exe()` from `bifrost` and use its returned path.
If a test writes artifacts, write them under the repo's `output/` folder (see
[`../../ARCHITECTURE.md`](../../ARCHITECTURE.md) entry `[20]`).
