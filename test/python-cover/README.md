# test/python-cover

Coverage of the Python → Julia bridge (juliacall / PythonCall) for BIFROST.
These tests exist because BIFROST is a Julia library that we expect to drive
from Python in real workflows, and the juliacall runtime has its own failure
modes (Julia version selection, multithreading, signal handling, project
environments) that the native Julia test suite cannot reach.

The single entry point is `test_julia_call.py`, which runs every probe in
this folder in order:

```bash
cd test/python-cover
uv run python test_julia_call.py
```

Individual probes can also be run directly (`uv run python <probe>.py`).
All scripts assume `julia` is on `PATH` — we pin the juliacall-managed
Julia to the system install via `PYTHON_JULIAPKG_EXE` so it matches the
version used by `julia --project=.`.

## Files

### `julia-call-demo.py` + `julia-call-demo.jl`

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
cd test/python-cover
uv run python julia-call-demo.py
```

### `julia-call-mcm.py` + `julia-call-mcm.jl` + `julia-call-mcm-native.jl`

Targeted probe for `MonteCarloMeasurements.jl` (MCM) under juliacall with
Julia started multi-threaded. BIFROST relies on MCM throughout
[`src/material-properties.jl`](../../src/material-properties.jl),
[`src/fiber/fiber-cross-section.jl`](../../src/fiber/fiber-cross-section.jl),
and the propagator — so any MCM thread-safety or juliacall numerical drift
matters.

The Julia module `MCMDemo` (in `julia-call-mcm.jl`) builds `M` independent
`T ~ Normal(T_nom, T_sigma)` `Particles` ensembles (default
`T_nom = 293 K`, `T_sigma = 10 K`, `N = 2000` samples each), calls
`refractive_index(PURE_SILICA, λ, T)` on each, and reduces the result to
`pmean` / `pstd` scalars. The Python driver passes `T_sigma` explicitly to
both the in-process call and the native subprocess, so both runners
provably use the same value. The same module is driven two ways:

1. From Python (`julia-call-mcm.py`) under juliacall, with
   `PYTHON_JULIACALL_THREADS=auto` so `Threads.nthreads() > 1`. Runs once
   serially and once with `Threads.@threads`; asserts both produce the
   same scalars; writes `output/julia-call-mcm.python.csv`.
2. From a native Julia process (`julia-call-mcm-native.jl`) shelled out by
   the Python driver, using the same module, same RNG seed, same iteration
   order. Writes `output/julia-call-mcm.julia.csv`.

The Python driver then `filecmp`s the two CSVs. They are written with
`@sprintf("%.17g", ...)` so every `Float64` round-trips exactly, and the
RNG is seeded deterministically — so a true match is byte-identical. On
mismatch the driver falls back to a numerical diff and prints the worst
delta before failing.

Outputs:

- `output/julia-call-mcm.python.csv` — written from juliacall.
- `output/julia-call-mcm.julia.csv` — written from native Julia.

Failure modes the probe surfaces:

- juliacall starts Julia single-threaded → `AssertionError`.
- MCM not resolvable in the repo project → `Pkg.instantiate()` raises.
- Serial vs threaded mismatch inside one Julia runtime → real MCM thread
  safety issue; the rows are printed.
- juliacall-driven CSV diverges from native-Julia CSV → juliacall is
  altering numerical behavior (signal handling, RNG state, library
  versions); worst |Δ| is logged.

## Adding a new python-cover test

Follow the file-naming pattern (`<purpose>.py` plus an optional sibling
`<purpose>.jl`). Keep environment setup minimal and portable: rely on
`shutil.which("julia")`, never hard-code paths to a specific machine's
Julia install. If a test writes artifacts, write them under the repo's
`output/` folder (see [`../../ARCHITECTURE.md`](../../ARCHITECTURE.md)
entry `[20]`).
