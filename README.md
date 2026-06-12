# BIFROST

BIFROST (Birefringence In Fiber: Research and Optical Simulation Toolkit) is a
Julia codebase for simulating polarization mode dispersion in optical fibers.
Silica-based fibers whose core and/or cladding are doped with germania can be
simulated.

The active implementation is a Julia refactor of the original Python
polarization model. Legacy Python code is retained under `test/legacy-python/`
as physics reference material and should not be edited during routine Julia
work.  

The major architectural change is that the optical fiber is not represented as
a pre-sliced list of Jones matrices. Instead, the code builds a continuous
centerline path, binds it to a transverse fiber cross section, and integrates a
local Jones generator:

```math
\frac{dJ}{ds}=K(s,\omega)J,\qquad J(s_0)=I.
```

Propagation uses an adaptive exponential-midpoint, Lie-group style integrator.
The adaptive controller never steps across path breakpoints, and its error
metric is insensitive to physically irrelevant global Jones phase.

📖 **[Documentation](https://brittonlab.github.io/BIFROST/stable/)**

## Installation

1. Install Juliaup, the julia version manager:

   ```bash
   curl -fsSL https://install.julialang.org | sh
   ```

   Follow the on-screen instructions to add Juliaup to your PATH.

2. Install julia 1.11 and set it as the system-wide default:

   ```bash
   juliaup add 1.11
   juliaup default 1.11
   ```

3. Setup the julia environment for BIFROST.
    ```bash
    cd bifrost
    julia --project=. -e "using Pkg; Pkg.instantiate()"
    ```
    You subsequently activate the environment using `julia --project=.`

    The resulting julia environment consists of the following.
    ```
    bifrost/
    ├── Project.toml    # declared dependencies
    └──  Manifest.toml  # exact dependency graph
    ~/.julia/           # global package cache
    ```

3. (optional) Setup python support using `uv`.  This supports the pythonic API for BIFROST.
    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    cd bifrost
    uv sync
    ```

    The resulting python environment consists of the following.
    ```
    bifrost/
    ├── pyproject.toml  # declared dependencies
    ├── uv.lock         # exact resolved versions
    └── .venv/          # python environment
    ```

Bifrost is a julia module. Add it to source files with `using Bifrost`. Plotting functionality
is segregated in a second module called `Bifrost.Plots`.

## Quick Start

From the repository root, run the test suite:

```bash
julia --project=. test/runtests.jl
```

Human-inspected visual demos live in the notebook
`test/human/bifrost-demos.ipynb` (its §1 is the smallest end-to-end example); the
notebook activates the environment pinned by `test/human/Project.toml`.

## Documentation

Full documentation is built with [Documenter.jl](https://documenter.juliadocs.org)
from the `docs/` tree. 

Build it locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

then open `docs/build/index.html`.
