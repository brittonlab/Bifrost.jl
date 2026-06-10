# Contributing

Read `AGENTS.md`, `CLAUDE.md`, and `ARCHITECTURE.md` before changing this repo. They
describe the agent workflow, project structure, architectural boundaries, tests, and
invariants that should guide implementation.

## Agent Skills

The `julia-docstrings` skill is registered at the claude.ai/code account level and is
available to all sessions on this repository. Use it when creating, revising, or
auditing inline documentation for Julia code.

## Tests

Run the Julia test suite with:

```bash
julia --project=. test/runtests.jl
```

## Documentation

Documentation is built with [Documenter.jl](https://documenter.juliadocs.org). Source
markdown lives under `docs/src/` (Home, Examples, Usage, Theory, and per-module API pages
generated via `@autodocs`); `docs/make.jl` defines the build and navigation. A GitHub
Actions workflow (`.github/workflows/Documenter.yml`) rebuilds and deploys the site to
GitHub Pages.

Build the docs locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The rendered site is written to `docs/build/` (open `docs/build/index.html`).
