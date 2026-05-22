"""
wrapper.py — minimal Python interface to the Bifrost Julia package.

This file exists so that ``import wrapper as bf`` works inside Python.

Once imported, every export of the `Bifrost` Julia module — and every export
of its sub-modules surfaced at top level (`FiberCrossSection`, `Fiber`,
`PathSpecBuilder`, `propagate_fiber`, …) — is accessible as `bf.<name>`.
When new files or functions are added to `Bifrost.jl`, they appear here
automatically the next time you import `wrapper`. No edits required.

Quick start
-----------
    pip install juliacall numpy
    # then, from inside the BIFROST repo (or set BIFROST_JL_PATH):
    import wrapper as bf

    xs = bf.FiberCrossSection(
        bf.GermaniaSilicaGlass(0.036),
        bf.GermaniaSilicaGlass(0.0),
        8.2e-6, 125e-6,
    )
    spec = bf.PathSpecBuilder()
    bf.straight_b(spec, length=0.5)      # `straight!` in Julia → `straight_b` in Python
    bf.bend_b(spec, radius=0.05, angle=3.14159 / 2)
    bf.straight_b(spec, length=0.5)

    fiber = bf.Fiber(bf.build(spec), cross_section=xs, T_ref_K=297.15)
    J, stats = bf.propagate_fiber(fiber, λ_m=1550e-9)
    J, G, _  = bf.propagate_fiber_sensitivity(fiber, λ_m=1550e-9)
    print(bf.output_dgd_2x2(J, G) * 1e12, "ps")
"""

from __future__ import annotations

import os
from pathlib import Path

try:
    from juliacall import Main as jl
except ImportError as e:
    raise ImportError(
        "wrapper.py requires `juliacall`. Install with: pip install juliacall"
    ) from e


# Locate the BIFROST repository
def _autodiscover_repo() -> Path:
    here = Path(__file__).resolve().parent
    for d in [here, *here.parents]:
        proj = d / "Project.toml"
        if proj.is_file() and 'name = "Bifrost"' in proj.read_text():
            return d
    raise FileNotFoundError(
        "Could not locate the Bifrost repo (no Project.toml with "
        'name = "Bifrost" found when walking up from this file). '
        "Set BIFROST_JL_PATH to the repo root."
    )


_REPO = Path(os.environ.get("BIFROST_JL_PATH") or _autodiscover_repo())


# Bootstrap Julia + the Bifrost package
jl.seval("using Pkg")
jl.seval(f'Pkg.activate(raw"{_REPO}")')
jl.seval("Pkg.instantiate()")
jl.seval("using Bifrost")


# ────────────────────────────────────────────────────────────────────────────
# Optionally include the nonlinear extension (Raman/Brillouin), which lives
# in `src/nonlinear/` which isn't yet wired into the Bifrost module proper. 
# Its three dependencies are installed on demand. THIS IS TEMPORARY.
# ────────────────────────────────────────────────────────────────────────────
_NONLIN = _REPO / "src" / "nonlinear"
has_nonlinear = False
if _NONLIN.is_dir():
    for pkg in ("DataInterpolations", "FFTW", "SpecialFunctions"):
        jl.seval(f'try; using {pkg}; catch; Pkg.add("{pkg}"); using {pkg}; end')
    for fname in ("raman.jl", "brillouin.jl"):
        path = _NONLIN / fname
        if path.is_file():
            jl.include(str(path))
            has_nonlinear = True

# Attribute forwarding (PEP 562). `bf.X` resolves to `Bifrost.X` first, then to
# `Main.X`. Adding new names on the Julia side requires no changes to this file.
def __getattr__(name: str):
    try:
        return getattr(jl.Bifrost, name)
    except AttributeError:
        pass
    try:
        return getattr(jl, name)
    except AttributeError as e:
        raise AttributeError(
            f"wrapper: '{name}' not defined in Bifrost or Main "
            f"(the running Julia session)."
        ) from e


def __dir__():
    """Expose names from Bifrost (and our own module-level bindings) for tab-completion."""
    own = {"jl", "has_nonlinear", "info"}
    try:
        return sorted(own | set(dir(jl.Bifrost)))
    except Exception:
        return sorted(own)


def info() -> None:
    """One-line diagnostic of what the wrapper loaded."""
    print(f"BIFROST repo:     {_REPO}")
    print(f"Bifrost loaded:   {bool(jl.seval('isdefined(Main, :Bifrost)'))}")
    print(f"has_nonlinear:    {has_nonlinear}")
    print(f"Top-level exports: {len(dir(jl.Bifrost))} names "
          f"(use `dir(wrapper)` for the full list)")
