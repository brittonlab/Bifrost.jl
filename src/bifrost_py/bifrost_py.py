"""Boot juliacall against the BIFROST Julia project.

[... your existing docstring ...]
"""

import os
import shutil
import sys
from pathlib import Path
from typing import Optional

REPO = Path(__file__).resolve().parents[2]

_jl = None
_start_kwargs = None
_julia_exe = None


def julia_exe() -> str:
    """Return the Julia executable selected for this Python process."""
    global _julia_exe
    if _julia_exe is not None:
        return _julia_exe

    julia = os.environ.get("PYTHON_JULIAPKG_EXE")
    if julia is None:
        julia = shutil.which("julia")
        if julia is None:
            sys.exit("No `julia` found on PATH. Install Julia (e.g. via juliaup).")
        os.environ["PYTHON_JULIAPKG_EXE"] = julia

    _julia_exe = julia
    return _julia_exe


def start(
    *,
    threads: Optional[str | int] = None,
    instantiate: bool = True,
    project: Optional[str | os.PathLike] = None,
):
    """Boot juliacall against this repo and return the Julia `Main` handle.
    
    [... your existing docstring ...]
    """
    print("Starting Julia environment...", end=" ")
    global _jl, _start_kwargs
    project_path = Path(project).resolve() if project is not None else REPO
    kwargs = {
        "threads": threads,
        "instantiate": instantiate,
        "project": str(project_path),
    }
    if _jl is not None:
        if kwargs != _start_kwargs:
            raise RuntimeError(
                f"bifrost.start already called with {_start_kwargs}; "
                f"second call with {kwargs} would be ignored. "
                "Call start() once per process."
            )
        return _jl

    if not (project_path / "Project.toml").is_file():
        sys.exit(f"No Project.toml found at {project_path}. Pass project=<path>.")

    julia_exe()
    os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
    if threads is not None:
        os.environ.setdefault("PYTHON_JULIACALL_THREADS", str(threads))

    from juliacall import Main as jl

    jl.seval(f'import Pkg; Pkg.activate(raw"{project_path}")')
    if instantiate:
        jl.seval("import Pkg; Pkg.instantiate()")

    jl.seval("using Bifrost")

    _jl = jl
    _start_kwargs = kwargs
    print("Complete.")
    return _jl


def get_jl():
    """Get the Julia Main handle, auto-starting if needed."""
    global _jl
    if _jl is None:
        _jl = start()
    return _jl


def info() -> None:
    """Print diagnostic info about the loaded environment."""
    if _jl is None:
        print("Julia environment not yet started. Call bifrost.start() or access")
        print("an attribute (which will trigger auto-start).")
        return
    
    proj = _start_kwargs["project"]
    print(f"BIFROST project:  {proj}")
    print(f"Julia exe:        {julia_exe()}")
    print(f"Bifrost loaded:   {bool(_jl.seval('isdefined(Main, :Bifrost)'))}")


def load_plots():
    """Load the optional Bifrost.Plots module."""
    jl = get_jl()
    try:
        jl.seval("using Bifrost.Plots")
    except Exception as e:
        raise RuntimeError(
            "Failed to load Bifrost.Plots. Is it installed? "
            "Check your Bifrost environment."
        ) from e


def __getattr__(name: str):
    """Dynamically forward attribute access to the Julia Bifrost module."""
    jl = get_jl()
    
    try:
        return getattr(jl.Bifrost, name)
    except AttributeError:
        pass
    try:
        return getattr(jl, name)
    except AttributeError as e:
        raise AttributeError(
            f"bifrost: '{name}' not defined in Bifrost or Main. "
            f"If this is from Bifrost.Plots, try bifrost.load_plots() first."
        ) from e


def __dir__():
    """Expose names from Bifrost for tab-completion."""
    if _jl is None:
        return ["start", "info", "load_plots"]
    
    own = {"start", "info", "load_plots"}
    try:
        return sorted(own | set(dir(_jl.Bifrost)))
    except Exception:
        return sorted(own)