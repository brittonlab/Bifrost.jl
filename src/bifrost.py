"""Boot juliacall against the BIFROST Julia project.

Call `start()` once per process before any other juliacall use. Sets the env
vars juliacall requires (PYTHON_JULIAPKG_EXE so juliacall uses the system
Julia on PATH, PYTHON_JULIACALL_HANDLE_SIGNALS, and optionally
PYTHON_JULIACALL_THREADS), activates and instantiates the repo's Project.toml,
and returns the Julia `Main` handle.

This file lives at `src/bifrost.py` next to the Julia sources so the wheel can
co-locate the juliacall entry point with the repo it activates. The default
project resolves to `parents[1]` of this file, which means an editable install
(`uv sync` / `pip install -e .` from the repo root) just works. Users who
install the wheel into an unrelated environment must pass `project=` to
`start()` to point at the Julia project they want to drive.
"""

import os
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

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
    threads: str | int | None = None,
    instantiate: bool = True,
    project: str | os.PathLike | None = None,
):
    """Boot juliacall against this repo and return the Julia `Main` handle.

    Idempotent: the first call boots Julia; subsequent calls return the cached
    handle. A second call with different arguments raises RuntimeError, since
    env vars from the first call have already taken effect and re-activating
    against a different project would be a foot-gun.

    Parameters
    ----------
    threads : str | int | None
        Value for PYTHON_JULIACALL_THREADS. None leaves it unset (Julia
        starts single-threaded). Pass "auto" or an int >= 1 to enable threads.
        Must be set before juliacall is imported, so this argument is honored
        only on the first call in a process.
    instantiate : bool
        Run `Pkg.instantiate()` after activation. Set False when the
        environment is known to be resolved already.
    project : str | Path | None
        Julia project to activate. Defaults to the repo root resolved
        from this file's location, which is correct for an editable install.
        Pass an explicit path when driving a different Julia project (e.g.
        when this package is wheel-installed outside the repo).

    Returns
    -------
    jl : juliacall.Main
        The Julia Main module handle, ready to use.

    Raises
    ------
    RuntimeError
        If start() was already called with different arguments.
    FileNotFoundError
        If the project path doesn't contain a Project.toml.
    """

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
    return jl

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
    """Load the optional Bifrost.plots module."""
    if _jl is None:
        start()
    
    try:
        _jl.seval("using Bifrost.plots")
    except Exception as e:
        raise RuntimeError(
            "Failed to load Bifrost.plots. Is it installed? "
            "Check your Bifrost environment or install the plots extra."
        ) from e

# ────────────────────────────────────────────────────────────────────────────
# PEP 562: Module-level __getattr__ and __dir__ for dynamic attribute access.
# This allows `import bifrost as bf; bf.some_julia_name(...)` to work naturally.
# On first attribute access, we auto-start Julia if it hasn't been started yet.
# ────────────────────────────────────────────────────────────────────────────

def __getattr__(name: str):
    """Dynamically forward attribute access to the Julia Bifrost module.
    
    Auto-starts the Julia environment on first access if not already started.
    """
    global _jl
    if _jl is None:
        # First access: auto-start with defaults
        _jl = start()
    
    # Try Bifrost module first, then Main
    try:
        return getattr(_jl.Bifrost, name)
    except AttributeError:
        pass
    try:
        return getattr(_jl, name)
    except AttributeError as e:
        raise AttributeError(
            f"bifrost: '{name}' not defined in Bifrost or Main "
            f"(the running Julia session)."
            f"If this is from Bifrost.plots, try bifrost.load_plots() first."
        ) from e


def __dir__():
    """Expose names from Bifrost for tab-completion."""
    if _jl is None:
        return ["start", "info"]
    
    own = {"start", "info"}
    try:
        return sorted(own | set(dir(_jl.Bifrost)))
    except Exception:
        return sorted(own)