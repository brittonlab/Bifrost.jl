"""Boot juliacall against the BIFROST Julia project.

Call `start()` once per process before any other juliacall use. Sets the env
vars juliacall requires (PYTHON_JULIAPKG_EXE so juliacall uses the system
Julia on PATH, PYTHON_JULIACALL_HANDLE_SIGNALS, and optionally
PYTHON_JULIACALL_THREADS), activates and instantiates the repo's Project.toml,
and returns the Julia `Main` handle.
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


def start(*, threads: str | int | None = None, instantiate: bool = True):
    """Boot juliacall against this repo and return the Julia `Main` handle.

    Idempotent: the first call boots Julia; subsequent calls return the cached
    handle. A second call with different `threads=` or `instantiate=` raises
    RuntimeError, since env vars from the first call have already taken effect
    and re-activating against a different project would be a foot-gun.

    threads: value for PYTHON_JULIACALL_THREADS. None leaves it unset (Julia
        starts single-threaded). Pass "auto" or an int >= 1 to enable threads.
        Must be set before juliacall is imported, so this argument is honored
        only on the first call in a process.
    instantiate: run `Pkg.instantiate()` after activation. Set False when the
        environment is known to be resolved already.
    """
    global _jl, _start_kwargs
    kwargs = {"threads": threads, "instantiate": instantiate}
    if _jl is not None:
        if kwargs != _start_kwargs:
            raise RuntimeError(
                f"bifrost.start already called with {_start_kwargs}; "
                f"second call with {kwargs} would be ignored. "
                "Call start() once per process."
            )
        return _jl

    julia_exe()
    os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
    if threads is not None:
        os.environ.setdefault("PYTHON_JULIACALL_THREADS", str(threads))

    from juliacall import Main as jl

    jl.seval(f'import Pkg; Pkg.activate(raw"{REPO}")')
    if instantiate:
        jl.seval("import Pkg; Pkg.instantiate()")

    _jl = jl
    _start_kwargs = kwargs
    return jl
