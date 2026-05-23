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

    threads: value for PYTHON_JULIACALL_THREADS. None leaves it unset (Julia
        starts single-threaded). Pass "auto" or an int >= 1 to enable threads.
        Must be set before juliacall is imported, so this argument is honored
        only on the first call in a process.
    instantiate: run `Pkg.instantiate()` after activation. Set False when the
        environment is known to be resolved already.
    project: Julia project to activate. Defaults to the repo root resolved
        from this file's location, which is correct for an editable install.
        Pass an explicit path when driving a different Julia project (e.g.
        when this package is wheel-installed outside the repo).
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

    _jl = jl
    _start_kwargs = kwargs
    return jl
