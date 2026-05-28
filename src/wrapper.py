"""Tab-completable wrappers around Julia modules via juliacall.

The model: Python is the workspace; Julia modules are guest libraries. Each
call to `wrap(name, [path])` returns a `JuliaModule` whose `dir()` and tab-
completion list only that module's `export`ed names — not the ~1000 names a
raw juliacall module inherits from `Base`.

Two modes:

- `wrap("Bifrost")` — module is available via the project's environment
  (Project.toml dependency, or already loaded by something else). We run
  `using Bifrost` once to ensure it's in scope, then fetch it.
- `wrap("Zoo", "docs/juliacall-demo.jl")` — ad-hoc module. The path is taken
  as absolute if absolute, otherwise resolved against the repo root. We
  `jl.include(...)` it, then fetch `jl.Zoo`.

No print side effects. Errors propagate unmodified from juliacall / Julia.

Python is the workspace; Julia is a guest. The booted `Main` handle is kept
private — only names exported by a wrapped module are reachable from Python.
"""

import os
import shutil
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]

_jl: Any = None
_julia_exe: str | None = None


def julia_exe() -> str:
    """Return the Julia executable selected for this Python process.

    Raises FileNotFoundError if no `julia` is on PATH and PYTHON_JULIAPKG_EXE
    isn't set — callers decide whether to exit.
    """
    global _julia_exe
    if _julia_exe is not None:
        return _julia_exe

    julia = os.environ.get("PYTHON_JULIAPKG_EXE")
    if julia is None:
        julia = shutil.which("julia")
        if julia is None:
            raise FileNotFoundError(
                "No `julia` found on PATH and PYTHON_JULIAPKG_EXE is not set. "
                "Install Julia (e.g. via juliaup)."
            )
        os.environ["PYTHON_JULIAPKG_EXE"] = julia

    _julia_exe = julia
    return _julia_exe


def _boot() -> Any:
    """Boot juliacall against this repo, returning the `Main` handle.

    Idempotent: caches the handle in `_jl`. Threads must be configured via the
    PYTHON_JULIACALL_THREADS env var before the first call; we don't expose a
    threads= argument because `wrap` is the primary entry point and threading
    is a process-wide property best set outside the call.
    """
    global _jl
    if _jl is not None:
        return _jl

    julia_exe()
    os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")

    from juliacall import Main as jl

    jl.seval(f'import Pkg; Pkg.activate(raw"{REPO}")')
    jl.seval("import Pkg; Pkg.instantiate()")

    _jl = jl
    return jl


def wrap(name: str, path: "str | os.PathLike[str] | None" = None) -> "JuliaModule":
    """Return a tab-completable wrapper around the Julia module named `name`.

    `wrap("Bifrost")` — no path. Runs `using <name>` then fetches `Main.<name>`.
        Works for any module the active project exposes (Project.toml dependency
        or otherwise reachable from `Main`).

    `wrap("Zoo", "docs/juliacall-demo.jl")` — with path. Absolute paths are
        used as-is; relative paths are resolved against the repo root. The file
        is `include`d into `Main` (defining `Zoo` there), then `Main.Zoo` is
        fetched. Re-runs the include on every call, so editing the .jl file and
        re-wrapping in a REPL reloads it.

    Errors propagate unmodified: a missing module raises whatever juliacall
    raises; a missing file raises FileNotFoundError from `Path.resolve(strict)`.
    """
    j = _boot()
    if path is not None:
        p = Path(path)
        if not p.is_absolute():
            p = REPO / p
        p = p.resolve(strict=True)
        j.include(str(p))
    else:
        j.seval(f"using {name}")
    return JuliaModule(getattr(j, name))


class JuliaModule:
    """Tab-completable view over a juliacall module, exposing only its exports.

    `dir()` on a raw juliacall module reflects all of `names(mod; all=true,
    imported=true)`, which inherits ~1000 names from `Base`. This wrapper's
    `__dir__` returns only the module's own `export`ed names. Attribute access
    forwards to the underlying module unchanged.

    The raw juliacall handle is available as `.jl_module` for escape hatches.
    """

    def __init__(self, module: Any):
        # Bypass our own __setattr__-by-default so these aren't mistaken for
        # exports during attribute lookup.
        object.__setattr__(self, "jl_module", module)
        names = [str(n) for n in _jl.names(module)]
        self_name = str(_jl.nameof(module))
        exports = sorted(n for n in names if n != self_name)
        object.__setattr__(self, "_exports", exports)

    def __getattr__(self, name: str) -> Any:
        # Called only when normal lookup fails, so it never shadows _exports
        # or jl_module. Raises straight through from juliacall on miss.
        return getattr(self.jl_module, name)

    def __dir__(self):
        return self._exports

    def __repr__(self) -> str:
        mod_name = str(_jl.nameof(self.jl_module))
        return f"<JuliaModule {mod_name}: {len(self._exports)} exports>"
