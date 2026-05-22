"""MCM-under-threads probe via juliacall, with a native-Julia parity check.

Runs MCMDemo.run() under juliacall (Julia started multi-threaded), once
serially and once threaded, asserts the two agree, writes the threaded
result to output/julia-call-mcm.python.csv, then shells out to native
`julia` to produce output/julia-call-mcm.julia.csv from the same module
and compares the two byte-for-byte.

See README.md."""

import filecmp
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "output"
OUT.mkdir(exist_ok=True)
HERE = Path(__file__).resolve().parent

_julia_exe = shutil.which("julia")
if _julia_exe is None:
    sys.exit("No `julia` found on PATH. Install Julia (e.g. via juliaup).")

os.environ.setdefault("PYTHON_JULIAPKG_EXE", _julia_exe)
os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
os.environ.setdefault("PYTHON_JULIACALL_THREADS", "auto")

from juliacall import Main as jl  # noqa: E402

jl.seval(f'import Pkg; Pkg.activate(raw"{REPO}"); Pkg.instantiate()')

n_threads = int(jl.seval("Threads.nthreads()"))
print(f"Julia threads (juliacall) = {n_threads}")
assert n_threads > 1, (
    f"juliacall started Julia with {n_threads} thread(s); "
    "set PYTHON_JULIACALL_THREADS=auto (or an integer > 1) "
    "to exercise threading."
)

jl.seval(f'include(raw"{HERE / "julia-call-mcm.jl"}")')
MCMDemo = jl.MCMDemo

T_SIGMA = 10.0  # K; keep python and native runners in sync
rows_serial = MCMDemo.run(threaded=False, T_sigma=T_SIGMA)
rows_threaded = MCMDemo.run(threaded=True, T_sigma=T_SIGMA)

# Exact scalar comparison: same seed, same input order, same arithmetic.
def to_tuples(rows):
    return [(int(r.i), float(r.T_mean), float(r.T_std),
             float(r.n_mean), float(r.n_std)) for r in rows]

s = to_tuples(rows_serial)
t = to_tuples(rows_threaded)
assert s == t, f"Serial vs threaded mismatch:\n  serial={s}\n  threaded={t}"
print(f"Serial vs threaded (same Julia) match: True ({len(s)} rows)")

py_csv = OUT / "julia-call-mcm.python.csv"
MCMDemo.write_table(str(py_csv), rows_threaded)
print(f"Wrote {py_csv}")

jl_csv = OUT / "julia-call-mcm.julia.csv"
proc = subprocess.run(
    [
        _julia_exe,
        f"--project={REPO}",
        f"--threads={n_threads}",
        str(HERE / "julia-call-mcm-native.jl"),
        str(T_SIGMA),
    ],
    cwd=REPO,
    check=False,
)
if proc.returncode != 0:
    sys.exit(f"native Julia runner failed (rc={proc.returncode})")

if filecmp.cmp(py_csv, jl_csv, shallow=False):
    print(f"Parity OK: {py_csv.name} == {jl_csv.name}")
else:
    # Fall back to a numerical diff so the user sees what diverged.
    import csv

    def read(path):
        with open(path, newline="") as f:
            return list(csv.reader(f))

    py_rows = read(py_csv)
    jl_rows = read(jl_csv)
    worst = 0.0
    for a, b in zip(py_rows[1:], jl_rows[1:]):
        for av, bv in zip(a[1:], b[1:]):
            worst = max(worst, abs(float(av) - float(bv)))
    sys.exit(
        f"PARITY FAIL: {py_csv} != {jl_csv} "
        f"(worst |Δ| = {worst:g}; see files for full diff)"
    )
