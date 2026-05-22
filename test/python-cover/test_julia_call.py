"""Entry point for the python-cover tests.

Runs every probe in this folder in order, fails fast on the first error.
Each probe is a standalone script invoked via `uv run python <script>`.

Usage:

    cd test/python-cover && uv run python test_julia_call.py
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

PROBES = [
    "julia-call-demo.py",
    "julia-call-mcm.py",
]


def main() -> int:
    for name in PROBES:
        script = HERE / name
        print(f"\n=== {name} ===")
        proc = subprocess.run(
            ["uv", "run", "python", str(script)],
            cwd=HERE,
        )
        if proc.returncode != 0:
            print(f"\nFAIL: {name} exited with rc={proc.returncode}")
            return proc.returncode
    print("\nAll python-cover probes passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
