#!/bin/bash
# SessionStart hook for Claude Code on the web.
#
# Installs the toolchain BIFROST needs so tests and demos run in a fresh
# remote container: Juliaup + Julia 1.11, all Julia package dependencies
# (228 entries in Manifest.toml), and the optional Python shim used by the
# juliacall wrapper. Mirrors the manual steps in README.md "Installation".
#
# NOTE: this requires the environment's network policy to allow outbound
# access to the Julia and uv download hosts:
#   install.julialang.org, julialang-s3.julialang.org, pkg.julialang.org,
#   cache.julialang.org, astral.sh. Without those, the downloads 403 and the
#   hook cannot install the toolchain.
set -euo pipefail

# Only run in Claude Code on the web (remote) sessions.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"
cd "$PROJECT_DIR"

# ----------------------------------------------------------------------
# 1. Install Juliaup + Julia 1.11 (idempotent).
# ----------------------------------------------------------------------
export JULIAUP_DIR="${JULIAUP_DIR:-$HOME/.juliaup}"
export PATH="$JULIAUP_DIR/bin:$PATH"

if ! command -v juliaup >/dev/null 2>&1; then
  curl -fsSL https://install.julialang.org \
    | sh -s -- --yes --default-channel 1.11
  export PATH="$JULIAUP_DIR/bin:$PATH"
fi

# Ensure 1.11 is present and the default, even if juliaup pre-existed.
juliaup add 1.11 || true
juliaup default 1.11 || true

# ----------------------------------------------------------------------
# 2. Instantiate the Julia project (downloads all package dependencies).
#    Pkg.instantiate resolves against the checked-in Manifest.toml, so it
#    is reproducible and benefits from the cached container layer.
# ----------------------------------------------------------------------
julia --project="$PROJECT_DIR" -e 'using Pkg; Pkg.instantiate()'

# ----------------------------------------------------------------------
# 3. Optional Python shim (juliacall wrapper) via uv. Non-fatal: the Julia
#    test suite does not depend on it, so a Python failure must not block
#    the session.
# ----------------------------------------------------------------------
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh || true
  export PATH="$HOME/.local/bin:$PATH"
fi
if command -v uv >/dev/null 2>&1; then
  uv sync || true
fi

# ----------------------------------------------------------------------
# 4. Persist toolchain on PATH for the rest of the session.
# ----------------------------------------------------------------------
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo "export PATH=\"$JULIAUP_DIR/bin:\$HOME/.local/bin:\$PATH\"" \
    >> "$CLAUDE_ENV_FILE"
fi
