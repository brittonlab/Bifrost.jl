#!/bin/bash
# SessionStart hook for Claude Code on the web.
#
# Installs the toolchain BIFROST needs so tests and demos run in a fresh
# remote container: Juliaup + Julia 1.11, all Julia package dependencies
# (228 entries in Manifest.toml), the optional Python shim used by the
# juliacall wrapper, and the GitHub CLI (gh). Mirrors the manual steps in
# README.md "Installation".
#
# NOTE: this requires the remote environment's network egress allowlist to
# permit outbound access to the Julia, uv, and gh download hosts. The
# allowlist is NOT a repo file — set it in the Claude Code on the web UI:
# environment settings -> Network access -> Custom -> Allowed domains
# (one per line; "*." matches any subdomain). Required entries:
#   *.julialang.org            (Juliaup, Julia, Pkg, cache mirrors)
#   *.astral.sh                (uv installer)
#   github.com                 (gh release tarball; in the defaults)
#   objects.githubusercontent.com  (gh release asset redirects; in defaults)
#   cli.github.com             (only if installing gh via apt, not defaults)
# Keep "Also include default list of common package managers" checked.
# Without these, the downloads 403 and the hook cannot install the toolchain.
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

# Persist the toolchain on PATH for the rest of the session NOW, before any
# step that could fail. The harness applies CLAUDE_ENV_FILE to every tool
# shell, including non-interactive ones where ~/.bashrc returns early at its
# "[ -z "$PS1" ] && return" guard and so never reaches juliaup's PATH block.
# Doing this early guarantees `julia` is on PATH even if a later step (e.g.
# Pkg.instantiate) errors out under `set -e`.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo "export PATH=\"$JULIAUP_DIR/bin:\$HOME/.local/bin:\$PATH\"" \
    >> "$CLAUDE_ENV_FILE"
fi

# ----------------------------------------------------------------------
# 2. Instantiate the Julia project (downloads all package dependencies).
#    The checked-in Manifest.toml can drift out of sync with Project.toml
#    (a stale project_hash makes Pkg.instantiate fail with errors like
#    "failed to find source of parent package"). To keep the session robust
#    against that drift, delete the manifest and let Pkg resolve a fresh,
#    self-consistent one against Project.toml before instantiating. The
#    deletion is local to the ephemeral container, so the committed
#    Manifest.toml is untouched.
# ----------------------------------------------------------------------
rm -f "$PROJECT_DIR/Manifest.toml"
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
# 4. GitHub CLI (gh). Used for PRs, issues, and gh api calls. Non-fatal:
#    the test suite does not depend on it. Prefer apt when available
#    (cli.github.com), otherwise fall back to the release tarball
#    (github.com / objects.githubusercontent.com) into ~/.local/bin so no
#    root is required.
# ----------------------------------------------------------------------
export PATH="$HOME/.local/bin:$PATH"
if ! command -v gh >/dev/null 2>&1; then
  if command -v apt-get >/dev/null 2>&1; then
    (
      set -e
      SUDO=""
      [ "$(id -u)" -ne 0 ] && command -v sudo >/dev/null 2>&1 && SUDO="sudo"
      curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
        | $SUDO dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
      $SUDO chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg
      echo "deb [arch=$(dpkg --print-architecture) \
signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] \
https://cli.github.com/packages stable main" \
        | $SUDO tee /etc/apt/sources.list.d/github-cli.list >/dev/null
      $SUDO apt-get update
      $SUDO apt-get install -y gh
    ) || true
  fi
fi
# Fallback: download the release tarball if apt was unavailable or failed.
if ! command -v gh >/dev/null 2>&1; then
  (
    set -e
    GH_VERSION="2.63.2"
    case "$(uname -m)" in
      x86_64|amd64) GH_ARCH="amd64" ;;
      aarch64|arm64) GH_ARCH="arm64" ;;
      *) GH_ARCH="amd64" ;;
    esac
    mkdir -p "$HOME/.local/bin"
    TMP="$(mktemp -d)"
    curl -fsSL \
      "https://github.com/cli/cli/releases/download/v${GH_VERSION}/gh_${GH_VERSION}_linux_${GH_ARCH}.tar.gz" \
      -o "$TMP/gh.tar.gz"
    tar -xzf "$TMP/gh.tar.gz" -C "$TMP"
    install -m 0755 "$TMP/gh_${GH_VERSION}_linux_${GH_ARCH}/bin/gh" \
      "$HOME/.local/bin/gh"
    rm -rf "$TMP"
  ) || true
fi

# ----------------------------------------------------------------------
# 5. PATH persistence already happened right after the Julia install above
#    (so it survives a later failure under `set -e`). Nothing to do here.
# ----------------------------------------------------------------------
