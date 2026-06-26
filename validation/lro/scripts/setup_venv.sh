#!/usr/bin/env bash
# Creates a Python virtual environment with all packages needed to run
# the NEPTUNE validation scripts (lro_validation.py, plot_kepler_elements.py,
# plot_variance.py).
#
# Usage:
#   ./setup_venv.sh            # creates .venv in the current directory
#   ./setup_venv.sh /some/path # creates the venv at /some/path

set -euo pipefail

VENV_DIR="${1:-.venv}"

# Prefer python3.11+ if available (spiceypy wheels exist for 3.8-3.13).
PYTHON=$(command -v python3 2>/dev/null || command -v python 2>/dev/null)

if [[ -z "$PYTHON" ]]; then
    echo "ERROR: No Python 3 interpreter found on PATH." >&2
    exit 1
fi

PY_VERSION=$("$PYTHON" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "Using $PYTHON  (Python $PY_VERSION)"

# Create venv
if [[ -d "$VENV_DIR" ]]; then
    echo "Virtual environment already exists at '$VENV_DIR' — skipping creation."
else
    echo "Creating virtual environment at '$VENV_DIR' ..."
    "$PYTHON" -m venv "$VENV_DIR"
fi

# Activate
# shellcheck source=/dev/null
source "$VENV_DIR/bin/activate"

echo "Upgrading pip ..."
pip install --quiet --upgrade pip

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Installing packages ..."
pip install --quiet -r "$SCRIPT_DIR/requirements.txt"

echo ""
echo "Done.  Activate the environment with:"
echo "  source $VENV_DIR/bin/activate"
