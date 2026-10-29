#!/usr/bin/env bash
# CCID-Epigenome environment bootstrap (WSL / Linux)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python3 -m venv .venv 2>/dev/null || true
source .venv/bin/activate
pip install --upgrade pip --quiet
pip install -r requirements.txt --quiet
export PYTHONPATH="$PWD/src"
pytest -q
echo
echo "Setup complete. Expect '72 passed'."
echo "Next: read START_HERE.md section 3 to regenerate the FAZ 0 inputs."
