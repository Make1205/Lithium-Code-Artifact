#!/usr/bin/env bash
set -euo pipefail
export LANG=C LC_ALL=C
cd "$(dirname "$0")/.."
exec python3 bench/run.py
