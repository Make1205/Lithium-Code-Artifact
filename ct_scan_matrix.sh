#!/bin/sh
# Full machine-code constant-time scan matrix.
set -eu
cd "$(dirname "$0")"
REPORT="${CT_SCAN_REPORT:-test/ct_scan_matrix.txt}"
exec python3 tools/ct_scan.py \
  --backends "${CT_BACKENDS:-all}" \
  --modes "${CT_MODES:-all}" \
  --ccs "${CT_CCS:-gcc clang}" \
  --opts="${CT_OPTS:--O3 -Os}" \
  --include-optional \
  --report "$REPORT"
