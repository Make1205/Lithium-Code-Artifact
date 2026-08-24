#!/usr/bin/env sh
set -eu

DUDECT_N="${DUDECT_N:-200000}"
OUT="${1:-test/dudect_all.txt}"
BACKENDS="${BACKENDS:-ref avx2 avx512}"
MODES="${MODES:-120 180 260}"

mkdir -p "$(dirname "$OUT")"
: > "$OUT"

log() {
  printf '%s\n' "$*" | tee -a "$OUT"
}

log "=== Lithium-Code dudect-sign all backends ==="
log "date: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
log "host: $(uname -sr)"
log "DUDECT_N: $DUDECT_N"
log "BACKENDS: $BACKENDS"
log "MODES: $MODES"
log ""

for backend in $BACKENDS; do
  for mode in $MODES; do
    log "--- backend=$backend mode=$mode start=$(date -u '+%Y-%m-%dT%H:%M:%SZ') ---"
    (
      cd "$backend"
      make -j1 -s MODES="$mode" DUDECT_N="$DUDECT_N" dudect-sign
    ) 2>&1 | tee -a "$OUT"
    log "--- backend=$backend mode=$mode end=$(date -u '+%Y-%m-%dT%H:%M:%SZ') ---"
    log ""
  done
done

log "=== done ==="
log "wrote: $OUT"
