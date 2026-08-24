#!/bin/sh
# One-shot TIMECOP driver: ensure a patched variable-latency Valgrind exists,
# then run the positive-control secret-division probe and the Lithium signing
# smoke for each runnable backend x parameter mode. Optional dynamic-taint
# evidence; NOT part of the hard gate (it needs a local patched-Valgrind build +
# network on first run).
#
#   tools/timecop/run_all.sh                       # build-if-needed, then probes
#   TIMECOP_BACKENDS="ref avx2" TIMECOP_MODES="120 180 260" tools/timecop/run_all.sh
#   TIMECOP_TRY_AVX512=1 tools/timecop/run_all.sh   # attempt avx512 (will SKIP: SIGILL)
#   TIMECOP_VALGRIND=/path/to/valgrind tools/timecop/run_all.sh
#
# Backends: ref and avx2 run under Valgrind; avx512 cannot (Valgrind 3.23 SIGILLs
# on the EVEX/AVX-512 instructions), so it is reported as SKIP and its
# variable-latency property is covered statically by ct_scan + dudect-sign.
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
WORK=${TIMECOP_WORKDIR:-$ROOT/external/timecop/valgrind-varlat}
VG=${TIMECOP_VALGRIND:-$WORK/install/bin/valgrind}

has_varlat() {
  [ -x "$1" ] && "$1" --help 2>&1 | grep -q -- '--variable-latency-errors'
}

if ! has_varlat "$VG"; then
  echo "== patched Valgrind not found; building (one-time, downloads + compiles) =="
  JOBS=${JOBS:-$(nproc 2>/dev/null || echo 2)} sh "$ROOT/tools/timecop/build_valgrind_varlat.sh"
  VG=$WORK/install/bin/valgrind
fi
has_varlat "$VG" || { echo "FAIL: no usable variable-latency Valgrind at $VG" >&2; exit 1; }

echo "== TIMECOP Valgrind: $("$VG" --version) ($VG) =="
sh "$ROOT/tools/timecop/prove_varlat.sh" --valgrind "$VG"
for be in ${TIMECOP_BACKENDS:-ref avx2}; do
  for mode in ${TIMECOP_MODES:-120 180 260}; do
    sh "$ROOT/tools/timecop/run_sign_smoke.sh" --valgrind "$VG" --backend "$be" --mode "$mode"
  done
done
# AVX-512 cannot run under Valgrind (SIGILL). Off by default; opt in to demonstrate.
if [ "${TIMECOP_TRY_AVX512:-0}" = 1 ]; then
  for mode in ${TIMECOP_MODES:-120 180 260}; do
    sh "$ROOT/tools/timecop/run_sign_smoke.sh" --valgrind "$VG" --backend avx512 --mode "$mode" \
      || { rc=$?; [ "$rc" = 3 ] || exit "$rc"; }
  done
else
  echo "timecop avx512: SKIP — Valgrind 3.23 cannot execute AVX-512 (SIGILL); covered statically by ct_scan + dudect-sign. Set TIMECOP_TRY_AVX512=1 to demonstrate."
fi
echo "timecop run_all: PASS"
