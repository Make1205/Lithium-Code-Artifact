#!/bin/sh
# Lithium signing variable-latency (KyberSlash/TIMECOP) smoke for one
# backend + parameter mode.
#
# Marks signing secret-key material and signing randomness as secret/undefined
# (ref/test/timecop_smoke.c, shared by every backend via -I../ref/test) and runs
# one keygen+sign+verify under the patched variable-latency Memcheck. The GATE is
# the KyberSlash property: FAIL iff any secret operand reaches a variable-latency
# instruction ("Variable-latency instruction operand ..."). Secret-dependent
# control flow is out of scope (owned by ct_scan_matrix.sh + dudect-sign); its
# intentional constant-time cmov/select idioms (Memcheck:Cond) are suppressed via
# lithium-ct.supp.
#
# Backends (--backend):
#   ref     scalar reference  -- always runnable under Valgrind.
#   avx2    AVX2 SIMD         -- runnable (pure AVX2, built with -mno-avx512f).
#   avx512  AVX-512 SIMD      -- Valgrind 3.23 cannot decode the EVEX/AVX-512
#           instructions and dies with SIGILL, so this backend is reported as
#           SKIP (exit 3), not FAIL. Its variable-latency property is covered
#           statically by ct_scan_matrix.sh + dudect-sign.
#
# The per-backend source lists mirror each backend Makefile's `test_sign` list,
# minus randombytes.c (the harness defines randombytes()); drift surfaces as a
# loud compile error.
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
SUPP=${TIMECOP_SUPP:-$ROOT/tools/timecop/lithium-ct.supp}
HARNESS=$ROOT/ref/test/timecop_smoke.c
VG=
BACKEND=ref
MODE=120
while [ $# -gt 0 ]; do
  case "$1" in
    --valgrind) VG=$2; shift 2;;
    --backend) BACKEND=$2; shift 2;;
    --mode) MODE=$2; shift 2;;
    *) echo "usage: $0 --valgrind PATH [--backend ref|avx2|avx512] [--mode 120|180|260]" >&2; exit 2;;
  esac
done
[ -n "$VG" ] || { echo "usage: $0 --valgrind PATH [--backend ref|avx2|avx512] [--mode 120|180|260]" >&2; exit 2; }
# Resolve VG and its bundled headers to ABSOLUTE paths before we cd into a backend.
case "$VG" in /*) ;; *) VG=$(CDPATH= cd -- "$(dirname "$VG")" && pwd)/$(basename "$VG");; esac
VG_INCLUDE=${VALGRIND_INCLUDE:-$(CDPATH= cd -- "$(dirname "$VG")/../include" 2>/dev/null && pwd || true)}
VG_INC_FLAG=
[ -n "$VG_INCLUDE" ] && [ -d "$VG_INCLUDE/valgrind" ] && VG_INC_FLAG="-I$VG_INCLUDE"

# Per-mode vendored NTT directory (shared by avx2/avx512).
case "$MODE" in
  120) qs=q57601n256;;
  180) qs=q30977n256;;
  260) qs=q64513n512;;
  *) echo "unknown --mode $MODE (expected 120|180|260)" >&2; exit 2;;
esac

cc=${CC:-gcc}
case "$BACKEND" in
  ref)
    DIR=$ROOT/ref
    CF="-O2 -g -fno-inline -fno-omit-frame-pointer -Wall -Wextra -std=c99 $VG_INC_FLAG -I. -Itest"
    SRCS="ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c test/cpucycles.c"
    ;;
  avx2)
    DIR=$ROOT/avx2
    CF="-O3 -g -fno-omit-frame-pointer -std=c99 $VG_INC_FLAG -I. -I../ref -I../ref/test -mavx2 -mbmi2 -mpopcnt -march=native -mtune=native -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X"
    SRCS="ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $qs/ntt.S $qs/ntt_consts.c test/cpucycles.c"
    ;;
  avx512)
    DIR=$ROOT/avx512
    CF="-O3 -g -fno-omit-frame-pointer -std=c99 $VG_INC_FLAG -I. -I../ref -I../ref/test -mavx512f -mavx512bw -mavx512vl -mavx512dq -mavx512ifma -mbmi2 -mpopcnt -march=native -mtune=native -DUSE_AVX512_NTT -DUSE_AVX512_SAMPLER -DUSE_AVX512_SHAKE8X"
    SRCS="ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x8.c keccakf1600x8.c $qs/ntt_avx512.S $qs/ntt_consts_avx512.c test/cpucycles.c"
    ;;
  *) echo "unknown --backend $BACKEND (expected ref|avx2|avx512)" >&2; exit 2;;
esac

cd "$DIR"
mkdir -p out ../docs/design-notes
bin=out/timecop_smoke_${BACKEND}_$MODE
# shellcheck disable=SC2086
$cc $CF -DLITHIUM_MODE="$MODE" "$HARNESS" $SRCS -lm -o "$bin"

# Full (unsuppressed) run is the audit artifact; the gate run applies the
# constant-time cmov suppression so stderr isolates the variable-latency signal.
set +e
"$VG" -q --variable-latency-errors=yes "./$bin" \
  > "$bin.full.out" 2> "$bin.full.err"
"$VG" -q --variable-latency-errors=yes --suppressions="$SUPP" "./$bin" \
  > "$bin.out" 2> "$bin.err"
rc=$?
set -e

# Valgrind cannot execute some SIMD (notably AVX-512): it aborts the guest with
# SIGILL on an instruction it cannot decode. Report SKIP (exit 3), not FAIL.
if grep -qiE 'SIGILL|Illegal opcode|unhandled instruction|disInstr' "$bin.full.err"; then
  why=$(grep -m1 -iE 'Illegal opcode|unhandled instruction|disInstr' "$bin.full.err" | sed 's/^==[0-9]*== *//')
  echo "timecop $BACKEND mode-$MODE smoke: SKIP — Valgrind cannot execute this backend (${why:-SIGILL}); covered statically by ct_scan + dudect-sign"
  exit 3
fi

varlat=$(grep -c 'Variable-latency instruction operand' "$bin.err" || true)
cond_full=$(grep -c 'Conditional jump or move depends' "$bin.full.err" || true)
value_full=$(grep -c 'Use of uninitialised value' "$bin.full.err" || true)

if [ "$varlat" -ne 0 ]; then
  echo "FAIL: Lithium $BACKEND mode-$MODE smoke reported variable-latency secret operand" >&2
  grep -A8 'Variable-latency instruction operand' "$bin.err" >&2
  exit 1
fi
if [ "$rc" -ne 0 ]; then
  echo "FAIL: $BACKEND mode-$MODE smoke binary exited rc=$rc" >&2
  cat "$bin.out" >&2
  cat "$bin.err" >&2
  exit 1
fi
{
  echo
  echo "## Lithium $BACKEND mode-$MODE Smoke"
  echo
  echo "- Valgrind: \`$VG\`"
  echo "- Binary: \`$BACKEND/$bin\` (\`$(cat "$bin.out")\`)"
  echo "- Flags: \`--variable-latency-errors=yes --suppressions=tools/timecop/lithium-ct.supp\`"
  echo "- Variable-latency secret operands: **$varlat** (gate: must be 0) -> PASS"
  echo "- Out-of-scope intentional-CT reports in the unsuppressed run: $cond_full \`Memcheck:Cond\`, $value_full \`Memcheck:Value\` (declassified output + branchless CT idioms; owned by ct_scan + dudect-sign)."
} >> ../docs/design-notes/TIMECOP.md
echo "timecop $BACKEND mode-$MODE smoke: PASS (variable-latency operands=$varlat; out-of-scope cond=$cond_full value=$value_full)"
