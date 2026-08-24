#!/usr/bin/env bash
#
# run_avxexp.sh -- generate a candidate binary64 (double) polynomial that
# approximates exp(x) on a narrow interval held *away from 0*, using Sollya's
# fpminimax.  This is the floating-point counterpart of run_polyapprox.sh
# (which produces integer-coefficient polynomials via the GALACTICS method).
#
# The "away from 0" interval mirrors the integer ApproxExp range reduction
#     exp(-a) = 2^{-m/B} * exp(r_l),   r_l in [ (B-1)/B*ln2 , ln2 ),
# so the rejection error  |p-f|/|1-f|  does not blow up.  Full rationale and
# error analysis: docs/design-notes/AVXExp.md.
#
# Usage:
#   ./tools/avxexp/run_avxexp.sh                       # defaults below
#   ./tools/avxexp/run_avxexp.sh --dmax 10 --target 43
#   ./tools/avxexp/run_avxexp.sh --lo '15/16*log(2)' --hi 'log(2)'
#   LO='31/32*log(2)' HI='log(2)' DMIN=3 DMAX=14 TARGET=43 ./tools/avxexp/run_avxexp.sh
#
# Defaults:  exp(x) on [31/32*log(2), log(2)], degrees 3..14, target 43 bits.
# Endpoints are Sollya expressions (e.g. log(2), 31/32*log(2)).
#
# Requires Sollya:   sudo apt-get install -y sollya
#               or:  micromamba install -c conda-forge sollya
#
set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLLYA_FILE="$SCRIPT_DIR/AVXPolyApprox.sollya"
LOG_DIR="$SCRIPT_DIR/LogAVXApprox"

# Defaults (overridable by environment, then by --flags).
LO="${LO:-31/32*log(2)}"
HI="${HI:-log(2)}"
DMIN="${DMIN:-3}"
DMAX="${DMAX:-14}"
TARGET="${TARGET:-43}"

while [ $# -gt 0 ]; do
    case "$1" in
        --lo)     LO="$2";     shift 2 ;;
        --hi)     HI="$2";     shift 2 ;;
        --dmin)   DMIN="$2";   shift 2 ;;
        --dmax)   DMAX="$2";   shift 2 ;;
        --target) TARGET="$2"; shift 2 ;;
        -h|--help) sed -n '2,30p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
        *) echo "[run_avxexp] unknown argument: $1" >&2; exit 1 ;;
    esac
done

if ! command -v sollya >/dev/null 2>&1; then
    echo "[run_avxexp] 'sollya' not found in PATH." >&2
    echo "[run_avxexp] Install it, then run again:" >&2
    echo "    sudo apt-get install -y sollya" >&2
    echo "    # or:  micromamba install -c conda-forge sollya" >&2
    exit 1
fi

if [ ! -f "$SOLLYA_FILE" ]; then
    echo "[run_avxexp] missing $SOLLYA_FILE" >&2
    exit 1
fi

mkdir -p "$LOG_DIR"
tag="$(printf 'exp_%s_%s_d%s-%s_t%s' "$LO" "$HI" "$DMIN" "$DMAX" "$TARGET" \
        | tr -c 'A-Za-z0-9._-' '_')"
LOG="$LOG_DIR/${tag}.log"

echo "[run_avxexp] exp(x) on [$LO, $HI], degrees $DMIN..$DMAX, target $TARGET bits"
echo "[run_avxexp] sollya: $(command -v sollya)"

# Predefine the parameters, then execute the body.  AVXPolyApprox.sollya
# fills a default only when a variable is NOT already bound (isbound), so the
# assignments below act as overrides.  Sollya's own exit code is unreliable
# (nonzero even on success), so we route through tee and inspect the log.
{
    echo "LO = $LO;"
    echo "HI = $HI;"
    echo "DMIN = $DMIN;"
    echo "DMAX = $DMAX;"
    echo "TARGET = $TARGET;"
    echo "execute(\"$SOLLYA_FILE\");"
} | sollya 2>&1 | tee "$LOG"

echo
if grep -q "CHOSEN degree" "$LOG"; then
    echo "[run_avxexp] success; log written to: $LOG"
elif grep -q "No degree" "$LOG"; then
    echo "[run_avxexp] no degree met the target (raise --dmax); log: $LOG" >&2
    exit 2
else
    echo "[run_avxexp] sollya produced no result (check the log): $LOG" >&2
    exit 3
fi
