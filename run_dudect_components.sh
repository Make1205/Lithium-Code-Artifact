#!/bin/sh
# Optional component-level dudect smoke. Records PASS/WARN/FAIL but does not hard-gate run_tests.sh.
set -eu
cd "$(dirname "$0")"

mkdir -p test
REPORT=${DUDECT_COMPONENTS_REPORT:-test/dudect_components_matrix.txt}
BACKENDS=${DUDECT_COMPONENTS_BACKENDS:-ref avx2 avx512}
MODES=${DUDECT_COMPONENTS_MODES:-120 180 260}
DUDECT_COMPONENT_N=${DUDECT_COMPONENT_N:-512}

status=PASS
{
    printf '=== Lithium-Code optional dudect components matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'DUDECT_COMPONENT_N: %s\n\n' "$DUDECT_COMPONENT_N"
} > "$REPORT"

for backend in $BACKENDS; do
    {
        printf '## backend=%s target=dudect_components\n' "$backend"
        if ! make --no-print-directory -C "$backend" MODES="$MODES" TESTS=dudect_components DUDECT_COMPONENT_N="$DUDECT_COMPONENT_N" check; then
            printf 'BACKEND %s FAIL\n' "$backend"
            exit 1
        fi
        printf '\n'
    } >> "$REPORT" 2>&1 || status=FAIL
    if rg -q 'status=FAIL|^WARN dudect_components' "$REPORT"; then
        [ "$status" = FAIL ] || status=WARN
    fi
done
printf 'OVERALL: %s\n' "$status" >> "$REPORT"
cat "$REPORT"
[ "$status" != FAIL ]
