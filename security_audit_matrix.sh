#!/bin/sh
# Full security audit gate: fault injection, negative parser corpus, and trace KAT.
set -eu
cd "$(dirname "$0")"

mkdir -p test
FAULT_REPORT=test/fault_injection_matrix.txt
PARSER_REPORT=test/parser_negative_matrix.txt
TRACE_REPORT=test/trace_kat_matrix.txt
MATRIX_REPORT=${SECURITY_AUDIT_REPORT:-test/security_audit_matrix.txt}
BACKENDS=${SECURITY_AUDIT_BACKENDS:-ref avx2 avx512}
MODES=${SECURITY_AUDIT_MODES:-120 180 260}

run_group()
{
    label=$1
    target=$2
    report=$3
    tmp=$(mktemp)
    : > "$tmp"
    status=PASS
    {
        printf '=== Lithium-Code security %s ===\n' "$label"
        printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
        printf 'backends: %s\n' "$BACKENDS"
        printf 'modes: %s\n\n' "$MODES"
    } > "$report"
    for backend in $BACKENDS; do
        {
            printf '## backend=%s target=%s\n' "$backend" "$target"
            make --no-print-directory -C "$backend" MODES="$MODES" TESTS="$target" check
            printf '\n'
        } >> "$report" 2>&1 || status=FAIL
    done
    if [ "$target" = trace_kat ] && [ "$status" = PASS ]; then
        awk '/PASS trace_kat/ { for (i = 1; i <= NF; i++) { if ($i ~ /^mode=/) { split($i, m, "="); mode=m[2] } if ($i ~ /^trace_hash=/) { split($i, h, "="); print mode, h[2] } } }' "$report" | sort > "$tmp"
        for mode in $MODES; do
            count=$(awk -v m="$mode" '$1 == m { print $2 }' "$tmp" | sort -u | wc -l | tr -d ' ')
            if [ "$count" != 1 ]; then
                printf 'TRACE-CONSISTENCY mode=%s FAIL unique_hashes=%s\n' "$mode" "$count" >> "$report"
                status=FAIL
            else
                hash=$(awk -v m="$mode" '$1 == m { print $2; exit }' "$tmp")
                printf 'TRACE-CONSISTENCY mode=%s PASS trace_hash=%s\n' "$mode" "$hash" >> "$report"
            fi
        done
    fi
    rm -f "$tmp"
    printf 'OVERALL: %s\n' "$status" >> "$report"
    [ "$status" = PASS ]
}

overall=PASS
run_group 'fault injection matrix' security_fault_injection "$FAULT_REPORT" || overall=FAIL
run_group 'parser negative matrix' parser_negative "$PARSER_REPORT" || overall=FAIL
run_group 'trace KAT matrix' trace_kat "$TRACE_REPORT" || overall=FAIL
{
    printf '=== Lithium-Code security audit matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'fault_report: %s\n' "$FAULT_REPORT"
    printf 'parser_report: %s\n' "$PARSER_REPORT"
    printf 'trace_report: %s\n\n' "$TRACE_REPORT"
    for report in "$FAULT_REPORT" "$PARSER_REPORT" "$TRACE_REPORT"; do
        printf '## %s\n' "$report"
        tail -n 12 "$report"
        printf '\n'
    done
    printf 'OVERALL: %s\n' "$overall"
} > "$MATRIX_REPORT"
cat "$MATRIX_REPORT"
[ "$overall" = PASS ]
