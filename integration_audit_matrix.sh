#!/bin/sh
# Full integration audit gate: cross-backend differential fuzz plus checklist plumbing.
set -eu
cd "$(dirname "$0")"

mkdir -p test
DIFF_REPORT=${DIFF_FUZZ_REPORT:-test/diff_fuzz_matrix.txt}
CHECKLIST_REPORT=${PERF_REVIEW_REPORT:-test/perf_review_checklist_status.txt}
MATRIX_REPORT=${INTEGRATION_AUDIT_REPORT:-test/integration_audit_matrix.txt}
BACKENDS=${INTEGRATION_AUDIT_BACKENDS:-ref avx2 avx512}
MODES=${INTEGRATION_AUDIT_MODES:-120 180 260}
DIFF_FUZZ_N=${DIFF_FUZZ_N:-32}

overall=PASS
: > "$DIFF_REPORT"
{
    printf '=== Lithium-Code differential fuzz matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'DIFF_FUZZ_N: %s\n\n' "$DIFF_FUZZ_N"
} >> "$DIFF_REPORT"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT INT HUP TERM

for backend in $BACKENDS; do
    out="$tmpdir/$backend.out"
    {
        printf '## backend=%s target=diff_fuzz\n' "$backend"
        if ! make --no-print-directory -C "$backend" MODES="$MODES" TESTS=diff_fuzz DIFF_FUZZ_N="$DIFF_FUZZ_N" check; then
            printf 'BACKEND %s FAIL\n' "$backend"
            exit 1
        fi
        printf '\n'
    } > "$out" 2>&1 || overall=FAIL
    cat "$out" >> "$DIFF_REPORT"
    awk '/^(CASE|PASS diff_fuzz)/ { print }' "$out" > "$tmpdir/$backend.norm"
done

ref_backend=
for backend in $BACKENDS; do
    ref_backend=$backend
    break
done
if [ -n "$ref_backend" ]; then
    for backend in $BACKENDS; do
        if ! cmp -s "$tmpdir/$ref_backend.norm" "$tmpdir/$backend.norm"; then
            printf 'DIFF-CONSISTENCY backend=%s vs=%s FAIL\n' "$backend" "$ref_backend" >> "$DIFF_REPORT"
            diff -u "$tmpdir/$ref_backend.norm" "$tmpdir/$backend.norm" >> "$DIFF_REPORT" || true
            overall=FAIL
        else
            lines=$(wc -l < "$tmpdir/$backend.norm" | tr -d ' ')
            printf 'DIFF-CONSISTENCY backend=%s vs=%s PASS normalized_lines=%s\n' "$backend" "$ref_backend" "$lines" >> "$DIFF_REPORT"
        fi
    done
fi
printf 'OVERALL: %s\n' "$overall" >> "$DIFF_REPORT"

check_status=PASS
{
    printf '=== Lithium-Code performance-review checklist status ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    for path in docs/design-notes/PERF_PATCH_REVIEW.md .github/pull_request_template.md; do
        if [ -f "$path" ]; then
            printf 'PATH %s PASS\n' "$path"
        else
            printf 'PATH %s FAIL missing\n' "$path"
            check_status=FAIL
        fi
    done
    if grep -q 'PERF_PATCH_REVIEW.md' README.md; then
        printf 'README_LINK PERF_PATCH_REVIEW.md PASS\n'
    else
        printf 'README_LINK PERF_PATCH_REVIEW.md FAIL\n'
        check_status=FAIL
    fi
    printf 'OVERALL: %s\n' "$check_status"
} > "$CHECKLIST_REPORT"
[ "$check_status" = PASS ] || overall=FAIL

{
    printf '=== Lithium-Code integration audit matrix ===\n'
    printf 'date_utc: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
    printf 'backends: %s\n' "$BACKENDS"
    printf 'modes: %s\n' "$MODES"
    printf 'DIFF_FUZZ_N: %s\n' "$DIFF_FUZZ_N"
    printf 'diff_report: %s\n' "$DIFF_REPORT"
    printf 'checklist_report: %s\n\n' "$CHECKLIST_REPORT"
    printf '## %s\n' "$DIFF_REPORT"
    tail -n 16 "$DIFF_REPORT"
    printf '\n## %s\n' "$CHECKLIST_REPORT"
    cat "$CHECKLIST_REPORT"
    printf '\nOVERALL: %s\n' "$overall"
} > "$MATRIX_REPORT"
cat "$MATRIX_REPORT"
[ "$overall" = PASS ]
