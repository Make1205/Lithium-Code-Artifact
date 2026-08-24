#!/bin/sh
# End-to-end Lithium-Code gate.  This script is intentionally a thin fail-fast
# wrapper around the Makefile targets, so the source lists have one owner.
set -eu

cd "$(dirname "$0")"

run()
{
    printf '\n== %s ==\n' "$*"
    "$@"
}

run make -C ref tables
run make -C ref check
run make -C ref check-rans-failrate
run ./ct_scan_matrix.sh
run ./security_audit_matrix.sh
run ./integration_audit_matrix.sh

run make -C avx2 check-consts
run make -C avx2 check
run make -C avx2 check-kat
run make -C avx2 check-kat-avx2exp
run make -C avx2 check-rans-failrate
run make -C avx2 ct-scan
run make -C avx2 check-no-avx512
run make -C avx2 check-refill

run make -C avx512 check-consts
run make -C avx512 AVX512IFMAEXP=1 check-consts
run make -C avx512 check
run make -C avx512 check-kat
run make -C avx512 check-kat-avx512ifma
run make -C avx512 check-rans-failrate
run make -C avx512 ct-scan
run make -C avx512 check-refill

# Optional KyberSlash/TIMECOP dynamic cross-check. Off by default because it
# builds a patched Valgrind on first use (download + ~3 min compile); the
# variable-latency-division property is already gated statically by ct-scan.
if [ "${TIMECOP:-0}" = 1 ]; then
    run sh tools/timecop/run_all.sh
fi

printf '\nrun_tests.sh: all gates passed\n'
