## Change Summary

Describe the patch and the intended performance or maintainability effect.

## Performance Patch Review

If this PR touches signing, verification, sampling, packing, parsing, generated constants, vector code, or approximation code, fill out `docs/design-notes/PERF_PATCH_REVIEW.md` and summarize the answers here.

- Distribution changed: yes/no, with evidence.
- Randomness consumption changed: yes/no, with trace evidence.
- Parser rejection conditions changed: yes/no, with negative-parser evidence.
- Machine-code CT envelope changed: yes/no, with `./ct_scan_matrix.sh` evidence.
- Default byte-for-byte KAT unchanged: yes/no, with backend KAT evidence.
- New public-data-only allowlist needed: yes/no, with classification.

## Commands Run

- `./integration_audit_matrix.sh`
- Add any optional or backend-specific commands here.
