# Performance Patch Review Checklist

Use this checklist for any optimization, vectorization, approximation, packing, parser, sampler, or generated-constant patch that can affect the signing or verification implementation. A reviewer should be able to answer each item from the pull request description, the linked reports, or a short code comment near the changed invariant.

## Required Statements

- Distribution: state whether the patch changes any sampled distribution, rejection distribution, acceptance bound, or retry distribution; if it does, link the design note and new KAT evidence.
- Randomness consumption: state whether the patch changes `randombytes()` call count, byte count, ordering, or retry-dependent consumption; if it does, update the trace KAT report.
- Parser rejection conditions: state whether the patch changes signature, key, rANS, padding, length, or malformed-input rejection conditions; if it does, add a negative parser case.
- Machine-code CT envelope: state whether the patch changes secret-dependent arithmetic, memory access, branches, division, modulo, gather/scatter, floating-point, or table lookups; if it does, link the `./ct_scan_matrix.sh` report and justify any public-data-only allowlist.
- Byte-for-byte KAT status: state whether the default build remains byte-for-byte identical to the recorded KAT hashes; if not, explain why the change is intentional and update all backend KAT records together.
- Public-data-only allowlist: state whether new divisions, modulo operations, gathers, variable shifts, early exits, or table lookups operate only on public data; cite `docs/design-notes/SECRET_PUBLIC_AUDIT.md` for the classification.

## Reviewer Gate

Before approval, run or inspect the relevant retained evidence. The default hard-gate command is `./run_tests.sh`, and the integration-audit command is `./integration_audit_matrix.sh`. Optional timing evidence for deterministic components is `./run_dudect_components.sh`; treat `WARN` as a request for a quieter rerun or a stronger static argument, not as an automatic rejection. Retained native-Ubuntu end-to-end signing dudect evidence is `test/dudect_all.txt`, which passed all backends and modes with `DUDECT_N=200000` on 2026-06-05.
