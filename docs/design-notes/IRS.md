# IRS Algorithm Notes

This file summarizes the compressed IRS signing API implemented by this repository. The public API follows the standard keygen, detached-signature, verification, and attached-signature wrappers in `ref/sign.c` and `ref/api.h`.

The scheme stores the NTWE public key as `b = (A0*s0 + e)*f^{-1} mod q`. Verification reconstructs the implicit matrix `A = [qj - 2b, 2A0, 2I] mod 2q` from the public key and the mod-2q lift. The implementation keeps `A0` expansion deterministic from the seed and uses the compact rANS signature encoding documented in `docs/design-notes/IRS_NIST_rANS.md`.

The algorithmic details that are executable in this repository are encoded in the C tests, Python reference checks, KAT generator, and reproducible-constant generators under `tools/`.
