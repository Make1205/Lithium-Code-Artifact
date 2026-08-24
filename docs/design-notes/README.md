# Design Notes

This directory contains project-local design references used by the generators, reproducibility checks, security review, and backend implementation notes. The notes are committed here so reviewers can audit the scheme from the repository alone, without relying on private planning files.

## Reading Guide

- `param.md`, `Trunc.md`, `IRS.md`, `IRS_NIST_rANS.md`, and `RCDT.md` explain the scheme parameters, rejection-sampling analysis, rANS signature layer, and Gaussian tables.
- `PreNTT.md`, `AVXExp.md`, `IRS_x64_ApproxExp.md`, `PK_UNPACK.md`, and `IRS_NIST_rANS.md` document implementation choices that affect generated constants or backend equivalence.
- `SECRET_PUBLIC_AUDIT.md`, `TIMECOP.md`, and `PERF_PATCH_REVIEW.md` are the review-facing security notes. New public-data-only allowlists, parser early exits, or performance-sensitive changes should cite these notes.

## Source Of Truth

Executable sources of truth for generated constants live under `tools/`, especially `tools/ntt/`, `tools/avxexp/`, `gen_zetas.py`, `gen_bounds.py`, `gen_rcdt.py`, `gen_rans_tables.py`, and `SigSize.py`. Earlier private derivations are intentionally summarized in public source comments and design notes rather than copied verbatim.
