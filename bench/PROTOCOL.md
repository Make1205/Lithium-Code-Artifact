# Shared byte-key benchmark rule (2026-11-10 revision)

This is the common benchmark rule. The user has authorized its rollout and
the FP-contraction amendment below. Historical COPA comparisons are now
report-only: there is no 3% acceptance gate, in QUICK or FULL.

1. Use system GCC with -O3 -march=native plus options actually needed by the
   repository backend, its headers/ABI and build. For LITHIUM use -std=c99,
   the namespace configuration, and its native AVX2/NTT/SHAKE selector macros.
2. A portable reference source does not by itself require disabling compiler
   auto-vectorization. Do not add -fno-tree-vectorize or
   -fno-tree-slp-vectorize unless the repository's own ref build requires them.
   Hand-written SIMD backend selection remains distinct from auto-vectorization.
3. Sign starts with byte-form sk and includes its decoding/expansion. Verify
   starts with byte-form pk and includes its decoding/expansion. Required message
   hashing, signature encoding/decoding and native API work stay inside the call.
4. Message length is 32 bytes. Expanded-key operations, if available, are separate
   supplemental rows and must not replace the byte-key rows.
5. The four Dilithium bba9534 timing files remain byte-identical. Record exact
   timing boundaries, counter instruction, overhead subtraction and sample count.
   Different timer implementations must not be described as identical.
6. Record native versus deterministic entropy, seed/corpus provenance, key reuse
   versus fresh-key sampling, warmups, CPU affinity and compiler/source identities.
   A single-key Sign median and a fresh-key-per-sample median are different
   statistics even if both start with serialized keys.
7. QUICK=1 uses 11 timing samples and 100 correctness/length signatures; only this
   mode runs in cloud. FULL uses 10001 timing samples and 10000 signatures and is
   executed locally. Upstream KAT case counts remain unchanged.
8. Historical COPA comparisons are report-only. No percentage, CPU/GCC mismatch
   or ratio may fail a run. CALIBRATION remains a backwards-compatible report
   option; controlled attribution defaults to the recorded BENCH_SEED. FULL
   same-machine runs separate source and driver/protocol factors; QUICK is only
   a smoke check. Build failures, KAT, correctness and FMA audit still fail.

Historical protocol differences (not a blocker): COPA 3c676bd does not use this workload verbatim: it uses a
59-byte message, deterministic default entropy, fresh key triples, 20 warmups,
serialized RDTSCP, benchmark source instrumentation and a modified AVX2 SHAKE.
See NOTES.md. The first three requested rules are now applied to LITHIUM; other
sampling/timer differences remain explicit pending a comparable calibration target.


## Floating-point contraction rule (2026-11-10)

All builds that execute floating-point arithmetic, ref and AVX2 alike, must add
`-ffp-contract=off` to system GCC `-O3 -march=native` and required backend options.
This prevents compiler-generated fused multiply/add and preserves separate
binary64 rounding. It does not disable or rewrite explicit upstream FMA
intrinsics. Those are retained, named and counted separately; never remove
algorithm intrinsics just to obtain a zero count. FALCON_FMA is not changed by
this compiler policy. Integer-only paths do not need a contraction flag.

Every run checks benchmark-owned object files and the linked driver with
objdump (apt package binutils), saves per-build fma-audit JSON and reports counts
in SUMMARY. Expected count is zero unless a named original source object uses
explicit FMA intrinsics. Shared-library internals are outside this audit. Any
unexplained FMA causes the run to fail before timing. HAWK KeyGen/Sign are
integer/fixed-point, and LITHIUM's selected default exp path is integer; neither
is silently changed into another algorithm variant by this policy.

The local c-fn-dsa case 3462 (zero-based index 3461) established FP contraction
as the reason for the old mismatch, not an implementation defect: both signatures
cross-verified; disabling contraction removed 15 compiler-generated FMA
instructions and restored equality with ref for that case. c-fn-dsa now requires
ref/AVX2 byte equality and also retains all cross-verification checks. c-fn-dsa
is the paper's primary Falcon-family measurement; VecFalcon is the backup.
