# Shared byte-key benchmark rule (2026-11-10 revision)

This is the requested common rule for every scheme bench. Rollout beyond LITHIUM
is paused at the user's calibration gate until the reference workload is resolved.
The existing results from other repositories are not relabeled as following it.

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
8. A 3% same-machine FULL calibration requires the actual reference sample workload
   and raw medians. QUICK, a different CPU, a changed source tree or a different
   message/RNG/key corpus cannot certify that acceptance condition. Do not tune
   the implementation or choose a favorable key just to force a historical ratio.

Current blocker: COPA 3c676bd does not use this workload verbatim: it uses a
59-byte message, deterministic default entropy, fresh key triples, 20 warmups,
serialized RDTSCP, benchmark source instrumentation and a modified AVX2 SHAKE.
See NOTES.md. The first three requested rules are now applied to LITHIUM; other
sampling/timer differences remain explicit pending a comparable calibration target.
