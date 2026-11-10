# LITHIUM — QUICK

QUICK is a workflow check only. Full-scale figures are measured locally.

Calibration status: Report-only historical COPA ratios; no 3% gate, machine mismatch or ratio can fail the run.

CPU: AMD EPYC 9V74 80-Core Processor; affinity: 0. GCC: gcc (Debian 14.2.0-19) 14.2.0.

Tested commit: `81389d1971ef64f6469167ab10855e3594ea7729`. Benchmark commit: `c9862dd7e387fc5860f10aa661cd97586ad91ed7`.

Base flags: `-O3 -march=native -std=c99 -ffp-contract=off`. Each operation: 11 iterations; correctness and lengths: 100 signatures per build/set. Message: 32 bytes. Median units: thousand RDTSC ticks.

Exact per-build options (in addition to base flags and local include paths):

* 120 ref: `-DLITHIUM_MODE=120 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 120 avx2: `-DLITHIUM_MODE=120 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 120 vendored_ref: `-DLITHIUM_MODE=120 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 120 vendored_avx2: `-DLITHIUM_MODE=120 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 120 vendored_ref_instrumented: `-DLITHIUM_MODE=120 -Itest -I../ref -I../ref/test -include config.h -DBENCH_COUNT_ATTEMPTS`; link `-lm`.
* 120 vendored_avx2_instrumented: `-DLITHIUM_MODE=120 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h -DBENCH_COUNT_ATTEMPTS`; link `-lm`.
* 180 ref: `-DLITHIUM_MODE=180 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 avx2: `-DLITHIUM_MODE=180 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 vendored_ref: `-DLITHIUM_MODE=180 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 vendored_avx2: `-DLITHIUM_MODE=180 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 vendored_ref_instrumented: `-DLITHIUM_MODE=180 -Itest -I../ref -I../ref/test -include config.h -DBENCH_COUNT_ATTEMPTS`; link `-lm`.
* 180 vendored_avx2_instrumented: `-DLITHIUM_MODE=180 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h -DBENCH_COUNT_ATTEMPTS`; link `-lm`.
* 260 ref: `-DLITHIUM_MODE=260 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 avx2: `-DLITHIUM_MODE=260 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 vendored_ref: `-DLITHIUM_MODE=260 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 vendored_avx2: `-DLITHIUM_MODE=260 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 vendored_ref_instrumented: `-DLITHIUM_MODE=260 -Itest -I../ref -I../ref/test -include config.h -DBENCH_COUNT_ATTEMPTS`; link `-lm`.
* 260 vendored_avx2_instrumented: `-DLITHIUM_MODE=260 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h -DBENCH_COUNT_ATTEMPTS`; link `-lm`.

| 参数集 | n | q | pk | sig（仓库实际） | 论文 sig | ref KeyGen/Sign/Verify | AVX2 KeyGen/Sign/Verify | ratio to previous ref / AVX2 |
|---|---:|---:|---:|---|---:|---|---|---|
| 120 | 256 | 57601 | 1056 | 1187; unpadded avg 1100.1600, max 1112 | 1187 | 352.091/446.835/111.851 | 132.105/199.133/50.569 | 0.481/0.598/0.657 ; 0.466/0.571/0.565 |
| 180 | 256 | 30977 | 1472 | 1618; unpadded avg 1523.6000, max 1535 | 1618 | 732.679/667.081/174.303 | 268.761/306.747/74.203 | 0.717/0.563/0.696 ; 0.695/0.590/0.613 |
| 260 | 512 | 64513 | 2080 | 2348; unpadded avg 2242.5600, max 2258 | 2348 | 759.693/835.821/226.381 | 278.849/399.931/93.521 | 0.533/0.549/0.730 ; 0.500/0.619/0.652 |



Expanded-key Sign (thousand ticks; `--` means no supported expanded-key API):

| 参数集 | ref | AVX2 |
|---|---:|---:|
| 120 | -- | -- |
| 180 | -- | -- |
| 260 | -- | -- |

Correctness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.

Historical COPA ratios are report-only; no 3% gate or machine check can reject this run. Historical messages: 59 bytes; current protocol: 32 bytes.


Source attribution with the SAME new driver (thousand cycles):

| Set | Build | Upstream 81389d1 | COPA vendor, instrumentation off | COPA vendor, native instrumentation on | Vendor/upstream | Instrumented/vendor |
|---|---|---|---|---|---|---|
| 120 | ref | 352.091/446.835/111.851 | 350.817/441.739/112.007 | 453.101/442.909/113.385 | 0.996/0.989/1.001 | 1.292/1.003/1.012 |
| 120 | avx2 | 132.105/199.133/50.569 | 134.627/198.223/50.569 | 132.885/196.507/50.439 | 1.019/0.995/1.000 | 0.987/0.991/0.997 |
| 180 | ref | 732.679/667.081/174.303 | 745.029/676.259/173.913 | 743.989/731.171/173.237 | 1.017/1.014/0.998 | 0.999/1.081/0.996 |
| 180 | avx2 | 268.761/306.747/74.203 | 269.151/307.163/75.165 | 271.361/305.759/75.191 | 1.001/1.001/1.013 | 1.008/0.995/1.000 |
| 260 | ref | 759.693/835.821/226.381 | 803.971/862.289/227.395 | 726.985/836.887/228.513 | 1.058/1.032/1.004 | 0.904/0.971/1.005 |
| 260 | avx2 | 278.849/399.931/93.521 | 320.605/398.579/93.234 | 279.083/399.021/94.015 | 1.150/0.997/0.997 | 0.870/1.001/1.008 |

All three source conditions use the identical driver, message, RNG provider, warmup, counter and flags (only the documented native instrumentation macro differs). This controls driver differences; FULL local results are needed to attribute the user-machine discrepancy quantitatively.

FMA audit: selected integer/fixed-point paths; zero FMA instructions in all benchmark-owned objects and linked drivers. See *-fma-audit.json.
Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider. Calibration audit: extra ref auto-vectorization suppression has been removed. Historical COPA ratios are descriptive only, not 3% acceptance evidence. See PROTOCOL.md and the detailed COPA comparison in NOTES.md. The revised timer uses fresh KeyGen/Sign/Verify triples and discards 20 warmup triples. CALIBRATION=1 sets the same default BENCH_SEED as COPA; without BENCH_SEED the native OS provider is retained. Current source attribution compares original 81389d1, uninstrumented COPA-vendored 3c676bd, and the same vendored source with native BENCH_COUNT_ATTEMPTS, all through one driver. CALIBRATION is a backwards-compatible seeded-RNG/report-only option, never an acceptance gate.

Build flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.
