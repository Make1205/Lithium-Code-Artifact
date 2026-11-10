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
| 120 | 256 | 57601 | 1056 | 1187; unpadded avg 1100.1600, max 1112 | 1187 | 355.965/448.265/113.619 | 133.457/201.031/50.517 | 0.486/0.600/0.667 ; 0.471/0.576/0.564 |
| 180 | 256 | 30977 | 1472 | 1618; unpadded avg 1523.6000, max 1535 | 1618 | 741.519/677.923/172.535 | 269.021/382.225/76.179 | 0.726/0.572/0.689 ; 0.696/0.735/0.630 |
| 260 | 512 | 64513 | 2080 | 2348; unpadded avg 2242.5600, max 2258 | 2348 | 827.397/927.627/227.759 | 277.809/402.141/95.393 | 0.580/0.609/0.734 ; 0.498/0.623/0.666 |



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
| 120 | ref | 355.965/448.265/113.619 | 354.665/511.419/114.035 | 351.181/447.043/111.876 | 0.996/1.141/1.004 | 0.990/0.874/0.981 |
| 120 | avx2 | 133.457/201.031/50.517 | 132.001/197.651/50.777 | 132.625/197.859/50.127 | 0.989/0.983/1.005 | 1.005/1.001/0.987 |
| 180 | ref | 741.519/677.923/172.535 | 736.579/719.861/173.601 | 745.575/666.379/172.119 | 0.993/1.062/1.006 | 1.012/0.926/0.991 |
| 180 | avx2 | 269.021/382.225/76.179 | 268.605/304.589/73.761 | 271.309/307.189/73.865 | 0.998/0.797/0.968 | 1.010/1.009/1.001 |
| 260 | ref | 827.397/927.627/227.759 | 763.749/828.749/228.356 | 724.853/994.837/227.446 | 0.923/0.893/1.003 | 0.949/1.200/0.996 |
| 260 | avx2 | 277.809/402.141/95.393 | 279.967/400.919/94.170 | 278.537/397.617/94.977 | 1.008/0.997/0.987 | 0.995/0.992/1.009 |

All three source conditions use the identical driver, message, RNG provider, warmup, counter and flags (only the documented native instrumentation macro differs). This controls driver differences; FULL local results are needed to attribute the user-machine discrepancy quantitatively.

Driver/protocol residual against historical COPA (ratios; valid for attribution only on the same local environment):

| Set | Build | Vendored instrumented / upstream (code) | Vendored new-driver / historical vendor old-driver | Upstream new-driver / historical vendor old-driver |
|---|---|---|---|---|
| 120 | ref | 0.987/0.997/0.985 | 0.480/0.599/0.657 | 0.486/0.600/0.667 |
| 120 | avx2 | 0.994/0.984/0.992 | 0.468/0.567/0.560 | 0.471/0.576/0.564 |
| 180 | ref | 1.005/0.983/0.998 | 0.730/0.562/0.688 | 0.726/0.572/0.689 |
| 180 | avx2 | 1.009/0.804/0.970 | 0.701/0.590/0.611 | 0.696/0.735/0.630 |
| 260 | ref | 0.876/1.072/0.999 | 0.508/0.653/0.733 | 0.580/0.609/0.734 |
| 260 | avx2 | 1.003/0.989/0.996 | 0.500/0.616/0.663 | 0.498/0.623/0.666 |

Decomposition: overall = driver/protocol factor divided by code factor. Historical messages were 59 bytes; the message change is part of the driver/protocol factor. Comparable FULL local environment: False. This label is informational and never a pass/fail gate.

FMA audit: selected integer/fixed-point paths; zero FMA instructions in all benchmark-owned objects and linked drivers. See *-fma-audit.json.
Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider. Calibration audit: extra ref auto-vectorization suppression has been removed. Historical COPA ratios are descriptive only, not 3% acceptance evidence. See PROTOCOL.md and the detailed COPA comparison in NOTES.md. The revised timer uses fresh KeyGen/Sign/Verify triples and discards 20 warmup triples. CALIBRATION=1 sets the same default BENCH_SEED as COPA; without BENCH_SEED the native OS provider is retained. Current source attribution compares original 81389d1, uninstrumented COPA-vendored 3c676bd, and the same vendored source with native BENCH_COUNT_ATTEMPTS, all through one driver. CALIBRATION is a backwards-compatible seeded-RNG/report-only option, never an acceptance gate.

Build flags and commands per variant are in config.json, metadata.json and build logs. Algorithm sources use the pinned tested commit and, when configured, the verified vendored snapshot. Algorithm source text is unchanged. See bench/NOTES.md.
