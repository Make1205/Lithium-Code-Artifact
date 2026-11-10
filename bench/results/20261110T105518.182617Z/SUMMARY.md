# LITHIUM — QUICK

QUICK is a workflow check only. Full-scale figures are measured locally.

Calibration status: QUICK is not numerical validation. CALIBRATION=1 FULL checks the supplied historical medians within 3% on the stated CPU/GCC and stops on failure; source/message/timer differences remain disclosed.

CPU: AMD EPYC 9V74 80-Core Processor; affinity: 0. GCC: gcc (Debian 14.2.0-19) 14.2.0.

Tested commit: `81389d1971ef64f6469167ab10855e3594ea7729`. Benchmark commit: `c334eac892e78b9fa5f45c1bcc355dd1b55ab3a9`.

Base flags: `-O3 -march=native -std=c99`. Each operation: 11 iterations; correctness and lengths: 100 signatures per build/set. Message: 32 bytes. Median units: thousand RDTSC ticks.

Exact per-build options (in addition to base flags and local include paths):

* 120 ref: `-DLITHIUM_MODE=120 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 120 avx2: `-DLITHIUM_MODE=120 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 ref: `-DLITHIUM_MODE=180 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 avx2: `-DLITHIUM_MODE=180 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 ref: `-DLITHIUM_MODE=260 -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 avx2: `-DLITHIUM_MODE=260 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.

| 参数集 | n | q | pk | sig（仓库实际） | 论文 sig | ref KeyGen/Sign/Verify | AVX2 KeyGen/Sign/Verify | ratio to previous ref / AVX2 |
|---|---:|---:|---:|---|---:|---|---|---|
| 120 | 256 | 57601 | 1056 | 1187; unpadded avg 1100.1600, max 1112 | 1187 | 358.123/447.563/112.631 | 132.625/199.731/50.309 | 0.489/0.599/0.662 ; 0.468/0.573/0.562 |
| 180 | 256 | 30977 | 1472 | 1618; unpadded avg 1523.6000, max 1535 | 1618 | 856.075/717.833/174.069 | 268.943/309.503/76.595 | 0.838/0.606/0.696 ; 0.695/0.595/0.633 |
| 260 | 512 | 64513 | 2080 | 2348; unpadded avg 2242.5600, max 2258 | 2348 | 820.455/868.321/228.019 | 277.887/399.515/94.925 | 0.575/0.570/0.735 ; 0.499/0.619/0.662 |



Expanded-key Sign (thousand ticks; `--` means no supported expanded-key API):

| 参数集 | ref | AVX2 |
|---|---:|---:|
| 120 | -- | -- |
| 180 | -- | -- |
| 260 | -- | -- |

Correctness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.

Calibration gate: **NOT_VALIDATED_QUICK**. Machine matches reference: False. Historical COPA messages are 59 bytes; this run follows the required 32-byte protocol.


FMA audit: selected integer/fixed-point paths; zero FMA instructions in all benchmark-owned objects and linked drivers. See *-fma-audit.json.
Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider. Calibration audit: extra ref auto-vectorization suppression has been removed. Historical COPA ratios are descriptive only, not 3% acceptance evidence. See PROTOCOL.md and the detailed COPA comparison in NOTES.md. The revised timer uses fresh KeyGen/Sign/Verify triples and discards 20 warmup triples. CALIBRATION=1 sets the same default BENCH_SEED as COPA; without BENCH_SEED the native OS provider is retained.

Build flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.
