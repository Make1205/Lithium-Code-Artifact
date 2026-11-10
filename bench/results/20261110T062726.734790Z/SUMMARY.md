# LITHIUM — QUICK

QUICK is a workflow check only. Full-scale figures are measured locally.

Calibration status: QUICK is not numerical validation. CALIBRATION=1 FULL checks the supplied historical medians within 3% on the stated CPU/GCC and stops on failure; source/message/timer differences remain disclosed.

CPU: AMD EPYC 9V74 80-Core Processor; affinity: 0. GCC: gcc (Debian 14.2.0-19) 14.2.0.

Tested commit: `81389d1971ef64f6469167ab10855e3594ea7729`. Benchmark commit: `36fe62a9f57bb76207a025eeb174470cfef6c7d3`.

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
| 120 | 256 | 57601 | 1056 | 1187; unpadded avg 1100.1600, max 1112 | 1187 | 351.441/447.407/113.021 | 132.001/198.119/50.335 | 0.480/0.599/0.664 ; 0.465/0.568/0.562 |
| 180 | 256 | 30977 | 1472 | 1618; unpadded avg 1523.6000, max 1535 | 1618 | 742.455/667.237/171.885 | 268.345/309.243/74.073 | 0.727/0.563/0.687 ; 0.694/0.594/0.612 |
| 260 | 512 | 64513 | 2080 | 2348; unpadded avg 2242.5600, max 2258 | 2348 | 767.545/828.957/228.721 | 338.623/398.397/93.157 | 0.538/0.544/0.737 ; 0.607/0.617/0.650 |



Expanded-key Sign (thousand ticks; `--` means no supported expanded-key API):

| 参数集 | ref | AVX2 |
|---|---:|---:|
| 120 | -- | -- |
| 180 | -- | -- |
| 260 | -- | -- |

Correctness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.

Calibration gate: **NOT_VALIDATED_QUICK**. Machine matches reference: False. Historical COPA messages are 59 bytes; this run follows the required 32-byte protocol.

Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider. Calibration audit: extra ref auto-vectorization suppression has been removed. Historical COPA ratios are descriptive only, not 3% acceptance evidence. See PROTOCOL.md and the detailed COPA comparison in NOTES.md. The revised timer uses fresh KeyGen/Sign/Verify triples and discards 20 warmup triples. CALIBRATION=1 sets the same default BENCH_SEED as COPA; without BENCH_SEED the native OS provider is retained.

Build flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.
