# LITHIUM — QUICK

QUICK is a workflow check only. Full-scale figures are measured locally.

Calibration status: NOT VALIDATED: COPA 3c676bd is a different 59-byte/deterministic/fresh-key/instrumented workload and contains a changed AVX2 SHAKE source. QUICK cannot establish a 3% local FULL reproduction.

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
| 120 | 256 | 57601 | 1056 | 1187; unpadded avg 1100.1600, max 1112 | 1187 | 397.929/444.027/84.317 | 192.347/199.705/41.859 | 0.544/0.594/0.496 ; 0.677/0.572/0.470 |
| 180 | 256 | 30977 | 1472 | 1618; unpadded avg 1523.6000, max 1535 | 1618 | 744.457/675.895/134.055 | 300.819/309.555/64.271 | 0.729/0.570/0.536 ; 0.777/0.595/0.531 |
| 260 | 512 | 64513 | 2080 | 2348; unpadded avg 2242.5600, max 2258 | 2348 | 1062.749/830.413/162.265 | 395.303/397.695/78.155 | 0.745/0.545/0.523 ; 0.710/0.616/0.547 |



Expanded-key Sign (thousand ticks; `--` means no supported expanded-key API):

| 参数集 | ref | AVX2 |
|---|---:|---:|
| 120 | -- | -- |
| 180 | -- | -- |
| 260 | -- | -- |

Correctness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.

Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider. Calibration audit: extra ref auto-vectorization suppression has been removed. Historical COPA ratios are descriptive only, not 3% acceptance evidence. See PROTOCOL.md and the detailed COPA comparison in NOTES.md.

Build flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.
