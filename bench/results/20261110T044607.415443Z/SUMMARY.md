# LITHIUM — QUICK

QUICK is a workflow check only. Full-scale figures are measured locally.

CPU: AMD EPYC 9V45 96-Core Processor; affinity: 0. GCC: gcc (Debian 14.2.0-19) 14.2.0.

Tested commit: `81389d1971ef64f6469167ab10855e3594ea7729`. Benchmark commit: `81389d1971ef64f6469167ab10855e3594ea7729`.

Base flags: `-O3 -march=native -std=gnu11`. Each operation: 11 iterations; correctness and lengths: 100 signatures per build/set. Message: 32 bytes. Median units: thousand RDTSC ticks.

Exact per-build options (in addition to base flags and local include paths):

* 120 ref: `-DLITHIUM_MODE=120 -fno-tree-vectorize -fno-tree-slp-vectorize -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 120 avx2: `-DLITHIUM_MODE=120 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 ref: `-DLITHIUM_MODE=180 -fno-tree-vectorize -fno-tree-slp-vectorize -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 180 avx2: `-DLITHIUM_MODE=180 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 ref: `-DLITHIUM_MODE=260 -fno-tree-vectorize -fno-tree-slp-vectorize -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.
* 260 avx2: `-DLITHIUM_MODE=260 -mavx2 -mbmi2 -mpopcnt -mno-avx512f -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X -Itest -I../ref -I../ref/test -include config.h`; link `-lm`.

| 参数集 | n | q | pk | sig（仓库实际） | 论文 sig | ref KeyGen/Sign/Verify | AVX2 KeyGen/Sign/Verify | ratio to previous ref / AVX2 |
|---|---:|---:|---:|---|---:|---|---|---|
| 120 | 256 | 57601 | 1056 | 1187; unpadded avg 1100.1600, max 1112 | 1187 | 610.271/448.213/80.261 | 105.897/158.469/32.187 | 0.834/0.600/0.472 ; 0.373/0.454/0.362 |
| 180 | 256 | 30977 | 1472 | 1618; unpadded avg 1523.6000, max 1535 | 1618 | 683.747/752.257/126.255 | 147.471/247.103/48.229 | 0.670/0.635/0.505 ; 0.381/0.475/0.399 |
| 260 | 512 | 64513 | 2080 | 2348; unpadded avg 2242.5600, max 2258 | 2348 | 773.083/1045.381/155.505 | 222.221/314.209/58.473 | 0.542/0.686/0.502 ; 0.399/0.486/0.409 |



Expanded-key Sign (thousand ticks; `--` means no supported expanded-key API):

| 参数集 | ref | AVX2 |
|---|---:|---:|
| 120 | -- | -- |
| 180 | -- | -- |
| 260 | -- | -- |

Correctness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.

Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider.

Build flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.
