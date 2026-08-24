# `ref/` profiling — per-component cycle breakdown

Pure-scalar reference. Benchmark context: i7-11700K (Rocket Lake, 2021), Turbo Boost disabled, Hyper-Threading disabled, native Ubuntu, gcc `-O3`, `taskset`-pinned; generated from `test/profiling.txt` on 2026-06-09 via `make profile_txt` (NKG=500 keygen, NSIG=2000 sign, 2000 verify). **cyc/op = sum of the instrumented top-level buckets** (rdtscp probes, about 50 cycles each); use it for time-share attribution, not absolute wall-clock (use `make speed_txt`). Rows are sorted by time-share descending; `↳` rows break down their parent. Lithium-120 is shown in detail; per-mode totals are at the end. The `avx2/` and `avx512/` builds vectorize these same components; see their `PROFILING.md`.

## keygen (Lithium-120)

| component | cyc/op | % of keygen |
|---|--:|--:|
| sample_noise (f0,s0,e) | 636,014 | 60.54% |
| compute_b (NTT/invert) | 277,655 | 26.43% |
| pack pk/sk | 79,628 | 7.58% |
| expand_A0 | 31,983 | 3.04% |
| norm window | 25,335 | 2.41% |
| **TOTAL (sum)** | **1,050,617** | 100% |

## sign (Lithium-120)

| component | cyc/op | % of sign |
|---|--:|--:|
| sample_y (Gaussian) | 571,022 | 57.30% |
| ↳ shake (keccak) | 367,731 | 64.40% of parent |
| ↳ ApproxExp | 67,592 | 11.84% of parent |
| ↳ a_q60+finalize+y | 42,865 | 7.51% of parent |
| ↳ BaseSampler (sigma2) | 77,577 | 13.59% of parent |
| sign_core (rejection) | 190,967 | 19.16% |
| ↳ approx_exp_wide/small | 3,782 | 1.98% of parent |
| pack (rANS) | 20,338 | 2.04% |
| expand_A0 | 28,020 | 2.81% |
| commitment (NTT mat-mul) | 143,801 | 14.43% |
| ↳ ntt forward | 84,463 | 58.74% of parent |
| ↳ ntt pointwise/basemul | 29,935 | 20.82% of parent |
| ↳ ntt inverse | 23,370 | 16.25% of parent |
| challenge (hash+sample) | 11,915 | 1.20% |
| setup (tr/mu/rhoprime) | 7,051 | 0.71% |
| highbits/lsb | 3,098 | 0.31% |
| normcheck | 12,853 | 1.29% |
| makehint | 7,420 | 0.74% |
| **TOTAL (sum)** | **996,489** | 100% |

## verify (Lithium-120)

| component | cyc/op | % of verify |
|---|--:|--:|
| unpack pk/sig | 23,600 | 13.41% |
| expand_A0 + b-hat NTT | 27,603 | 15.69% |
| commitment (NTT mat-mul) | 82,096 | 46.65% |
| setup (tr/mu) | 14,035 | 7.98% |
| usehint/lsb | 5,723 | 3.25% |
| challenge (hash+sample) | 11,925 | 6.78% |
| normcheck (z reconstruct) | 10,980 | 6.24% |
| **TOTAL (sum)** | **175,964** | 100% |

## Per-mode totals (cyc/op, profiler sum, default integer exp)

| primitive | Lithium-120 | Lithium-180 | Lithium-260 |
|---|--:|--:|--:|
| keygen | 1,050,617 | 1,493,461 | 2,168,163 |
| sign | 996,489 | 1,618,278 | 2,080,676 |
| verify | 175,964 | 278,389 | 326,588 |

See the top-level [README](../README.md#performance) for the current performance summary and retained verification snapshot.
