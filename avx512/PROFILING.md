# `avx512/` profiling — per-component cycle breakdown

AVX-512 build. Benchmark context: i7-11700K (Rocket Lake, 2021), Turbo Boost disabled, Hyper-Threading disabled, native Ubuntu, gcc `-O3`, `taskset`-pinned; generated from `test/profiling.txt` on 2026-06-09 via `make profile_txt` (NKG=500 keygen, NSIG=2000 sign, 2000 verify). **cyc/op = sum of the instrumented top-level buckets** (rdtscp probes, about 50 cycles each); use it for time-share attribution, not absolute wall-clock (use `make speed_txt`). Rows are sorted by time-share descending; `↳` rows break down their parent. Lithium-120 is shown in detail; per-mode totals are at the end. Default integer exp is shown first; the IFMA exp variant follows. Component design and the port-model campaign: [DESIGN.md §16](DESIGN.md).

## keygen (Lithium-120)

`sample_noise` remains the largest mode-120 keygen bucket at 136.5k cyc/op; `compute_b` is 53.9k cyc/op after the shared Montgomery batch inversion. The profile reflects Rocket Lake with Turbo Boost and Hyper-Threading disabled, so AVX-512 frequency effects are visible in end-to-end speed numbers.

| component | cyc/op | % of keygen |
|---|--:|--:|
| sample_noise (f0,s0,e) | 136,499 | 54.19% |
| compute_b (NTT/invert) | 53,941 | 21.41% |
| pack pk/sk | 55,212 | 21.92% |
| expand_A0 | 5,137 | 2.04% |
| norm window | 1,099 | 0.44% |
| **TOTAL (sum)** | **251,890** | 100% |

## sign (Lithium-120)

| component | cyc/op | % of sign |
|---|--:|--:|
| sample_y (Gaussian) | 130,157 | 55.93% |
| ↳ shake (keccak) | 45,833 | 35.21% of parent |
| ↳ ApproxExp | 53,397 | 41.03% of parent |
| ↳ a_q60+finalize+y | 9,946 | 7.64% of parent |
| ↳ BaseSampler (sigma2) | 8,884 | 6.83% of parent |
| sign_core (rejection) | 52,967 | 22.76% |
| ↳ approx_exp_wide/small | 4,019 | 7.59% of parent |
| pack (rANS) | 19,763 | 8.49% |
| expand_A0 | 3,996 | 1.72% |
| commitment (NTT mat-mul) | 8,137 | 3.50% |
| ↳ ntt forward | 1,207 | 14.84% of parent |
| ↳ ntt pointwise/basemul | 1,258 | 15.46% of parent |
| ↳ ntt inverse | 812 | 9.99% of parent |
| challenge (hash+sample) | 9,875 | 4.24% |
| setup (tr/mu/rhoprime) | 5,965 | 2.56% |
| highbits/lsb | 1,025 | 0.44% |
| normcheck | 460 | 0.20% |
| makehint | 381 | 0.16% |
| **TOTAL (sum)** | **232,731** | 100% |

## verify (Lithium-120)

The shared O(1) rANS-decode slot LUT keeps `unpack` to 23.6k cyc/op in mode 120; verify is now spread across unpack, setup/usehint/challenge, and the public-matrix NTT work.

| component | cyc/op | % of verify |
|---|--:|--:|
| unpack pk/sig | 23,575 | 37.43% |
| expand_A0 + b-hat NTT | 3,917 | 6.22% |
| commitment (NTT mat-mul) | 4,338 | 6.89% |
| setup (tr/mu) | 11,841 | 18.80% |
| usehint/lsb | 5,740 | 9.11% |
| challenge (hash+sample) | 10,041 | 15.94% |
| normcheck (z reconstruct) | 3,523 | 5.59% |
| **TOTAL (sum)** | **62,978** | 100% |

## IFMA exp variant (`make AVX512IFMAEXP=1`, opt-in)

The default build still uses the scalar integer exp (`approx_exp_4x`), so ApproxExp is the largest `sample_y` child above. The 8-lane IFMA kernel replaces it; keygen/verify are effectively unchanged apart from run noise. In mode 120 profiling, the `ApproxExp` child drops **53,397 -> 11,015 cyc/op**, `sample_y` drops **130,157 -> 87,464 cyc/op**, and the sign profiler sum drops **232,731 -> 188,721 cyc/op**. End-to-end median `speed_sign` drops **211k -> 166k (Lithium-120)**, **344k -> 284k (Lithium-180)**, and **442k -> 361k (Lithium-260)** in the 2026-06-11 `make speed_txt` run; see [`test/speed.txt`](test/speed.txt).

## Per-mode totals (cyc/op, profiler sum, default integer exp)

| primitive | Lithium-120 | Lithium-180 | Lithium-260 |
|---|--:|--:|--:|
| keygen | 251,890 | 281,682 | 435,226 |
| sign | 232,731 | 375,950 | 488,001 |
| verify | 62,978 | 95,955 | 115,138 |

vs the scalar [`ref/PROFILING.md`](../ref/PROFILING.md) and [`avx2/PROFILING.md`](../avx2/PROFILING.md).
