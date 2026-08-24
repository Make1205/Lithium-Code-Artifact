# `avx2/` profiling — per-component cycle breakdown

Pure-AVX2 build (`-mno-avx512f`). Benchmark context: i7-11700K (Rocket Lake, 2021), Turbo Boost disabled, Hyper-Threading disabled, native Ubuntu, gcc `-O3`, `taskset`-pinned; generated from `test/profiling.txt` on 2026-06-09 via `make profile_txt` (NKG=500 keygen, NSIG=2000 sign, 2000 verify). **cyc/op = sum of the instrumented top-level buckets** (rdtscp probes, about 50 cycles each); use it for time-share attribution, not absolute wall-clock (use `make speed_txt`). Rows are sorted by time-share descending; `↳` rows break down their parent. Lithium-120 is shown in detail; per-mode totals are at the end. Default integer exp is shown first; the AVX2Exp float variant follows. Component design and the port-model audit: [DESIGN.md §16](DESIGN.md).

## keygen (Lithium-120)

Two keygen bottlenecks remain visible in this Rocket Lake run: `sample_noise` is still the largest bucket after the 8-way SIMD CDT vectorization, while `compute_b` is the residual NTT/inversion cost after Montgomery batch inversion. The live mode-120 totals are `sample_noise` 207.1k cyc/op and `compute_b` 49.3k cyc/op; see DESIGN §16.

| component | cyc/op | % of keygen |
|---|--:|--:|
| sample_noise (f0,s0,e) | 207,083 | 65.51% |
| compute_b (NTT/invert) | 49,267 | 15.58% |
| pack pk/sk | 48,029 | 15.19% |
| expand_A0 | 10,198 | 3.23% |
| norm window | 1,541 | 0.49% |
| **TOTAL (sum)** | **316,120** | 100% |

## sign (Lithium-120)

| component | cyc/op | % of sign |
|---|--:|--:|
| sample_y (Gaussian) | 229,198 | 64.71% |
| ↳ shake (keccak) | 120,612 | 52.62% of parent |
| ↳ ApproxExp | 56,180 | 24.51% of parent |
| ↳ a_q60+finalize+y | 26,005 | 11.35% of parent |
| ↳ BaseSampler (sigma2) | 13,832 | 6.04% of parent |
| sign_core (rejection) | 66,585 | 18.80% |
| ↳ approx_exp_wide/small | 3,891 | 5.84% of parent |
| pack (rANS) | 19,795 | 5.59% |
| expand_A0 | 9,412 | 2.66% |
| commitment (NTT mat-mul) | 10,611 | 3.00% |
| ↳ ntt forward | 2,001 | 18.87% of parent |
| ↳ ntt pointwise/basemul | 2,363 | 22.27% of parent |
| ↳ ntt inverse | 1,189 | 11.21% of parent |
| challenge (hash+sample) | 9,944 | 2.81% |
| setup (tr/mu/rhoprime) | 5,896 | 1.66% |
| highbits/lsb | 1,212 | 0.34% |
| normcheck | 858 | 0.24% |
| makehint | 699 | 0.20% |
| **TOTAL (sum)** | **354,213** | 100% |

## verify (Lithium-120)

The shared O(1) rANS-decode slot LUT keeps `unpack` to 23.6k cyc/op in mode 120; verify is now spread across unpack, `A0` expansion, commitment NTT, setup/usehint, and challenge.

| component | cyc/op | % of verify |
|---|--:|--:|
| unpack pk/sig | 23,627 | 33.30% |
| expand_A0 + b-hat NTT | 9,348 | 13.18% |
| commitment (NTT mat-mul) | 6,537 | 9.21% |
| setup (tr/mu) | 11,729 | 16.53% |
| usehint/lsb | 5,783 | 8.15% |
| challenge (hash+sample) | 10,220 | 14.40% |
| normcheck (z reconstruct) | 3,706 | 5.22% |
| **TOTAL (sum)** | **70,952** | 100% |

## AVX2Exp variant (`make AVX2EXP=1`, opt-in float exp)

The masked-OR-tree + B=8 degree-6 float kernel replaces only the Gaussian-sampler exp; keygen/verify are not designed to benefit and their measured medians vary only by run noise. In mode 120 profiling, the `ApproxExp` child drops **56,180 -> 26,808 cyc/op** and `sample_y` drops **229,198 -> 198,841 cyc/op**. End-to-end median `speed_sign` drops **337k -> 303k (Lithium-120)**, **547k -> 494k (Lithium-180)**, and **689k -> 639k (Lithium-260)** in the 2026-06-11 `make speed_txt` run; see [`test/speed.txt`](test/speed.txt).

## Per-mode totals (cyc/op, profiler sum, default integer exp)

| primitive | Lithium-120 | Lithium-180 | Lithium-260 |
|---|--:|--:|--:|
| keygen | 316,120 | 449,138 | 685,160 |
| sign | 354,213 | 567,160 | 729,147 |
| verify | 70,952 | 109,580 | 128,201 |

vs the scalar [`ref/PROFILING.md`](../ref/PROFILING.md); the AVX-512 build is [`avx512/PROFILING.md`](../avx512/PROFILING.md).
