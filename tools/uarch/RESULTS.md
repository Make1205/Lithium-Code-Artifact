# AVX-512-capable benchmark system microarchitecture reference

Core = i7-11700K (Rocket Lake, 2021; Ice-Lake-family client core backported to 14nm). Measured on-device via `tools/uarch/ubench.c` (rdtsc + `vpaddq`-latency anchor = 1 core cycle; the test environment had no usable PMU). The current benchmark corpus was recorded on native Ubuntu with Turbo Boost disabled and Hyper-Threading disabled, so the cycle counts in `test/speed.txt` and `test/profiling.txt` use a stable fixed-frequency setup.

## Measured table (512-bit zmm unless noted)

| instruction            | latency (cyc) | throughput (ops/cyc) | port(s) |
|------------------------|:-------------:|:--------------------:|---------|
| `vpaddq`/`vpsubd`/`vpsubq` (int add/sub) | 1 | **2.0** | p0+p5 |
| `vpternlogq`/`vpandq`/`vpxorq` (logic)   | 1 | **2.0** | p0+p5 |
| `vpmadd52luq` (IFMA lo) | **4** | **1.0** | **p0** |
| `vpmadd52huq` (IFMA hi) | **4** | **1.0** | **p0** |
| `vpmullw` (16b mul lo)  | 5 | 1.0 | **p0** |
| `vpmulhw` (16b mul hi)  | 5 | 1.0 | **p0** |
| `vpsllq`/`vpsrlq` (imm shift)   | 1 | 1.0 | **p0** |
| `vpsllvq`/`vpsrlvq` (var shift) | 1 | 1.0 | **p0** |
| `vprolq` (rotate)       | 1 | 1.0 | **p0** |
| `vpermt2q` (2-src perm) | 3 | 1.0 | p5 |
| `vpermq` (1-src perm)   | 3 | 1.0 | p5 |
| `vpackssdw` (pack)      | 3 | 1.0 | p5 |
| scalar `mul r64` to hi     | 4 | — | p1 |
| scalar `mulx`           | 4 | 1.0 | p1 |
| scalar `add r64`        | 1 | (4 ALU ports) | p0156 |

## Mixed-port probes (the decisive ones)

| interleave                 | total ops/cyc | conclusion |
|----------------------------|:-------------:|------------|
| IFMA-lo + `vpsllq`         | **1.0** | **shift shares p0 with IFMA -> they CONTEND** |
| IFMA-lo + IFMA-hi          | **1.0** | lo and hi BOTH p0 -> `mulsh52`'s two IFMA serialise |
| IFMA-lo + `vpermq`         | **2.0** | shuffle is p5 -> **table lookup is FREE vs IFMA** |
| IFMA-lo + `vpaddq`         | **2.0** | add takes p5 -> **the `+C` adds are FREE vs IFMA** |

## Derived port model (512-bit execution domain)

- **Port 0 — the scarce resource.** IFMA (`vpmadd52*`), all vector **shifts** (`vpsll/srl/sllv/srlv/rol`), 16-bit/32-bit **multiplies** (`vpmull*`, `vpmulhw`), and **FP mul/FMA/convert** all issue *only* on p0 at 1/cyc. There is **one** 512-bit IFMA/FMA unit on this part.
- **Port 5** — shuffles/permutes/packs/broadcasts (`vperm*`, `vpack*`, `vpunpck*`, `vpbroadcast*`). 1/cyc.
- **Port 0 + Port 5 (2/cyc)** — simple int add/sub, bitwise logic (`vpternlogq`), blends, and compare-to-mask.
- For 512-bit work **p1 is fused into p0** (no separate 256-bit issue), so the classic "p0+p1+p5 = 3 vector ports" of 256-bit code collapses to **2 ports (p0,p5)** at 512-bit width.

**Optimisation rule of thumb for any 512-bit kernel here:** count p0 µops (IFMA + shift + mul + FP) vs p5 µops (perm/pack) vs p0p5 µops (add/logic). The floor is `max(p0_count, p5_count, (p0+p5+p0p5)/2)` cycles. Move work off p0 (into perm/add) or delete it; adds/perms are nearly free until p5 fills.

## Component cycle budget (sign, default x64-exp build, per-sig)

Updated 2026-06-09 mode-120 profile (`avx512/test/profiling.txt`, i7-11700K Rocket Lake, Turbo Boost off, Hyper-Threading off): total **232,731 cyc/op**. `sample_y` is **130,157** (55.9%) with ApproxExp **53,397**, shake8x **45,833**, sigma2 **8,884**, and a_q60+finalize+y **9,946**; `sign_core` is **52,967** (22.8%); rANS pack is **19,763** (8.5%); commitment NTT is **8,137** (3.5%); expand_A0 is **3,975** (1.7%); challenge is **9,875**; setup is **5,965**. Mode 260 total is **488,001 cyc/op**, with `sample_y` **247,807** (50.8%) and ApproxExp **107,262**.

With `AVX512IFMAEXP=1`, mode-120 sign profile total is **188,721 cyc/op**: ApproxExp drops to **11,015** and `sample_y` drops to **87,464**. The 2026-06-11 end-to-end median `speed_sign` is **210,944 -> 166,090** cycles in mode 120.

## exp kernel measurements (mode 120)

- scalar x64 batch (`approx_exp_4x`, B=8 deg-6): **19.56 cyc/exp** (4-wide); 1-lane latency-bound 42.4, 1-lane throughput 20.8.
- IFMA (B=32 deg-5): scalar-Q63 ref 5.927 ns; **W8 1.548 ns (3.83x), W16 1.288 ns (4.60x), W32 1.292 ns (4.59x)**. Accuracy max-rel $2^{-48.67}$ (gate $2^{-43}$).
- **W16 == W32 => 2 independent zmm streams already hide the 4-cyc IFMA latency; the kernel is p0-THROUGHPUT bound, not latency bound.** In `mulsh52` the 2 shifts cost as much as the 2 IFMA (all p0); the table (`vpermt2q`, p5) and `+C` adds (p5) are free.
