# one AVX2-capable benchmark system — 256-bit AVX2 port model

Companion to `RESULTS.md` (which covers the 512-bit domain). Source: **uops.info** per-instruction pages, **the recorded CPU generation (RKL)** rows (identical to Ice Lake / the tested core family — same core). Latency in cycles; throughput in **instructions/cycle**; ports per µop.

## Measured table (256-bit ymm unless noted)

| instruction | lat | TP (insn/cyc) | µops | ports |
|---|:--:|:--:|:--:|---|
| `vfmadd231pd/ps` (FMA) | 4 | **2.0** | 1 | **p01** |
| `vmulpd` / `vaddpd` | 4 | 2.0 | 1 | p01 |
| `vroundpd` (ceil) | 8 | 1.0 | 2 | 2xp01 |
| `vcvtpd2dq` (ymm->xmm i32) | 7 | 1.0 | 2 | p01 + **p5** |
| `vpmovsxdq` (i32->i64 widen) | 3 | 1.0 | 1 | **p5** |
| `vblendvpd` (variable blend) | 2 | 1.0 | **2** | 2xp015 |
| `vpcmpeqd` / `vpcmpgtd` | 1 | 2.0 | 1 | **p01** |
| `vpaddd`/`vpsubd`/`vpand`/`vpor`/`vpxor` | 1 | **3.0** | 1 | **p015** |
| `vpmuldq` (32x32->64) | 5 | 2.0 | 1 | p01 |
| `vpmulhuw`/`vpmullw` (16-bit mul) | 5 | 2.0 | 1 | p01 |
| `vpermd`/`vpermps` (cross-lane perm) | 3 | 1.0 | 1 | **p5** |
| `vpackssdw`/`vperm2i128` (pack/perm) | ~3 | 1.0 | 1 | p5 |

## Derived port model (256-bit) — THREE vector ports

- **p0, p1** — the two FP/"complex-int" ports: FMA/mul/add (**2/cyc**), 16/32-bit integer **multiply** (`vpmul*`, 2/cyc), **compares** (`vpcmpeqd`/`vpcmpgtd`, 2/cyc), vector shifts, and part of the FP converts.
- **p5** — the shuffle port: `vperm*`, `vpack*`, `vpunpck*`, **`vpmovsxdq`/widening converts**, broadcasts (1/cyc), plus the p5 half of `vcvtpd2dq`.
- **p0+p1+p5 (p015, 3/cyc)** — simple int add/sub, bitwise logic, immediate blend, and the 2-µop variable `vblendvpd` (spreads over all three).
- Loads p2+p3 (2/cyc), stores p4(+p7).

**Contrast with the 512-bit domain (`RESULTS.md`):** at 512-bit, p0+p1 fuse into a single unit -> only **2** vector ports and **1 FMA/cyc**; at 256-bit there are **2 FMA ports (2/cyc)** and **3 int/logic ports (3/cyc)**. So AVX2 has more issue parallelism per op (at half the width), and the rule shifts: the scarce resources are **p5** (every shuffle/permute/widen) and the **p01 FMA/mul/compare** pair — *not* a single fused port.

**Optimization rule for AVX2 kernels here:** count µops by port — p01 (FMA/mul/cmp/shift), p5 (perm/pack/widen), p015 (add/logic/blend). Floor = `max(p01/2, p5/1, ceil(total/3))` cycles, given enough ILP to hide the 4–5-cyc FMA/mul latency. Shed **p5** traffic (per-element widens, permutes) and **p01** traffic (compares/muls/shifts); plain adds/logic are cheap (3/cyc).

## Application — AVX2Exp (`avx2/approx_exp_avx2.c`), per 4-lane stream, B=32

The 32-entry constant-time `vblendvpd` scan (lines 50-56) per stream:
- 32x `vpcmpeqd` -> **32 p01** µops
- 32x `vpmovsxdq` (widen the 32-bit compare mask to 64-bit) -> **32 p5** µops <- a hard p5 bottleneck
- 32x `vblendvpd` (2 µops, sequential dependency on `s`) -> **64 p015** µops, **64-cyc latency chain**

vs the essential work: range reduction (`vmulpd`+`vroundpd`+`vaddpd`+`vcvtpd2dq`+`vfmsub` approx 6 p01/p5 µops) + degree-5 Horner (**5 `vfmadd`**, p01) + scale (~6 µops). So the scan is **~5–6x the essential cost**, p5-bound (the 32 widens) and latency-bound (the 32-deep blend chain) — confirming "far from ideal IPC."

**Fixes (Phase 3):**
- **F-a (KAT-neutral): masked-OR tree.** Widen `j` to 64-bit **once**, then per entry `vpcmpeqq` (p01) + `vpand` (p015) into an independent accumulator, combined by a depth-5 `vorpd` tree. Removes the 32 per-entry `vpmovsxdq` (frees p5) and the 32-deep `vblendvpd` chain (-> depth-5 tree). Bit-identical (`OR_k(T[k]&mask_k) == T[j]`). ~1.3x kernel.
- **F-b (re-records the opt-in avx2-float KAT): B=32 -> 8, deg-6.** AVX2 has no cheap 32-double permute (unlike AVX-512's `vpermt2q`), so the scan length is the lever: B=8 cuts it 4x (8 entries) for one extra Horner FMA. Brings the kernel near the FP-bound floor (~2.4 cyc/exp). The big throughput win (~2.8x); keep whichever the `speed_exp_avx2` sweep shows best.
