# Lithium-Code AVX-512 Implementation — Design Decisions

This document records the design of the AVX-512 IRS implementation in `avx512/`, for review. It covers what was built, why, how it stays KAT-compatible and constant-time, and where every constant comes from. It mirrors `avx2/DESIGN.md`; read that first for the shared scheme background.

## 1. Scope & goals

Optimize the three IRS NIST modes (120/180/260) with AVX-512 while:
- **KAT-compatible** with `ref/` and `avx2/` (byte-identical pk/sk/sig) on the default build;
- **constant-time core arithmetic** (no division, no long-term-secret-indexed branches or memory accesses; keygen retry is one-time provisioning, sign/rANS retry probability is near $2^{-100}$ by parameter design, and Gaussian sampler variance follows isochronous Gaussian sampling rather than exposing attacker-usable information);
- maximizing throughput via five components: AVX-512 NTT, a 16-lane σ2/noise BaseSampler (two independent zmm streams), 8-way SHAKE (keccak8x), 512-bit scheme-glue, and an opt-in **IFMA** ApproxExp.

Benchmark validation host: **i7-11700K Rocket Lake benchmark system (Turbo Boost off, Hyper-Threading off, native Ubuntu)**. The default `avx512` build (integer exp) reproduces the canonical KAT hashes exactly (120 -> `16073584959312876858`, 180 -> `18204602727810236833`, 260 -> `4980115632847006284`) — **shared by ref, avx2 and avx512**.

These hashes are the **re-recorded** canonical set (see §13): the AVX-512-driven sampler rework lives in shared `ref/` and re-records all three backends together. ref/avx2/avx512 are byte-identical.

## 2. Directory layout & file ownership

Same model as `avx2/`: symlink the shared scalar sources from `../ref/` (single source of truth), keep real forks for the files that need vectorization, plus the AVX-512-only files, the per-mode NTT assembly/constants are VENDORED locally in `<qset>/` (copied from `tools/ntt`).

| Category | Files |
|---|---|
| symlink -> `../ref/` (shared scalar) | `reduce.{c,h}`, `ntt.{c,h}`, `rounding.{c,h}`, `packing.{c,h}`, `irs_rans.{c,h}`, `approx_exp.{c,h}`, `fips202.{c,h}`, `symmetric*.{c,h}`, `randombytes.{c,h}`, `params.h`, `polyvec.h`, `sampler.h`, `sign.h`, `irs.h`, `api.h`, `namespace.h`; `test/*.c` shared cases |
| real override (`config.h`) | `config.h` (namespace `_avx512`) |
| **real forked** (vectorized, `#if __AVX512F__ / #elif __AVX2__ / #else`) | `irs.c`, `sign.c`, `polyvec.c`, `sampler.c`.  `sampler.c` now also holds the 16-lane σ2/noise CDT + IFMA `gauss_a_q60` **and** the keccak8x gauss/noise vector backends; `polyvec.c` the AVX-512-NTT routing **and** the shake128x8 uniform — all behind `USE_AVX512_*` |
| real AVX-512 NTT decls + keccak8x | `ntt_avx512_decls.h`, `fips202x8.{c,h}`, `keccakf1600x8.c` |
| real IFMA exp (opt-in) | `approx_exp_avx512ifma.c` (its Q51 constants are an `@@AUTOGEN:approx_exp_avx512ifma@@` region in the same `.c`) |
| per-mode NTT (vendored locally) | `<qset>/{ntt_avx512.S,ntt_consts_avx512.c}` for `qset` in `{q57601n256,q30977n256,q64513n512}` (copied from `tools/ntt`; self-contained build) |
| tests (new) | `test/test_ntt_avx512.c`, `test/test_sigma2_diff.c` (16-lane σ2+noise), `test/test_keccakx8.c`, `test/speed_exp_avx512ifma.c` |

Assembly symbols (`ntt_avx512`, `ntt512_qdata`, …) are global (one mode per binary). `fips202x8` symbols are namespaced `lithium_fips202x8_avx512_`.

## 3. The `ref/` <-> `avx512/` fork contract

`ref/` is the pure scalar master. Each forked `avx512/` `.c` is its vectorized twin, KAT-identical by construction (exact-integer SIMD; the IFMA exp is the only opt-in deviation, §7). The forks (`irs.c`/`sign.c`/`polyvec.c`/`sampler.c`) use a **three-way ladder** `#if defined(__AVX512F__) / #elif defined(__AVX2__) / #else` (scalar). The ladder is mandatory: `-march=native` on the i7-11700K Rocket Lake benchmark system (Turbo Boost off, Hyper-Threading off, native Ubuntu) defines BOTH `__AVX512F__` and `__AVX2__`, and `#if/#elif` selects the 512-bit arm; the AVX2 arm is the verbatim avx2 fork (so a future `diff avx2 avx512` is trivial), the `#else` is the portable reference. Routing macros: `USE_AVX512_NTT` / `_SAMPLER` / `_SHAKE8X` (set in the avx512 CFLAGS) exclude the scalar fallbacks in `polyvec.c`/`sampler.c` and route to the SIMD overrides. The shared opt-in exp switch `AVX512IFMA_EXP` lives in shared `approx_exp.c` (dormant unless `make AVX512IFMAEXP=1`).

## 4. AVX-512 NTT integration

The NTT routing in `polyvec.c` sends the three NTT-domain primitives to the in-tree AVX-512 assembly (`ntt_avx512_decls.h` carries the per-mode prototypes). The AVX-512 forward/inverse transforms use a **uniform** signature `(poly, qdata, ztab, scale)` for all modes (no `z0`/`cross` args, unlike avx2), and `int16_t[...][32]` full-ZMM constant tables (`ntt512_*`).

| mode | q | NTT | pointwise call | reduce |
|---|---|---|---|---|
| 120 | 57601 | incomplete | `basemul_avx512(c, qdata, a, b, gamma)` | none (strict) |
| 180 | 30977 | incomplete | `basemul_avx512(c, a, b, qdata, gamma)` | **none — basemul self-canonicalizes** |
| 260 | 64513 | complete | `pointwise_avx512(c, a, b, qdata)` | none (strict) |

**Mode-180 has NO reduce site** (the key divergence from avx2, which calls `reduce_avx`): the AVX-512 `basemul_avx512`'s `montmul512`/`addmod512`/`dmul512` macros each end in a branchless conditional subtract (`vpcmpuw` + masked `vpsubw`), so the stored output is already `[0,q)` — verified directly from `tools/ntt/q30977n256/ntt_avx512.S` and end-to-end by the KAT. `nttunpack_avx512` lays the Pre-NTT A0 seed-to-NTT samples into the avx512 NTT slot order (the forward shuffle ladder, no butterflies); A0h's internal layout differs from ref/avx2 but is invisible to pk/sk/sig (A0 is regenerated from the seed in keygen+verify), so the KAT is preserved. Validation: `test_ntt_avx512` (poly_mul vs schoolbook + a*a^{-1}==1, 3000+1000 polys/mode) and `t_poly`.

## 5. 16-lane BaseSampler + noise (in `sampler.c`)

`sampler_sigma2` and `noise_magnitude_batch` are 16-lane (zmm, 32-bit lanes) CDT scans over `RCDT_Z` (14/16/20 entries) and `RCDT_NOISE` (10 entries). `GAUSS_BATCH = NOISE_BATCH = 32`, processed as **two independent 16-lane streams** (dual accumulators `za`/`zb`) so the borrow chains fill both ALU ports — the user-specified 2-stream IPC design. Each 16-lane `v0/v1/v2` zmm is assembled (hoisted out of the RCDT loop) from two `vmovdqu ymm` + `vinserti64x4` over the grouped byte layout (`group*96 + lane*4`; v0/v1/v2 at +0/+32/+64), so the layout is byte-identical to the scalar/avx2 path.

The borrow chain is the user-prescribed **`sub,srli / sub,sub,srli / sub,sub,srli / add`** — `cc = (v - Z - borrow_in) >> 31` per limb. With the v limbs masked to 31 bits this is a *true* borrow (bit 31 is the borrow-out; no 32-bit wrap, no sign-flip), which is **exactly what the scalar `ref/sampler.c` already computes** — so the 16-lane kernel is bit-identical to ref by construction. Final: `vpackssdw` (saturation-safe; max count <= 20) + `vpermq {0,2,4,6,1,3,5,7}` to land samples 0..31 in order.

Validation: `test_sigma2_diff` compares both kernels against an in-test scalar reference over 460k random + crafted corner inputs (forced `v0=v1=0` borrow-into-zero, `v2==Z2` boundary) — 0 mismatches, all modes.

## 6. 8-way SHAKE (keccak8x) + the no-waste unified streams

**Keccak-x8** (`fips202x8.{c,h}` + `keccakf1600x8.c`): vendored from XKCP. `keccakf1600x8.c` is the XKCP `KeccakP-1600-times8-SIMD512` permutation **stripped to the 24-round `PermuteAll`** (vpternarylogic+vprolq, fully unrolled), with its local round-constant table deleted and redirected to ref's single `KeccakF_RoundConstants` (so the whole link has exactly one Keccak RC table — `check-consts` asserts this). The byte-accessors, 12/6/4-round variants, FastLoop/K12/Kravatte helpers, and the times4 code were removed; the 5 helper headers folded to none. Re-namespaced `lithium_fips202x8_avx512_`. Validation: `test_keccakx8` (shake128x8/shake256x8 == scalar SHAKE, byte-identical).

**Partition by total coefficients, never idle a lane.** Partitioning the vector samplers by polynomial would idle `8 - ((1+ℓ+m) mod 8)` of the 8 Keccak lanes (K = 5/7/5, never a multiple of 8). Instead both the Gaussian y-vector AND the keygen noise vector are sampled as ONE continuous stream of `N·(1+ℓ+m)` coefficients carved back-to-back from a round-robin of **8 domain-separated sub-streams** `SHAKE256(seed‖nonce‖j)`, j=0..7 (shared `sample_gauss_vec`/`gauss_stream_poly` and the new `sample_noise_vec`/`noise_stream_poly` in `ref/sampler.c`). The 8-block round-robin unit is emitted byte-identically by ref (8 scalar SHAKE), avx2 (2x keccak4x, perf-neutral) and avx512 (1x keccak8x) — so all three KATs converge; the only waste is the final partial 8-block unit (~one unit per vector).

The unified **noise** stream is new (vs the old per-poly `sample_polys_4x`): keygen calls `sample_noise_polys` -> `sample_noise_vec`, and `noise_stream_poly` carves each poly via the grouped, SIMD-friendly `noise_magnitude_batch` (so the 16-lane noise kernel applies). `poly_uniform_8x` (A0) is per-position domain-separated (a public-matrix requirement), so it batches M·ELL positions through keccak8x in groups of 8 (final partial group idles — bounded, keygen/verify-cold) and is KAT-neutral wrt the 4->8 widening.

## 7. IFMA ApproxExp (opt-in, `AVX512IFMAEXP=1`)

Per `docs/design-notes/AVXExp.md` §8: an 8-lane (1 zmm = 8x u64) **Q51 integer** exp via `vpmadd52`. B=32 away-from-0 range reduction (`m=ceil(32a/ln2)+31`, `n=m>>5`, `j=m&31`, `r_l=(m/32)ln2-a`), a degree-5 Q51 Horner (`mulsh52 = (p·X)>>51` = two vpmadd52), a constant-time 32-entry `2^{-j/32}` table scan (mask-blend, no gather), and the `2^{-n}` scale. The polynomial+table are pure integer (IFMA); the range reduction uses binary64 mul/ceil/fmsub (no division; data-operand-independent timing per Intel DOIT / Arm DIT). The accept-test value is produced **pure-integer**: `ps in [2^51, 2^52)` => `exp_q63 = ps << (12 - n)`, exact for `n <= 12`; n <= 11 over the whole sampler domain (asserted by the consts gate), so the `n>12` right-shift / `n>=76` zero branches are an unreachable total-function guard.

**Hard accuracy requirement: max relative error `< $2^{-43}$`** (the IRS accuracy bound). Measured **$2^{-48.67}$** on the actual kernel output (a_q60 -> Q63) vs libm over [0, 7.369] — ~5.7 bits of margin — double-gated: (a) the consts generator's `--check`, and (b) the runtime `speed_exp_avx512ifma` accuracy gate; both build-fail below $2^{-43}$.

**KAT**: separate, opt-in (`KAT_avx512ifma.txt`, `make check-kat-avx512ifma`). The Q51 deg-5 result is ~$2^{-48}$ from the scalar Q63 deg-6 kernel (can flip a rare boundary accept); for the recorded N=10 vectors it matches the shared KAT exactly. **TODO**: full bit-exact reconciliation with the scalar Q63 kernel (the integer-output structure already makes that the only remaining gap). Default `avx512` build (no flag) uses the x64 integer Q63 batch exp (matches the shared KAT) — NOT the avx2 float.

**Table lookup is a `vpermt2q` register permute, not a linear scan.** `s = TAB_Q51[j]` ($j \in [0,31]$) is two `vpermt2q` over the 4 zmm holding the 32 entries, blended by bit 4 of j (~4 ops) — vs the naive 32x(`vpcmpeqq`+`vpblendmq`) = ~64-op constant-time scan, which was the kernel's throughput wall.  `vpermt2q` is a register shuffle (data-independent timing, no memory gather), so it stays constant-time; j is public-range.  This dropped the kernel from ~2x to ~4.8x the scalar Q63 batch and let the multi-stream ILP show.

Microarch sweep (`speed_exp_avx512ifma`, the i7-11700K Rocket Lake benchmark system (Turbo Boost off, Hyper-Threading off, native Ubuntu)): W8/W16/W32 = 1/2/4 independent 8-lane streams = 1.54 / 1.31 / **1.26 ns/exp**, vs 6.00 for the scalar Q63 batch (W32 ~ **4.76x**).  Production is W32 (the GAUSS_BATCH=32 call loops four 8-lane chunks; the OoO engine overlaps them).  In the updated `crypto_sign` profile (mode 120, no Turbo/HT) this cuts ApproxExp 53.4k -> 11.0k cyc/op.

## 8. Scheme-glue 512-bit widening (KAT-neutral)

The avx2 `#ifdef __AVX2__` exact-integer fast paths were widened to 512-bit zmm arms (`#if __AVX512F__`), all KAT-neutral (exact integer; verified — KAT unchanged) and constant-time (mask arithmetic, public loop bounds, N divisible by 16 => no masked tail):
- `irs.c`: `inner` (`vpmuldq` 16xi32->8xi64 + `_mm512_reduce_add_epi64`), `step` (single `vpternlogd` 0xE4 select), `poly_cyclic_shift` (16-wide load/store/negate, scalar tail on the public index `j`).
- `sign.c`: `highbits16` (`vpcmpgt`-mask + masked sub), the makehint finalize loop (masked +/-2q + centermod).
- `polyvec.c`: `poly_sqnorm` (centering + `vpmuldq` squared norm), `mat_mul_2q`/`mat_mul_z1_2q` finalize (parity via `vptestmd`, masked +/-2q). rANS (`irs_rans.c`) and bit-packing (`packing.c`) stay scalar by design (serial state / horizontal).

## 9. Build system

`make` (default `AVX512IFMAEXP=0`, integer exp, shared KAT). Flags: `-O3 -mavx512f -mavx512bw -mavx512vl -mavx512dq -mavx512ifma -mbmi2 -mpopcnt -march=native`, `-DUSE_AVX512_{NTT,SAMPLER,SHAKE8X}`. Per-mode NTT asm appended by explicit path. Targets: `check` (14 tests x 3 modes = 42 PASS), `check-kat` (assert the shared hashes), `check-kat-avx512ifma`, `check-refill` (`-DGAUSS8X_INIT=8 -DNOISE8X_INIT=8`), `check-consts`, `ct-scan`, `profile`, `make AVX512IFMAEXP=1 …`. The build is self-contained (the per-mode NTT is vendored in `<qset>/`); only `check-consts` reaches into `tools/ntt` to re-validate the vendored copies.

## 10. Constant-time stance

- **No division/sqrt** (`make ct-scan`, authoritative): disassembles the secret-handling objects (`sampler.c`, `polyvec.c`, `approx_exp`, `fips202x8`, `keccakf1600x8`, the forked `irs/sign`, the IFMA exp, and all three `ntt_avx512.S`) and asserts no `idiv`/`div*`/`vdiv*`/`vsqrt*` (and no `vgather` in the IFMA table scan) — passes.
- **No secret-dependent branches/gathers**: 16-lane σ2/noise and the NTT assembly are straight-line SIMD; accept/zero-fold/sign are branchless masks; the IFMA `2^{-j/32}` table is a full 32-entry mask-blend scan (not an indexed gather); per-lane shifts use public-range `n`. The IFMA range-reduction float ops (mul/ceil/fmsub on normals) are data-operand-independent per Intel DOIT / Arm DIT.

## 11. Performance (2026-06-11 i7-11700K Rocket Lake, Turbo Boost off, Hyper-Threading off, native Ubuntu)

End-to-end median cycles from `test/speed.txt` (`make speed_txt`). The AVX2 build uses `-O3 -march=native -mno-avx512f`; the AVX-512 build uses `-O3 -march=native` with AVX-512F/BW/VL/DQ/IFMA enabled. Speedups are AVX-512 over default AVX2 from the same no-Turbo/no-HT native-Ubuntu run.

| mode | op | avx2 | avx512 | speedup |
|---|---|---:|---:|---:|
| 120 | keygen / sign / verify | 304k / 337k / 65k | 192k / 211k / 57k | 1.59x / 1.60x / 1.14x |
| 180 | keygen / sign / verify | 302k / 547k / 98k | 214k / 344k / 83k | 1.41x / 1.59x / 1.18x |
| 260 | keygen / sign / verify | 620k / 689k / 118k | 403k / 442k / 103k | 1.54x / 1.56x / 1.14x |

With the opt-in IFMA exp (`AVX512IFMAEXP=1`), sign improves further to **166k / 284k / 361k** for Lithium-120/180/260. That is **2.03x / 1.93x / 1.91x** faster than default AVX2 sign, and **5.85x / 5.59x / 5.66x** faster than scalar `ref` sign in the same benchmark run. Keygen/verify are effectively unchanged by IFMA exp apart from run noise because they do not use the Gaussian exp.

The updated profile (`test/profiling.txt`) shows the mode-120 default sign path dominated by `sample_y` at 130.2k cyc/op (55.9%), followed by `sign_core` at 53.0k (22.8%) and rANS pack at 19.8k (8.5%). Inside `sample_y`, scalar integer ApproxExp is 53.4k, Keccak8x is 45.8k, `a_q60+finalize+y` is 9.9k, and BaseSampler is 8.9k. Keygen is dominated by `sample_noise` 136.5k plus `pack pk/sk` 55.2k and `compute_b` 53.9k; verify is spread across unpack 23.6k, setup 11.8k, challenge 10.0k, usehint 5.7k, commitment NTT 4.3k, and `A0` expansion 3.9k. Full tables are in [`PROFILING.md`](PROFILING.md).

Compared with AVX2, AVX-512 wins most clearly on sign because the 8-way Gaussian SHAKE, 16-lane CDT/finalize glue, and AVX-512 NTT have enough work to amortize the wider vectors. Verify has less Gaussian work and remains about 1.14--1.18x faster than AVX2 in this Rocket Lake run, even with Turbo Boost and Hyper-Threading disabled.

## 12. As-built deviations from the plan (all improvements / corrections)

- **Mode-180 needs no reduce site** — verified the AVX-512 `basemul_avx512` self-canonicalizes (montmul/addmod/dmul each end in a conditional subtract), so `poly_pointwise_montgomery` for mode 180 drops the `reduce_avx` the avx2 path needs. Cleaner; validated end-to-end.
- **`keccakf1600x8.c` reuses ref's RC** rather than carrying the XKCP local table — keeps the single-`KeccakF_RoundConstants` `check-consts` gate green with no exception, and is one fewer constant table to drift.
- **`step` uses a single `vpternlogd` (0xE4)** for the branchless `m ? z-v : z+v` select (vs the avx2 and/andnot/or triple).
- **IFMA exp output is pure-integer** (`ps << (12-n)`) rather than a double->Q63 conversion — keeps the rejection-test value an integer function of the Q51 result; the only float is the range reduction.
- **avx2 stays perf-neutral** under the shared rework: σ2 wraps a 2-iteration outer loop (no extra unroll), and the unified gauss/noise streams emit the 8-block unit via 2x keccak4x — same per-block cost as the old 4-block unit.

## 13. KAT re-record (shared, user-approved)

The AVX-512-driven sampler rework is in shared `ref/sampler.{c,h}` + `ref/sign.c`, so ref/avx2/avx512 re-record together to one canonical KAT. Three layout changes drive it: (1) `GAUSS_BATCH` 16->32 (the two-16-lane-stream σ2 width; re-groups the mini-batch bytes), (2) the unified Gaussian y-stream widened 4->8 sub-streams, (3) the unified no-waste key-noise stream + grouped `noise_magnitude_batch`. New canonical hashes recorded in `ref/vector/KAT.txt`; `avx2/KAT_avx2exp.txt` updated. The re-record was validated step-by-step with ref<->avx2 convergence at each change, then avx512 confirmed byte-identical.

## 14. Provenance / regeneration

- AVX-512 NTT asm + constants: `tools/ntt/<qset>/gen_ntt.py` (one run emits both the avx2 and avx512 artifacts; `make check-consts` regenerates + byte-exact-validates).
- Keccak-x8: XKCP (https://github.com/XKCP/XKCP); banner + modification list in `keccakf1600x8.c`.
- IFMA exp Q51 coeffs/table: `tools/avxexp/gen_avxexp_ifma_consts.py` (Sollya fpminimax fit, same as the avx2 float coeffs; `--check` asserts max-rel-error < $2^{-43}$ and n <= 12).

## 15. Post-review optimizations (sign bottleneck round)

Profiling sign (mode 120) after the initial port showed the only big lever — SHAKE (keccak8x, 2.5x) — was largely offset by an unaccelerated ApproxExp, an expand_A0 regression, and serial tails.  Three targeted fixes:

1. **IFMA exp table lookup (the big one).** The Q51 kernel's hidden cost was the 32-entry constant-time table scan (~64 ops/8-lane).  Replaced by two `vpermt2q` register permutes + a blend (§7): kernel 2x -> **4.76x** the scalar batch; with `AVX512IFMAEXP=1`, the current no-Turbo/no-HT run shows IFMA sign **1.86--2.00x** faster than default AVX2 (§11).

2. **expand_A0 as 8 equal coefficient chunks (no idle lanes).**  The per-position 8-way left modes 120/260 (M·ELL=4 < 8) falling to scalar — a regression vs avx2.  Now A0 is sampled as 8 EQUAL chunks of M·ELL·N/8 coeffs from `SHAKE128(rho‖j)`, j=0..7 (`poly_uniform_chunkx8`; ref=8 scalar, avx2=2xkeccak4x, avx512=1xkeccak8x), so the 8-way engine is ALWAYS full for every mode.  expand_A0 23.5k -> **10.4k** cyc; keygen 1.7–1.95x -> **1.9–2.2x**, and verify (A0-dominated) flipped from ~0.95x to ~1.0–1.07x.  Re-records the KAT (the per-position->chunk byte layout changes A0; all three backends converge — see §13).

3. **IFMA `gauss_a_q60` (gauss_finalize exploration).**  Microbench: `gauss_a_q60` is 1.316 ns/call vs `gauss_finalize` 0.405 — the finalize bucket is dominated by the `num·I80>>20` 64x64->128 multiply, not the compaction.  Since `SAMP_I_HI==0` (I80 = 64-bit I_LO), splitting I80 = H·$2^{52}$+L gives an exact `a_q60 = (mid<<32)+(loL>>20)` (mid = hiL + num·H) in ~5 `vpmadd52` ops (`gauss_a_q60_batch`, 8-lane, KAT-neutral).  AVX2 kept this scalar (it has no 52x52->104 multiply — its prior SIMD attempt was a wash); IFMA does it natively.  In context the win is modest (~3k off sample_y): the scalar `mulx` is already pipelined and OoO-hidden, so the real floor of the finalize bucket is the serial per-candidate loop (rand-tail load + the data-dependent `r[coefcnt++]` compaction), not `a_q60` throughput.  Further gains there would need a `vpcompressw` (AVX512-VBMI2) compaction restructure — intricate, KAT-sensitive, and uncertain — left as future work.

KAT after this round (re-recorded for the chunk-A0 change, shared ref==avx2==avx512): 120 -> `16073584959312876858`, 180 -> `18204602727810236833`, 260 -> `4980115632847006284`.

## 16. Microarchitecture-driven optimization campaign (port-model round)

Instruction latency/throughput audit on one i7-11700K Rocket Lake native-Ubuntu test system with Turbo Boost and Hyper-Threading disabled via `tools/uarch/ubench.c` (rdtsc + `vpaddq`-latency anchor = frequency-normalised core cycles; the test environment had no usable PMU). The decisive facts (`tools/uarch/RESULTS.md`): **one** 512-bit IFMA unit on **port 0**, and IFMA + every vector **shift** + 16/32-bit **multiply** + FP convert all issue 1/cyc on p0 and *contend* (proven by mixed probe: IFMA-lo+`vpsllq` = 1.0 ops/cyc, while IFMA+`vpermq`/`vpaddq` = 2.0). Shuffles are p5; add/logic are p0+p5 (2/cyc). So the optimization rule is **shed port-0 work** (delete shifts, mask-compares instead of sub+srli, `vpmaddwd` instead of `vpmuldq`+shifts) or **vectorize scalar serial tails**. A per-component IPC audit (ideal port-floor vs measured) drove six changes:

1. **ApproxExp IFMA** — `mulsh52` kept hi-only (`madd52hi<<1`, drops the bit-51 low-product carry; acc $2^{-48.2}$, still >5 bits past the $2^{-43}$ gate) + deleted the provably-dead `n>12` output guard. Halves the kernel's p0 footprint: W16 1.29->0.80 ns/exp (1.62x), W32 0.74 ns (8.0x vs scalar). KAT-neutral.
2. **sigma2/noise CDT** — `vpcmpltud` + masked `vpaddd` mask-borrow chain (no p0 shift) replacing sub+srli. 8,530->6,556 cyc/sig. Bit-identical.
3. **expand_A0** — `vpcompressd` uniform-mod-q rejection (16-wide) replacing the scalar branchful `o[ctr++]`. 9,818->7,258 cyc/sig (m120); also speeds keygen/ verify. Bit-identical.
4. **gauss finalize + sampler_y** — `vpcompressw` accept/compaction (prefix-stop + output-indexed signs via pdep/tzcnt) and the 7-bit `sampler_y` unpack as `vpermb`+`vpsrlvw`. a_q60+finalize+y 26,944->8,040 cyc/sig (3.35x). Bit-identical.
5. **NTT input conversion — fixed; forward scheduling still deferred.** The large measured gap around forward NTT was partly a C-side input-conversion artifact, not just shuffle-to-butterfly scheduling inside the NTT body: `poly_ntt` used a scalar `freeze()` prefix and `poly_ntt_small` used a scalar small-signed mask-add prefix. The AVX-512 qset `ntt_avx512.S` files now expose three standard-order conversion helpers: full `poly_freeze_to_u16_avx512` (exact Barrett `freeze()` with the mode-180 33-bit reciprocal split), small-signed `poly_freeze_small_to_u16_avx512`, and `[0,q)` `poly_cast_modq_avx512`. `poly_ntt`/`poly_ntt_canonical` use the full helper, `poly_ntt_small` uses the small helper, and `compute_b` uses cast for `finv` plus small for the reduced `(A0s0+e)` accumulator. These helpers do not change AVX slot layout or canonical `bn = NTT_ref(b)` serialization; `test_freeze_avx512`, KAT, and ct-scan cover the change. Further butterfly scheduling work remains deferred because the assembly body itself is near its port floor.
6. **rANS — 2-way interleaved** (`irs_rans.c` + `tools/rans.py` mirror). The single-state encode is a per-symbol reciprocal-multiply dependency chain; 2 interleaved states hide the ~5-cyc multiply (isolated bench mode120/1280 symbols: 10080->7304 cyc, 1.38x; N>=4 adds nothing but flush bytes, so N=2). **Wire-format change** -> KATs were re-recorded at that time; the later true-rANS-PMF update re-recorded the current canonical KATs listed at the top of this file. +4 flush bytes at that time; the current true-rANS-PMF reserve uses the two-state flush in `tools/SigSize.py`.

Current end-to-end mode-120 medians are 210,944 cycles for the default AVX-512 sign path and 166,090 cycles with IFMA exp; this is 1.60x and 2.03x faster than default AVX2 sign in the same run. Keccak-x8 was measured **at** its port floor (47.9/48 cyc/round) and left untouched.

7. **rANS *decode* — O(1) slot LUT (shared, verify ~3x).** Separate from the encode interleave above: the decoder's per-symbol `find_slot` was an O(n) linear CDF scan — ~74 % of verify, and (being scalar in shared `irs_rans.c`) the reason avx512 verify barely beat avx2. Replaced by a `SLOT[1024]` `val->slot` table (`gen_rans_tables.py`; KAT-neutral, no wire change). **current mode-120 avx512 `unpack` is 23.6k cyc/op and verify median is 57.5k** — and ref/avx2 gain the same. With the rANS floor gone, verify is now NTT-mat-mul / setup bound.

8. **`poly_invert` — Montgomery batch inversion (shared, keygen).** `poly_invert` inverted `f` via per-element Fermat (`a^{q-2}` modexp): 512 residues (full NTT, mode 260) or 128 pair-determinants (incomplete, 120/180) — the bulk of `compute_b`. Replaced by Montgomery's trick (`batch_modq_inv` in shared `polyvec.c`): one Fermat + ~3n modmuls for the whole vector, each result the unique inverse -> KAT-neutral; ref/avx2/avx512 all gain. **current mode-120 `compute_b` is 53.9k cyc/op** (mode 120); largest at mode 260. Residual: the 2 scalar NTTs in `poly_invert` (scalar Cooley-Tukey coeff order).

Full per-component cycle profiles (keygen / sign / verify, sorted by time-share): [PROFILING.md](PROFILING.md).
