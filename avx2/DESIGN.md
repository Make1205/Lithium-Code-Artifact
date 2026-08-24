# Lithium-Code AVX2 Implementation — Design Decisions

This document records the design of the AVX2-optimized IRS implementation in `avx2/`, for review. It covers what was built, why, how it stays KAT-compatible and constant-time, and where every constant comes from.

## 1. Scope & goals

Optimize the three IRS NIST modes (120/180/260) with AVX2 while:
- **KAT-compatible** with `ref/` by default (byte-identical pk/sk/sig);
- **constant-time core arithmetic** (no division on secret data, no long-term-secret-indexed branches or memory accesses; keygen retry is one-time provisioning, sign/rANS retry probability is near $2^{-100}$ by parameter design, and Gaussian sampler variance follows isochronous Gaussian sampling rather than exposing attacker-usable information);
- maximizing throughput via four components: AVX2 NTT, 8-way BaseSampler, 4-way SHAKE, and an opt-in float ApproxExp.

AVX-512 is intentionally out of scope (AVX2 only).

Result: the default AVX2 build reproduces the three current reference KAT hashes exactly (120 $\to$ `16073584959312876858`, 180 $\to$ `18204602727810236833`, 260 $\to$ `4980115632847006284`). These are the hashes after the shared sampler rework and the true-rANS-PMF reserve update; ref, avx2, and avx512 are byte-identical in the default integer-exp build.

## 2. Directory layout & file ownership

**`ref/` is a pure scalar reference — no SIMD, no `#ifdef __AVX2__`, no `USE_AVX2_*` routing macros.** `avx2/` **symlinks the shared scalar sources from `../ref/`** (single source of truth), keeps **real forked copies** of the files that need vectorization, plus the AVX2-only files, and references the per-mode NTT assembly/constants **by explicit path** from `../tools/ntt/`.

| Category | Files |
|---|---|
| symlink $\to$ `../ref/` (shared scalar) | `reduce.{c,h}`, `ntt.{c,h}`, `rounding.{c,h}`, `packing.{c,h}`, `irs_rans.{c,h}`, `approx_exp.{c,h}`, `fips202.{c,h}`, `symmetric*.{c,h}`, `randombytes.{c,h}`, `params.h`, `polyvec.h`, `sampler.h`, `sign.h`, `irs.h`, `api.h`, `namespace.h`; `test/*.c` (incl. `../ref/test/prof.{c,h}`) |
| real override (`config.h`) | `config.h` (namespace `_avx2`) |
| **real forked** (pure scalar in `ref/`, **vectorized here**) | `irs.c`, `sign.c`, `sampler.c`, `polyvec.c`.  `sampler.c` now also holds the 8-way σ2 + key-noise CDT **and** the keccak4x 4-way Gaussian/noise backends; `polyvec.c` the SIMD-NTT routing **and** the shake128x4 uniform — all behind `USE_AVX2_*` |
| real AVX2 NTT decls + XOF | `ntt_avx_decls.h`, `fips202x4.{c,h}`, `f1600x4.S` |
| real AVX2 exp (opt-in) | `approx_exp_avx2.c`, `approx_exp_avx2.h` (the float-exp constants are an `@@AUTOGEN:approx_exp_avx2@@` region in the `.h`) |
| per-mode NTT (vendored locally) | `<qset>/{ntt.S,ntt_consts.c}` for `qset` $\in$ `{q57601n256,q30977n256,q64513n512}` (copied from `tools/ntt`; self-contained build) |
| tests (new) | `test/test_ntt_avx.c`, `test/test_sigma2_diff.c`, `test/speed_exp_avx2.c` |

The four forked `.c` files keep the `USE_AVX2_NTT` / `USE_AVX2_SAMPLER` / `USE_AVX2_SHAKE4X` guards (set in the avx2 `CFLAGS`) to exclude the scalar fallbacks and route to the SIMD overrides; their `ref/` twins have those branches stripped.  The one AVX2 switch that lives in *shared* code is the opt-in float-exp dispatch `AVX2_EXP` in `approx_exp.c` — dormant unless `make AVX2EXP=1`.

The assembly symbols (`ntt_avx`, `ntt_qdata`, …) are **global** (no namespace); a binary links exactly one mode, so the shared names never collide. A single process therefore cannot host multiple modes — acceptable for the per-mode-binary model the test suite uses.

## 3. The `ref/` $\leftrightarrow$ `avx2/` fork contract

`ref/` is byte-for-byte the scalar reference (its KAT is the master); each forked `avx2/` `.c` is its vectorized twin, KAT-identical by construction (exact-integer SIMD; the float exp is the only opt-in deviation). What each twin changes:
- `polyvec.c`: the three NTT-domain primitives (`poly_ntt`, `poly_invntt_tomont`, `poly_pointwise_montgomery`) routed to the in-tree AVX2 assembly under `USE_AVX2_NTT`; `poly_uniform_4x` to shake128x4 under `USE_AVX2_SHAKE4X`; `poly_sqnorm` to an exact int64 `vpmuldq` path under `__AVX2__`. `poly_invert` stays on the **scalar** transforms (its residue inversion `invert_pair16` / `zetas[half+i]` is tied to the scalar Cooley-Tukey order).
- `sampler.c`: scalar `sampler_sigma2` excluded under `USE_AVX2_SAMPLER` (the 8-way version, now in the same file below, replaces it); the 4-way wrappers + ref combined backend excluded under `USE_AVX2_SHAKE4X` (keccak4x versions, also in-file). The 8-way σ2/noise kernels and the keccak4x backends were consolidated **into** `sampler.c`. The RCDT tables are `extern const` and shared from `sampler.c`.
- `irs.c`: `inner` (exact `vpmuldq` dot product), the masked `z+/-v` update, and `poly_cyclic_shift` get `#ifdef __AVX2__` fast paths.
- `sign.c`: the HighBits/MakeHint loops get `#ifdef __AVX2__` fast paths.
- `approx_exp.{c,h}` (shared, NOT forked): `exp_batch()` is the single switch point between the integer backend (default, bit-identical to ref) and the float backend (`approx_exp_avx2`, under `AVX2_EXP`).
- 4-way batch scheduling (`sign.c` y-loop / keygen noise-loop, `polyvec.c` `expand_A0`) is identical in both twins — same seeds/nonces, pure grouping.

## 4. AVX2 NTT integration

The NTT primitives in `polyvec.c` are per-mode `#if`-switched wrappers (`ntt_avx_decls.h` carries the exact NTT prototypes so the compiler type-checks the **mode-specific basemul argument order**):

| mode | q | NTT | pointwise call | reduce |
|---|---|---|---|---|
| 120 | 57601 | incomplete | `basemul_avx(c, qdata, a, b, gamma)` | none (strict) |
| 180 | 30977 | incomplete | `basemul_avx(c, a, b, qdata, gamma)` | `reduce_avx(c)` |
| 260 | 64513 | complete | `pointwise_avx(c, a, b, qdata)` | none (strict) |

**Lazy-reduction contract (q30977)**: forward/basemul outputs are lazy `[0,2q)` (2q = 61954 fits uint16). The only elementwise ops on NTT-domain data in `polyvec.c` (`subm16(0,·)` negation, `poly16_add`) require `[0,q)` and run *only* on basemul outputs — so a **single** `reduce_avx` inside the pointwise wrapper canonicalizes them. q57601/q64513 are strict (2q overflows uint16), so values are always `[0,q)`; no reduce needed. `poly16` (`uint16_t[N]`) and the assembly's `int16_t*` are bit-compatible; no buffer alignment is required (the asm uses unaligned loads; only the `aligned(32)` constant tables use aligned loads).

Validation: `test_ntt_avx` (poly_mul vs schoolbook, 3000 polys/mode + 1000 inversions) and `t_poly` pass for all three modes.

## 5. 8-way BaseSampler (in `sampler.c`)

8-way `sampler_sigma2` over each 16-sample unit (two groups of 8, dual accumulators for ILP); the canonical `GAUSS_BATCH=32` mini-batch runs two such units and reuses the byte layout the scalar ref already uses.

**Comparison chain**: the RCDT limbs are $\le 31$ bits (loaded `v` masked to `0x7FFFFFFF`, `RCDT_Z` limbs 31-bit), so the exact multi-precision borrow is a plain signed subtract whose sign bit (bit 31) is the borrow-out: `borrow0 = (v0 <u Z0)` (signed `cmpgt` suffices, no flip), then `borrow_k = bit31(v_k - Z_k - borrow_{k-1})`. `v - Z - borrow` lands in $[-2^{31}, 2^{31})$ so bit 31 is the true sign — no uint32 wrap, no sign-bit flip. This is the **exact** 93-bit unsigned compare `(v0,v1,v2) <u (Z0,Z1,Z2)`; per group per RCDT row it is `cmpgt, add, sub, srli, sub, sub, srli` (about half the ALU ops of a sign-flip unsigned-`<` chain). The ref scalar chain (`sampler.c`) is the identical borrow. RCDT values are broadcast straight from the canonical `RCDT_Z` table, so there is no separate broadcast table to drift.

History: an earlier design used the scalar `ct_lt_u32(v_mid - borrow, Z_mid)` form, which has a uint32 wrap at `v_mid=0, borrow=1` (`0-1 = 0xFFFFFFFF` is not `< Z`, so the borrow is dropped) — a $\sim 2^{-31}$ under-borrow vs the exact compare. Both ref and avx2 now use the exact bit-31 borrow above; this is KAT-neutral on the committed vectors (no sample hits the `v_mid=0`-with-borrow corner) while removing the quirk.

Validation: `test_sigma2_diff` compares the AVX2 kernel against an in-test scalar reference (same exact borrow) over 200k random + 200k crafted corner inputs (forced `v_mid=0` with borrow) plus a `v2==Z2` equality-boundary case — 0 mismatches, all modes. `sampler_y` / `gauss_finalize` stay scalar (cheap, branchless).

## 6. 4-way SHAKE (`fips202x4`; 4-way samplers consolidated into `sampler.c` / `poly_uniform_4x` into `polyvec.c`)

`fips202x4.c` + `f1600x4.S` vendored from IRS, namespace `lithium_fips202x4_avx2_`; they link directly against ref's `KeccakF_RoundConstants` (standard 24-round order). The 4-way samplers run 4 independent keccak4x sponges with a lockstep squeeze + per-lane cursors; each lane reproduces the scalar per-lane byte schedule exactly (40-byte `seed‖nonce_LE8` absorb, 0x1F pad, $\text{signs} \to \sigma2 \to y \to \text{tails}$), so per-lane output is byte-identical. Idle/finished lanes have their buffers reset so the unconditional 4-way writes stay in bounds. `poly_uniform_4x` is the standard shake128x4 rejection expansion, one block $\times 4$ per round.

Validation: `test_sampler` (4x==1x), `test_keys`, and the full KAT pass.

### 6.1 Unified single-stream Gaussian `y` (near-zero squeeze waste)

The Gaussian `y`-vector is $K = 1+\ell+m$ polys (5/6/5). Naively running keccak4x as "4 lanes $\to$ 4 polys" wastes SHAKE two ways: (a) idle padding lanes for the K-mod-4 remainder, and (b) lockstep **uneven completion** — the 4 lanes consume different amounts (rejection variance) but keccak4x squeezes all 4 to the slowest, wasting the rest. Both are eliminated by sampling **all K polys from ONE continuous stream**: `sample_gauss_vec` runs a single keccak4x over four domain-separated sub-streams lane *j* = `SHAKE256(seed‖nonce_LE8‖j)`, j=0..3 (41-byte input $\to$ domain-separated from the 40-byte single-nonce streams, so no collision with any other poly/attempt), and samples the K polys **back-to-back** from a persistent buffer (`gauss_stream`) whose leftover carries from one poly to the next. The only waste is the final partial 4-block unit after the last poly — squeeze is $\sim 1.006\times$ consumed (effectively zero).

The KAT-critical loop is the shared `gauss_stream_poly` (`sampler.c`): it reads round-robin 4-block units (one keccak4x perm = `[lane0 blk][lane1][lane2] [lane3]`) from a `squeeze4` callback and runs the exact `sample_gauss_N` math, carrying the cursor across polys. Only the byte source differs — ref uses 4 scalar SHAKE states, avx2 uses keccak4x — and both emit the **identical** stream, so ref and avx2 are bit-identical by construction. (Noise/keygen keeps the simple per-nonce 4-way + padding — it's cold and KAT-unaffected; `sample_gauss_N_4x` remains for `test_sampler`'s 4x==1x check.)

This scheme changed the canonical KAT (each poly's randomness is now a slice of the 4-lane combine, not its own single stream) — an explicit, user-approved re-record; ref and avx2 share the new hashes. Net effect vs the original per-poly-single-stream avx2: Gaussian SHAKE squeeze cut ~45 % and sign ~11–19 % (§11). No `GAUSS4X_INIT`-style tuning is needed for `y` (the stream refills in 4-block units on demand); `sample_gauss_N_4x`'s per-mode `GAUSS4X_INIT` (56/68/104) only affects the now-test-only per-nonce path.

### 6.2 Pre-NTT: `A0` seed-to-NTT + cached `b̂` (Dilithium-style)

`A0` is a uniform, public matrix used in every `A·y` (sign) and `A1·z1` (verify) product, so — exactly like Dilithium's `ExpandA` — it is generated **directly in the NTT domain** (`expand_A0_ntt`): the uniform `[0,q)` rejection samples ARE the forward-NTT coefficients (no `ntt16`), with `A0 := INTT(samples)`. Uniformity is preserved because NTT is a bijection $R_q \leftrightarrow T_q$ (incl. the incomplete-NTT product rings of modes 120/180). `b` is relational (`b = (A0 s0 + e) f^{-1}`), so it is computed normally and then **cached** as `b̂ = NTT(b)` once per keygen/sign/verify. The structural `qj`/`2I_m`/parity terms stay in the coefficient domain (`HighBits`/`LSB` are coefficientwise); only the convolution part `-b·x0 + A0·xs` runs in NTT (`docs/design-notes/PreNTT.md`).

This removes **all** forward NTTs of `A0` (it is born transformed) and, crucially, the **per-attempt** re-NTT of `A0`/`b` that the rejection loop previously paid each try (the transforms used to live inside `mat_mul_2q`; they now live once, before the loop). `mat_mul_2q` / `mat_mul_z1_2q` / `compute_b` take NTT-domain `A0h`/`bh` inputs.

Because `A0`'s values change (now `INTT(samples)`), the KAT was re-recorded (user-approved, §1). **ref vs avx2 layout:** the avx `ntt.S` uses the **same primitive root** as the scalar `ntt16` (`gen_zetas.py`'s smallest primitive root), so the two are **transpose-equivalent** — their NTT-domain outputs differ only by the avx shuffle network, a pure bit-permutation of the coefficient index. So, exactly like Kyber `gen_matrix` / Dilithium `ExpandA`, the avx2 `expand_A0_ntt` samples in standard order and lays each poly into the avx slot order with `nttunpack_avx` — the forward shuffle ladder (`shuf8/4/2(/1)`) with **no butterflies**, emitted into `ntt.S` and narrowing the `int32` samples to `uint16` on load. Transpose-equivalence is asserted per-mode by the `tools/gen_ntt_perm.py` guard (run by `make check-consts`: it derives the avx<->ref permutation `P` from both Python NTT models and checks `P` is a bit-permutation — a root mismatch would make it non-affine, which no shuffle ladder can realize). The ladder's runtime correctness is pinned by `tools/ntt/<qset>/test_ntt.c` (`nttunpack(ntt_ref(p)) == ntt(p)`) and end-to-end by the avx2-vs-ref KAT. The key bytes now store canonical `bn = NTT_ref(b)` in the scalar `ref/` order, so AVX2 imports `bn` with the same `nttunpack_avx` ladder before using it as native-layout `b̂`; this keeps pk/sk portable instead of serializing an AVX-private slot order.

## 7. Float ApproxExp (opt-in, `AVX2EXP=1`)

Integer `approx_exp` uses `mulh64` ($64\times 64\to 128$ high), which AVX2 cannot vectorize, so it stays scalar (ILP-batched, KAT-identical to ref) by default. The opt-in float path (`approx_exp_avx2.c`) uses a degree-6 binary64 minimax (B=8 range reduction). **Route (ii)**: the kernel emits `exp(-a)` as a double, converted to a Q63 uint64 so the accept test remains the exact integer compare `rand_rej_63 < exp_q63` in the shared scalar `gauss_finalize` — isolating the float approximation to a single value. The production wrapper is now fully vectorized around the W16 kernel: Q60 inputs are split into uint32 hi/lo halves and converted by `vcvtdq2pd` plus masked `$2^{32}$` repair, while Q63 outputs are recovered from the binary64 exponent/mantissa with AVX2 variable shifts; this removes the old scalar `uint64_t`/`double` casts, `vcomisd`/`jnb` unsigned-conversion branch, and `vcvttsd2si` tail from the secret-dependent path.

KAT: this is a **separate, opt-in KAT** (float $\ne$ Q63-integer rounding). In practice the $\sim 2^{-50}$ approximation error is far below the rejection granularity, so for the recorded vectors it matches the reference KAT exactly (see `KAT_avx2exp.txt`); this is **not guaranteed** in general (a $\sim 2^{-43}$-measure set of inputs could flip an accept decision).

**Microarchitecture sweep** (`speed_exp_avx2.c`, the user-requested exploration): W4 (1 stream) / W8 (2) / W16 (4 `__m256d` Horner chains). The production path is **W16**, called twice by the `GAUSS_BATCH=32` wrapper, filling the two FMA ports across the ~4-cycle FMA latency. The kernel was subsequently rebuilt in the §16 port-model round — a masked-OR-tree table lookup (replacing the 32-deep `vblendvpd` scan) plus a shorter B=8 deg-6 range reduction — taking it from $\approx 4.84$ to **1.06 ns/exp at W16 (4.55$\times$)**; see §16 for the full audit. Determinism: explicit FMA intrinsics + `-ffp-contract=off -fno-fast-math`; FTZ/DAZ set and restored around the kernel.

## 8. Build system

`make` (default `AVX2EXP=0`, integer exp, ref-KAT) builds all tests $\times$ 3 modes.  Flags: `-O3 -mavx2 -mbmi2 -mpopcnt -march=native`.  Symbol namespacing (`LITHIUM_NAMESPACE`, all renames centralized in `ref/namespace.h`) is on by default; build with `-DDISABLE_NAMESPACE=1` to keep the plain un-mangled names (single-library / debugging / indexing).  Targets:
- `make check` — build + run all tests (13 tests $\times$ 3 modes = 39 PASS).
- `make check-kat` — assert the three recorded FNV hashes.
- `make check-kat-avx2exp` — lock the float-exp build's KAT against `KAT_avx2exp.txt`.
- `make check-refill` — force the 4-way SHAKE continuation-squeeze path (`-DGAUSS4X_INIT=4 -DNOISE4X_INIT=4`) and assert the KAT is unchanged.
- `make check-consts` — reproducible-constants gate (Section 9), incl. the single-`KeccakF_RoundConstants`-symbol `nm` check.
- `make ct-scan` — assert no division/sqrt in secret-handling code, and additionally reject scalar FP conversion/compare opcodes in the opt-in AVX2Exp wrapper.
- `make AVX2EXP=1 …` — float exp variant (own KAT in `KAT_avx2exp.txt`).
- `make AVX2EXP=1 TESTS=speed_exp_avx2` — float-exp microarch sweep (W4/W8/W16) + accuracy gate + dudect (test/control).
- `make TESTS=speed_sign` — end-to-end keygen/sign/verify cycle benchmark.
- `make profile` / `make AVX2EXP=1 profile` — per-component time + randomness breakdown of crypto_sign (Section 14).
- `make TESTS="…"` — build a subset.

### Test inventory

| test | what it checks |
|---|---|
| `test_ntt`, `t_poly` | scalar NTT + poly_mul/poly_invert (AVX path via `t_poly`) |
| `test_ntt_avx` | AVX `poly_mul` == schoolbook + `a*a^{-1}==1` (3000+ polys) |
| `test_sigma2_diff` | AVX 8-way σ2 == scalar exact-borrow chain over 400k random + forced `v_mid=0` and `v2==Z2` corner inputs |
| `test_sampler` | sampler stats + `4x==1x` lock |
| `test_irs`, `test_mod2q`, `test_rans`, `test_keys` | scheme components |
| `test_sign` | sign/verify accept-5/5 reject-5/5 |
| `test_kat` | full deterministic KAT FNV hash (the gate) |
| `speed_sign` | end-to-end cycle benchmark (ref vs avx2) |
| `speed_exp_avx2` | float exp W4/W8/W16 sweep + accuracy + dudect (AVX2EXP only) |

## 9. Reproducible constants (`make check-consts`)

Per the project principle (no magic numbers; every constant reproducible + double-checked). Generators rewrite a delimited `@@AUTOGEN:<tag>@@` region in place (via `tools/autogen.py`), leaving the hand-written parts intact, so `make tables` regenerates and the constants still diff byte-for-byte:

| Constant | Region (file) | Generator | Check |
|---|---|---|---|
| RCDT `RCDT_Z`/`RCDT_NOISE` | `@@rcdt@@` in `sampler.c` (ref **and** avx2 copy) | `tools/gen_rcdt.py` | `cmp` both byte-stable |
| NTT/Montgomery (`MONT`/`zetas`/…) | `@@ntt@@` in `ntt.h` | `tools/gen_zetas.py` | (`make tables`) + `ntt_ref.py` model |
| norm bounds (`BK_SQ`/…) | `@@bounds@@` in `params.h` | `tools/gen_bounds.py` | `make tables` |
| exp coeffs/tables (`EXP_POLY`/…) | `@@approx_exp@@` in `approx_exp.h` | `tools/gen_approx_exp.py` | runtime acc gate (`t_exp`) |
| rANS tables / reserve / irs consts | `@@rans_tables@@`/`@@rans_meta@@`/`@@irs_consts@@` in `irs_rans.h` | `gen_rans_tables.py`/`SigSize.py`/`irs_ref.py` | `check_rans.py` |
| per-mode NTT `ntt_consts.c` + `ntt.S` | vendored in `<qset>/` | `tools/ntt/<qset>/gen_ntt.py` (bit-exact model) | regenerate + `cmp`; `test_ntt_avx`/`t_poly` |
| float exp `C[7]`, `T[8]`, derived scalars | `@@AUTOGEN:approx_exp_avx2@@` region in `approx_exp_avx2.h` | `tools/avxexp/gen_avxexp_consts.py` (Sollya log under `LogAVXApprox/`) | regenerate region + `cmp` + max-rel-error $< 2^{-43}$ |
| AVX2Exp KAT vectors | `KAT_avx2exp.txt` | `make AVX2EXP=1` capture | `check-kat-avx2exp` |

`avx2 make check-consts` regenerates the AVX2-side constants (float exp, in-tree NTT asm, the `@@rcdt@@` region in both sampler.c copies) and fails on any drift or accuracy miss; the remaining `@@AUTOGEN@@` regions are covered by `ref make tables` + their runtime golden-reference gates.

## 10. Constant-time stance

- **No division on secret data** (primary, authoritative): `make ct-scan` disassembles the secret-handling objects (`sampler`, `polyvec`, `approx_exp`, the float `approx_exp_avx2`, and all three vendored `ntt.S`) and asserts no `idiv`/`div*`/`vdiv*`/`vsqrt*`; for `approx_exp_avx2` it also rejects the scalar FP conversion/compare pattern (`cvtsi2sd`, `cvttsd2si`, `comisd`) that would reintroduce the GCC unsigned-cast branch — passes.
- **No secret-dependent branches/gathers**: the 8-way σ2 and the NTT assembly are straight-line SIMD; accept/zero-fold/sign are branchless masks (reused from ref); the float-exp table is a full masked-OR scan (not an indexed gather) and `2^{-n}` is formed by a constant-distance exponent-field shift (52) after computing the secret-derived but range-bounded `n`. FTZ/DAZ (saved/restored) removes the only data-dependent FP-latency hazard (denormals, which never arise here — exp $\in$ [6e-4, 1]). Float mul/add/fma/round on normals are data-independent per Intel DOIT / Arm DIT, which the kernel relies on by using only those instructions.
- **Empirical (supplementary)**: the reference `tools/avxexp/bench.c` dudect on the i7-11700K Rocket Lake benchmark system (Turbo Boost off, Hyper-Threading off, native Ubuntu) gives |t|$\approx$1.5 (AVX2) — no leakage. The in-tree `speed_exp_avx2` dudect prints a **test (fixed/random)** *and* a **control (random/random2)**: on the i7-11700K Rocket Lake no-Turbo/no-HT native-Ubuntu test environment the control |t| ($\approx$200–300) is far above the test |t| ($\approx$1–2) because of a memory-placement measurement bias, i.e. the harness noise floor exceeds the (zero) real signal — so the in-tree number is informational only and the static ct-scan + bench.c are authoritative. (A real data-dependence would show control$\approx$0, test high — the opposite of what we observe.)

## 11. Performance (2026-06-11 i7-11700K Rocket Lake, Turbo Boost off, Hyper-Threading off, native Ubuntu)

End-to-end median cycles from `test/speed.txt` (`make speed_txt`). The scalar reference uses gcc `-O3`; the AVX2 build uses `-O3 -march=native -mno-avx512f`. Parenthesized factors are speedups over `ref` from the same native-Ubuntu run.

| mode | op | ref (scalar) | avx2 | speedup |
|---|---|---:|---:|---:|
| 120 | keygen / sign / verify | 1029k / 972k / 171k | 304k / 337k / 65k | $3.38\times$ / $2.88\times$ / $2.61\times$ |
| 180 | keygen / sign / verify | 1463k / 1585k / 272k | 302k / 547k / 98k | $4.84\times$ / $2.90\times$ / $2.78\times$ |
| 260 | keygen / sign / verify | 2083k / 2045k / 307k | 620k / 689k / 118k | $3.36\times$ / $2.97\times$ / $2.61\times$ |

The updated profile (`test/profiling.txt`) shows the mode-120 sign path dominated by `sample_y` at 229.2k cyc/op (64.7%), followed by `sign_core` at 66.6k (18.8%) and rANS pack at 19.8k (5.6%). Inside `sample_y`, Keccak4x is 120.6k, integer ApproxExp is 56.2k, `a_q60+finalize+y` is 26.0k, and BaseSampler is 13.8k. Keygen is dominated by `sample_noise` 207.1k plus `compute_b` 49.3k; verify is spread across unpack 23.6k, setup 11.7k, challenge 10.2k, `A0` expansion 9.3k, commitment NTT 6.5k, and usehint 5.8k. Full tables are in [`PROFILING.md`](PROFILING.md).

**Float ApproxExp (`AVX2EXP=1`).** Opt-in, separate KAT, and only the Gaussian-sampler exp is intended to change the hot path. Median sign cycles improve from **337k -> 303k** (Lithium-120, 10.1% faster), **547k -> 494k** (Lithium-180, 9.7% faster), and **689k -> 639k** (Lithium-260, 7.3% faster). In the mode-120 profile, the `ApproxExp` child drops **56.2k -> 26.8k cyc/op** and `sample_y` drops **229.2k -> 198.8k cyc/op**. Keygen/verify medians in the float build are reported in `test/speed.txt` for completeness, but the exp variant is a sign-path optimization and those changes should be treated as run noise plus different benchmark ordering rather than an algorithmic win.

## 12. As-built deviations from the plan (all improvements)

- **8-way σ2 uses `RCDT_Z` directly** (vpbroadcastd) instead of a generated `RCDT_Z_AVX2` broadcast table — eliminates a drift-prone duplicate; the canonical RCDT table is the single source. (Plan had proposed generating a broadcast copy.)
- **σ2 borrow is the exact bit-31 chain** (`cmpgt, add, sub, srli, sub, sub, srli` per group/row), not a sign-flip unsigned `<` chain — limbs are <=31 bits so a signed subtract's sign bit is the borrow-out, halving the ALU ops and giving the exact 93-bit compare. Both ref and avx2 use it; KAT-neutral on the committed vectors (Section 5). Locked by `test_sigma2_diff`.
- **`poly_invert` inlines its scalar transforms** rather than calling a separate `poly_invert_scalar.c` — fewer files, same isolation guarantee.
- **Refill path verified** by `make check-refill` (forced tiny init), since the production init sizes usually avoid it.
- **Unified single-stream Gaussian `y` (§6.1)** — all K y-polys sampled from one continuous keccak4x combine (persistent buffer, leftover carries across polys); removes both padding-lane and lockstep-uneven-completion waste. Gaussian SHAKE squeeze $\approx$1.006$\times$ consumed (near-zero), sign ~11–19 % faster than the original per-poly-single-stream avx2. User-approved KAT re-record (ref==avx2).
- **GAUSS4X_INIT per mode** (56/68/104) — now only affects the test-only per-nonce `sample_gauss_N_4x`; the unified `y` stream refills on demand.
- **Vectorized sign_core + rounding (Tier 1, KAT-NEUTRAL — no re-record)** — exact integer SIMD guarded by `#ifdef __AVX2__` directly in `irs.c`/`sign.c`/`polyvec.c` (the `#else` is the portable scalar reference; ref built without `-mavx2` stays scalar). Both paths are bit-identical (integer add is associative), so the KAT is unchanged. Constant-time: only add/sub/mul (`vpmuldq`/`vpmullo`)/shift/and/or/cmp, no div/sqrt, no secret branches (loop bounds are public: the challenge index `j`, and `M`/`N`); verified by objdump on irs/sign/poly.
- `irs.c`: `inner` (dot product via exact `vpmuldq`), the masked `z$\pm$v` update, and `poly_cyclic_shift` (negacyclic rotate). sign_core $117k \to 55k$ (~2.1$\times$).
- `sign.c`: HighBits (buckets) and makehint loops — inputs are pre-bounded (`w`$\in$[0,2q), `|w-2z2|`<2q) so the scheme Barrett is a no-op $\to$ shifts/masks only. makehint $16k \to 0.7k$, highbits $11k \to 6.6k$.
- `poly.c`: `poly_sqnorm` (norm coeffs are `|·|≪q`, so freeze+center is masked centering, no Barrett; exact int64 sum via `vpmuldq`). normcheck $10k \to 0.6k$; also speeds keygen/verify norms.
- Net: sign ~22 % faster than the unified-scheme baseline (2.6–2.9$\times$ vs ref at this point; the branchless packing below then raised it to 3.2–3.6$\times$).
- **Branchless bit-packing (profiling-driven, KAT-NEUTRAL, portable — helps ref too).** The profiler flagged `challenge` (~10 %) and `pack` (~17 %) as scalar hot spots. Both were dominated by the per-bit `if (x & (1<<b)) buf[..] |= ..` packing idiom, whose data-dependent branch mispredicts ~50 % on random data. Replaced with a bit-accumulator writer/reader whose byte-emit schedule depends only on the loop index (no data branches, no `memset`), producing byte-identical output (`N·W1_BITS` and the pack widths are byte-aligned, the streams are continuous $\to$ exact):
- `sign.c hash_challenge`: challenge **$44.6k \to 7.7k$** cyc/sig (mode 120).
- `packing.c bitpack`/`bitunpack`: drives the z-lo / pk / sk / hint packing — `pack` **$71.7k \to 23.9k$**; **keygen ~2$\times$** (pk/sk pack was branch-bound) and **verify ~1.6$\times\to$** faster unpack. End-to-end: sign **+~8 %**, keygen **~2$\times$**. This is not SIMD but a portable scalar win; it lives in `ref/` (no `#ifdef`), so the scalar reference gets the same speedup. A genuinely-SIMD 8-wide `a_q60` (exact, since `SAMP_I_HI==0`) was prototyped and **reverted** — measured a wash (a_q60 runs ~100$\times$/sig and the compiler already pipelines the scalar mulq; its `…+finalize` bucket is dominated by the ~1640 serial finalize-compaction calls, which are not SIMD-friendly).
- **rANS encoder divide $\to$ multiply-by-reciprocal (KAT-NEUTRAL, portable).** `enc_put`'s `x/freq` + `x%freq` (one hardware divide per encoded symbol, ~1280/sig) is replaced by the ryg rANS reciprocal: `q = (x*RCP >> 32) >> RSH`, `x = x + BIAS + q*(PSCALE-freq)`. By Granlund–Montgomery (RCP = $\lceil 2^{shift+31}/freq\rceil$, shift = $\lceil\log_2 freq\rceil$, RSH = shift-1) this is `q == floor(x/freq)` EXACTLY for the encoder's range ($x < x_{max} \le 2^{31}$), so the rANS byte stream — hence the KAT — is unchanged. RCP/RSH/BIAS are generated into the `@@rans_tables@@` region (`tools/gen_rans_tables.py`) with a per-slot self-check (reciprocal state-update == plain-divide over a thorough x sample) and are byte-exactly validated end-to-end by `test_rans` (C vs `rans.py`). pack **$23.9k \to 16.9k$** cyc/sig (mode 120, avx2; ref $21.9k$). Shared `irs_rans.c`, so ref benefits too. (The decoder has no divide.)
- **§12.1 Post-Pre-NTT finalize/LSB cleanup (profiling-driven, KAT-NEUTRAL).** Two small wins the post-Pre-NTT profile exposed (avx2 sign $\approx2$–$4\%$, verify $\approx2$–$4\%$): **(1)** the challenge `LSB(w[0])` loop (`sign.c`) called `lsb()` (a full Barrett mod-2q) although `w[0]` is already in $[0,2q)$ from `mat_mul_2q`, so `lsb()==\&1` — replaced with the direct mask (portable, helps ref too: highbits $6.8k\to3.4k$ avx2, $12.1k\to8.7k$ ref). **(2)** the `mat_mul_2q`/`mat_mul_z1_2q` finalize (`avx2/polyvec.c`) — `freeze`d `t` gives `v = 2t + q\cdot\text{parity}` in $[0,3q)$ (sign) / $[-q,3q)$ (verify), so the per-coeff `reduce_mod_2q` Barrett collapses to a branchless conditional $\pm 2q$; vectorized `__m256i` under `#ifdef __AVX2__` (scalar `#else` is the ref). commitment finalize $16.9k\to13.6k$; speeds both sign and verify (verify uses `mat_mul_z1_2q`). Both are bit-identical (KAT unchanged) and constant-time (mask arithmetic, no secret branch).
- **AVX-512 dropped** from scope (AVX2 only), per decision. (Caveat: `-march=native` on an AVX-512 host lets the compiler emit AVX-512 for auto-vectorized C — confirmed for the *opt-in* float-exp kernel's `vcvtuqq2pd`; the default integer/KAT build stays pure AVX2.)

- **§12.2 NTT input conversion in assembly (profiling-driven, KAT-NEUTRAL).** The old `poly_ntt` boundary paid a scalar `freeze()` loop before entering the hand-written NTT, and `poly_ntt_small` still had a scalar mask-add prefix after the small-input split. The AVX2 qset `ntt.S` files now expose three standard-order conversion helpers: full `poly_freeze_to_u16_avx` (exact Barrett `freeze()` with the mode-180 33-bit reciprocal split), small-signed `poly_freeze_small_to_u16_avx`, and `[0,q)` `poly_cast_modq_avx`. `poly_ntt`/`poly_ntt_canonical` use the full helper, `poly_ntt_small` uses the small helper, and `compute_b` uses cast for `finv` plus small for the reduced `(A0s0+e)` accumulator. The helpers only convert representation; they do not apply the AVX slot shuffle, so canonical `bn = NTT_ref(b)` serialization remains ref-order and AVX imports still go through `nttunpack_avx`. Differential `test_freeze_avx` covers full int32/corner, small, and cast inputs for all modes; KAT is unchanged.
## 13. Provenance / regeneration

- NTT assembly + constants: `tools/ntt/<qset>/gen_ntt.py` (run `make gen`/`make test` there to regenerate and bit-exact-validate). Linked by explicit path.
- Float exp coefficients: `tools/avxexp/run_avxexp.sh` (Sollya fpminimax, interval `[7/8 ln2, ln2)`, target $2^{-43}$), logged under `tools/avxexp/LogAVXApprox/`; embedded + checked by `gen_avxexp_consts.py`.
- RCDT / other IRS tables: `tools/*.py` (`make tables` in `ref/`).

## 14. Profiling infrastructure

Two independent, macro-gated probe systems (`prof.h`/`prof.c`, shared via the ref symlinks).  Both are **no-ops unless their macro is defined**, so normal / KAT / CT builds are byte-identical and zero-cost — verified: `make check` is warning-free and KAT-unchanged, and a KAT build *with* `-DPROF_TIME -DPROF_RAND` reproduces the same hashes (probes never touch crypto state).

- **`-DPROF_TIME`** — per-component cycle accounting.  Disjoint regions are timed with `PROF_START/PROF_STOP` (rdtscp); top-level buckets sum to the sign total, sub-buckets are a breakdown of their parent.  Buckets: setup, expand_A0, sample_y {shake, BaseSampler(sigma2), ApproxExp, a_q60+finalize}, commitment {ntt fwd, pointwise/basemul, inv}, highbits, challenge, sign_core {approx_exp}, normcheck, makehint, pack.
- **`-DPROF_RAND`** — randomness accounting.  `PROF_SQ` in the SHAKE squeeze functions (fips202.c keccak_squeeze/squeezeblocks; fips202x4.c for the 4-lane path) tallies bytes squeezed against the current `PROF_CTX` (set around each sign region: setup/A0/gauss/challenge/irs); `PROF_USE` tallies fine Gaussian consumption (signs/sigma2/y/rej_tail).

Run: **`make profile`** (avx2/, integer exp) or **`make AVX2EXP=1 profile`** (float exp); **`make profile` in `ref/`** for the scalar baseline.  Runs N=2000 real signatures (incl. rejection retries) and prints per-sig cycles/% and SHAKE bytes/%.  Probe overhead is ~4–5 % of the sign total (lands in the hot inner buckets); use it for *relative* shares, not absolute cycle counts.

**Findings on the updated benchmark system (mode 120; i7-11700K Rocket Lake, Turbo Boost off, Hyper-Threading off, native Ubuntu):**

(profiling-on cycles; current source plus the O(1) rANS decode LUT, AVX2 key-noise CDT, and shared Montgomery batch inversion):

| component | ref (scalar) cyc | avx2 cyc | speedup |
|---|---:|---:|---:|
| setup | 7.1k | 5.9k | 1.2x |
| expand_A0 | 28.0k | 9.4k | 3.0x |
| sample_y (Gaussian) | 571.0k | 229.2k | 2.5x |
|   - shake (keccak) | 367.7k | 120.6k | 3.0x |
|   - BaseSampler | 77.6k | 13.8k | 5.6x |
|   - ApproxExp (int) | 67.6k | 56.2k | 1.2x; AVX2Exp profile is 26.8k |
|   - a_q60+finalize+y | 42.9k | 26.0k | 1.6x |
| commitment (NTT mat-mul) | 143.8k | 10.6k | 13.6x |
| sign_core | 191.0k | 66.6k | 2.9x |
| highbits / makehint / normcheck | 23.4k | 2.8k | 8.4x |
| challenge (hash) | 11.9k | 9.9k | 1.2x |
| pack (rANS + z-lo) | 20.3k | 19.8k | 1.0x |
| **TOTAL** (sum) | **996k** | **354k** | **2.8x** |

Post-optimization AVX2 sign is about 65% `sample_y` and 19% `sign_core`; commitment is down to about 3%. Mode-120 keygen is dominated by `sample_noise` 207.1k and `compute_b` 49.3k, while verify is spread across unpack 23.6k, setup 11.7k, challenge 10.2k, `A0` expansion 9.3k, commitment NTT 6.5k, and usehint 5.8k. The live per-mode tables are in [`PROFILING.md`](PROFILING.md).

Randomness (mode 120, bytes/sig): AVX2 sign squeezes about **39.3k** total SHAKE bytes, of which Gaussian sampling accounts for **36.2k**. Gaussian consumption is about **35.5k** (sigma2 58.3%, rej_tail 37.0%, y 4.3%, signs 0.5%). Going below this would require byte-granular scalar squeezing, which gives up keccak4x throughput and is a net loss.


## 16. Microarchitecture audit round (256-bit port model)

On-device-grounded 256-bit port model from uops.info (the recorded CPU generation), recorded in `tools/uarch/RESULTS_avx2.md`: **three** vector ports — FMA/mul/cmp/shift on **p0+p1 (2/cyc, lat 4)**, shuffles/permutes/widens on **p5 (1/cyc)**, int add/logic on **p0+p1+p5 (3/cyc)**.  (Contrast the 512-bit domain's *two* ports + 1 FMA/cyc.)  The build was first made a **pure AVX2 path** (`-mno-avx512f` + `make check-no-avx512`): `-march=native` had been letting GCC auto-vectorize the symlinked scalar C into AVX-512 (`vpermt2w` appeared in a `-march=native` binary), so the "AVX2" numbers were mildly contaminated; the flag makes it genuinely AVX2-only (perf-neutral).

A 7-component IPC audit (ideal port-floor vs measured) drove these:

1. **AVX2Exp — 4.55x (the flagged component).**  The `2^{-j/32}` lookup was a 32-entry *sequential* `vblendvpd` scan with a per-entry `vpmovsxdq` (p5): p5-bound **and** a 32-deep latency chain (the inefficiency AVX-512 killed with `vpermt2q`, which AVX2 lacks).  Fix in two steps, both bit-identical to the float KAT: (a) **masked-OR tree** — widen `j` once, OR-mask all entries into 4 accumulators (`vpcmpeqq`+`vpand`+`vpor`), dropping the 32 `vpmovsxdq` and the blend chain -> W16 4.84->2.30 ns/exp; (b) **B=32->8, deg-6** (Sollya refit on `[7/8 ln2, ln2)`) — AVX2 has no cheap 32-double permute, so the scan *length* is the lever -> W16 2.30->**1.06 ns/exp** (tsc/exp 18->4.1), near the FP floor.  In the updated mode-120 profile, the AVX2Exp bucket is 56.2k -> 26.8k cyc/sig; float-exp KAT unchanged.
2. **gauss a_q60 + sampler_y — vectorized.**  `a_q60 = (num*L1)<<12 + (num*L0)>>20` (exact split of the 64-bit `I_LO`, 4-wide `vpmuludq`); the 7-bit `sampler_y` unpack via `vbroadcasti128`+`vpshufb`+`vpsrlvd` (AVX2 has no `vpermb`/`vpsrlvw`).  Both bit-identical -> shared KAT unchanged. a_q60+finalize+y bucket **35.7k->22.6k** cyc/sig.

### Follow-up round (profiling-driven; both KAT-neutral)

3. **rANS decode `find_slot` -> O(1) `val->slot` LUT (the verify elephant).**  The decoder's per-symbol slot lookup ([`irs_rans.c`](../ref/irs_rans.c)) was an **O(n) linear CDF scan** (n = 121–165 z-head / 63–85 hint) over `(ℓ+1+M)·N ~ 1280` symbols/sig — **~73 % of verify**, scalar and shared (the `irs_rans.c` symlink), so identical in ref/avx2/avx512.  Replaced by a direct `SLOT[2^prob_bits=1024]` table (`val->slot` in one load, `uint8`, emitted by `gen_rans_tables.py`; decode result identical -> KAT-neutral, no wire change). **current mode-120 AVX2 `unpack` is 23.6k cyc/op and verify is 71.0k; the shared LUT keeps decode from returning to the old linear-scan floor** — and it lifts `ref/` and avx512 too.  (Encode already uses indexed start/freq, no scan.)
4. **avx2 key-noise CDT — vectorized (keygen).**  `noise_magnitude_batch` had stayed **scalar** in ref/avx2 (only avx512 had a SIMD version).  Added the 8-way kernel to `sampler.c` — the exact 31-bit-limb borrow chain already proven for `sampler_sigma2`, over the 10-entry `RCDT_NOISE` table and `NOISE_BATCH` samples, reading the identical grouped layout -> bit-identical, KAT unchanged.  **current mode-120 AVX2 `sample_noise` is 207.1k cyc/op and keygen median is 304k cycles in `test/speed.txt`.**  Residual: noise SHAKE + the per-candidate accept tail.
5. **`poly_invert` — Montgomery batch inversion (keygen).**  `poly_invert` inverted `f` in the residue domain with a **per-element Fermat** modexp (`a^{q-2}`): `N` elements for the full NTT (mode 260, 512x) or `N/2` pair-determinants for the incomplete NTT (modes 120/180, 128x) — the bulk of `compute_b`, shared-scalar in all three builds.  Replaced by Montgomery's trick (`batch_modq_inv` in shared `polyvec.c`): **one** Fermat + ~3n modmuls for the whole vector.  Each result is the unique modular inverse -> bit-identical, KAT-neutral; ref/avx2/avx512 all gain.  **current mode-120 `compute_b` is 49.3k (avx2) and 53.9k (avx512)**; saving scales with the inversion count (largest at mode 260).  Residual: the 2 scalar NTTs in `poly_invert` (they need the scalar Cooley-Tukey coeff order).

**Audit verdicts (does it achieve ideal IPC?):** NTT — **yes** (at its p015 modular-reduction floor; q>2^15 allows no lazy headroom).  shake/keccak4x — **no, 4.8x off, latency/ILP bound**, and ~2x of that is a fundamental **AVX2 ISA wall**: no `vprolq` (rho = 3 ops) and no `vpternlog` (chi = 2 ops), the exact ops the avx512 keccak8x fuses to 1 each — uncloseable without AVX-512.  finalize compaction — the `r[ctr++]` left-pack has **no AVX2 `vpcompress`**, so it stays scalar.  Smaller items (sigma2 mask-borrow, `inner` `vpmuldq`, rANS prep-loop) are ~1% each, below the i7-11700K Rocket Lake benchmark system (Turbo Boost off, Hyper-Threading off, native Ubuntu)'s ~7% `speed_sign` run-to-run noise; analyzed, deferred.  Full audit: the workflow report (per-component port tallies + ranked plan).

Full per-component cycle profiles (keygen / sign / verify, sorted by time-share): [PROFILING.md](PROFILING.md).
