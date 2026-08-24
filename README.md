# Lithium-Code

A from-scratch C reference implementation of the **IRS** post-quantum lattice signature scheme — Gärtner, *Compact Lattice Signatures via Iterative Rejection Sampling* (the **compressed** variant `CSign`/`CVerify`, the paper appendix section "A More Compact Signatures") — for three NIST security levels.

The scheme is built on the **NTWE** problem (not plain LWE): the public key is `b = (A0·s0 + e)·f^{-1} mod q` with secret `s = [f, s0, e]`, and the verification matrix `A = [qj − 2b, 2A0, 2I] mod 2q` is kept implicit (reconstructed from `b` via the mod-2q lift). Parameters are the target-B no-abort set documented in `docs/design-notes/param.md`, not the example parameters in the paper.

> Current snapshot: the scalar, AVX2, and AVX-512 builds share the same default KATs across all three modes, and the retained local audit reports under `test/` record passing constant-time scans, security regression tests, and integration checks on 2026-06-06.

## Quick Artifact Evaluation

Use the scalar reference as the fastest portability smoke test, then run the full in-tree audit gate on an AVX-512-capable x86-64 machine when available. The AVX2 backend needs AVX2, BMI2, and POPCNT; the AVX-512 backend additionally needs AVX-512F/BW/VL/DQ/IFMA. If the host lacks AVX-512, run the scalar and AVX2 commands and treat the AVX-512 rows as architecture-specific evidence; the upstream Wikipedia [CPUs with AVX-512](https://en.wikipedia.org/wiki/AVX-512#CPUs_with_AVX-512) list is a convenient starting point for finding compatible machines.

```bash
make -C ref check                  # portable scalar build and functional/KAT gates
./run_tests.sh                     # full ref + avx2 + avx512 verification gate
./reproduce.sh                     # optional third-party benchmark reproduction; network + longer runtime
```

The full `./run_tests.sh` gate writes retained evidence under `test/`; optional timing evidence is available through `./run_dudect_components.sh` and `sh tools/timecop/run_all.sh`.

## Performance

Median `speed_sign` cycles were re-recorded on 2026-06-11 on an i7-11700K (Rocket Lake, 2021) running native Ubuntu with Turbo Boost disabled and Hyper-Threading disabled. The scalar build uses gcc `-O3`; AVX2 and AVX-512 builds additionally use `-march=native` (AVX2 also forces `-mno-avx512f`). The default AVX2 and AVX-512 builds are byte-identical to the reference KAT. Parenthesized factors are speedups over `ref`; `keygen` includes a rejection loop and is noisier than `sign`/`verify`.

**sign** -- the hot path:

| security | ref | avx2 | avx2 + AVX2Exp | avx512 | avx512 + IFMA exp |
|----------|-----|------|----------------|--------|-------------------|
| Lithium-120 | 972k | 337k ($2.9\times$) | 303k ($3.2\times$) | 211k ($4.6\times$) | 166k ($\mathbf{5.9\times}$) |
| Lithium-180 | 1585k | 547k ($2.9\times$) | 494k ($3.2\times$) | 344k ($4.6\times$) | 284k ($\mathbf{5.6\times}$) |
| Lithium-260 | 2045k | 689k ($3.0\times$) | 639k ($3.2\times$) | 442k ($4.6\times$) | 361k ($\mathbf{5.7\times}$) |

The opt-in exp variants (`avx2 + AVX2Exp`, `avx512 + IFMA exp`) carry separate KATs and accelerate only the Gaussian sampler's exp, so they are relevant to `sign`. The 2026-06-11 Rocket Lake measurements show AVX2Exp reducing median sign cycles by about 7--10% and IFMA exp reducing AVX-512 sign cycles by about 18--21%; each build writes the machine-stamped snapshot to its `test/speed.txt` ([ref](ref/test/speed.txt), [avx2](avx2/test/speed.txt), [avx512](avx512/test/speed.txt)).

**keygen**:

| security | ref | avx2 | avx512 |
|----------|-----|------|--------|
| Lithium-120 | 1029k | 304k ($3.4\times$) | 192k ($\mathbf{5.4\times}$) |
| Lithium-180 | 1463k | 302k ($4.8\times$) | 214k ($\mathbf{6.8\times}$) |
| Lithium-260 | 2083k | 620k ($3.4\times$) | 403k ($\mathbf{5.2\times}$) |

**verify**:

| security | ref | avx2 | avx512 |
|----------|-----|------|--------|
| Lithium-120 | 171k | 65k ($2.6\times$) | 57k ($\mathbf{3.0\times}$) |
| Lithium-180 | 272k | 98k ($2.8\times$) | 83k ($\mathbf{3.3\times}$) |
| Lithium-260 | 307k | 118k ($2.6\times$) | 103k ($\mathbf{3.0\times}$) |

With the O(1) rANS decode LUT in place, verify is no longer dominated by the old scalar CDF scan; it is spread across rANS unpack, `A0`/public-matrix NTT work, setup/usehint, challenge, and norm reconstruction. In the updated no-Turbo/no-HT native-Ubuntu run, AVX-512 verify is ahead of AVX2 by about 1.14--1.18x, while AVX-512 sign is about 1.56--1.60x faster than default AVX2 and about 1.91--2.03x faster with IFMA exp.


**Third-party comparison snapshot (BLISS excluded)**: optional `./reproduce.sh` timing evidence was re-recorded on the same 2026-06-11 Rocket Lake host from the upstream revisions listed in `third_party/`. The table reports median cycles and uses each upstream implementation's own harness/default flags; sizes are `pk+sig` for the same 59-byte benchmark message, so these rows are comparison context rather than a uniform security/format claim.

| scheme / implementation | pk+sig | keygen | sign | verify |
|---|---:|---:|---:|---:|
| Lithium-120 avx512 + IFMA exp | 2243B | 192k | 166k | 57k |
| Lithium-180 avx512 + IFMA exp | 3090B | 271k | 284k | 83k |
| Lithium-260 avx512 + IFMA exp | 4428B | 406k | 361k | 103k |
| Dilithium2 Xu avx512 | 3732B | 39k | 115k | 49k |
| Dilithium3 Xu avx512 | 5245B | 65k | 185k | 49k |
| Dilithium5 Xu avx512 | 7187B | 91k | 212k | 76k |
| Falcon-512 opt avx512 | 1563B | 11625k | 645k | n/a |
| Falcon-1024 opt avx512 | 3073B | 43339k | 1312k | n/a |
| HAWK-512 avx2 | 1579B | 5097k | 94k | 122k |
| HAWK-1024 avx2 | 3661B | 29058k | 189k | 253k |
| HAETAE2 avx2 | 2466B | 413k | 1163k | 92k |
| HAETAE3 avx2 | 3821B | 623k | 1537k | 147k |
| HAETAE5 avx2 | 5028B | 666k | 1804k | 173k |
| NTRU+SIGN512 ref | 1519B | 621k | 1332k | 71k |
| NTRU+SIGN1024 ref | 3215B | 1301k | 3015k | 134k |
| Raccoon-128-1 ref | 13780B | 3094k | 5501k | 1863k |
| Raccoon-192-1 ref | 17704B | 4697k | 7924k | 2978k |
| Raccoon-256-1 ref | 24394B | 7325k | 11715k | 5028k |

**AVX2 per-component profile** (mode 120, `make profile_txt`, default integer exp): `sample_y` is 229.2k cyc/op (64.7% of sign), with Keccak4x 120.6k, integer ApproxExp 56.2k, and `a_q60+finalize+y` 26.0k; `sign_core` is 66.6k; rANS pack is 19.8k. Keygen is dominated by `sample_noise` 207.1k and `compute_b` 49.3k. Verify is 71.0k total, with unpack 23.6k, setup 11.7k, challenge 10.2k, `A0` expansion 9.3k, commitment NTT 6.5k, and usehint 5.8k; see [`avx2/PROFILING.md`](avx2/PROFILING.md).

**AVX-512 per-component profile** (mode 120, `make profile_txt`, default integer exp): `sample_y` is 130.2k cyc/op (55.9% of sign), with scalar integer ApproxExp 53.4k, Keccak8x 45.8k, and `a_q60+finalize+y` 9.9k; `sign_core` is 53.0k; rANS pack is 19.8k. Keygen is dominated by `sample_noise` 136.5k, `pack pk/sk` 55.2k, and `compute_b` 53.9k. Verify is 63.0k total, with unpack 23.6k, setup 11.8k, challenge 10.0k, usehint 5.7k, commitment NTT 4.3k, and `A0` expansion 3.9k; see [`avx512/PROFILING.md`](avx512/PROFILING.md).

**Optional float / IFMA ApproxExp**: `make AVX2EXP=1` swaps the Gaussian-sampler exp for a 4-wide binary64 kernel, while `make AVX512IFMAEXP=1` uses the 8-lane AVX-512 IFMA integer exp with the same relative-accuracy target. In the updated profiles, AVX2 mode-120 `ApproxExp` drops 56.2k -> 26.8k cyc/op and `sample_y` drops 229.2k -> 198.8k; AVX-512 mode-120 `ApproxExp` drops 53.4k -> 11.0k and `sample_y` drops 130.2k -> 87.5k. Full per-component cycle profiles live in [`ref/PROFILING.md`](ref/PROFILING.md), [`avx2/PROFILING.md`](avx2/PROFILING.md), and [`avx512/PROFILING.md`](avx512/PROFILING.md).


## Implementation lessons from KyberSlash and ML-DSA bugs

Two recent implementation studies are now part of the engineering checklist for Lithium-Code: Bernstein et al. (2025) show that source code that looks constant-time can become variable-time when a compiler emits hardware division for secret-derived operands, and Bernstein (2026) shows that small, interoperable signature-code changes can pass typical tests while enabling rapid key recovery or universal forgery. The practical lesson is stricter than "keep KATs green": every optimization must preserve the exact distribution, transcript binding, parsing rejection behavior, and machine-code constant-time envelope.

For KyberSlash-style risks, the rule is **no runtime division, modulus, square root, gather/scatter, or secret-indexed table load on secret-derived data**. This includes apparently harmless divisions by compile-time constants: `-Os`, embedded compiler backends, and non-x86 targets can choose variable-time division instructions where `-O3` on the developer machine does not. In this tree, secret reductions use explicit reciprocal or power-of-two code (`reduce.c`, `rounding.c`, `irs.c`), secret table lookups use full scans (`approx_exp.c`), and SIMD CT gates disassemble objects with `make ct-scan`. When adding a reduction or bit-packing shortcut, inspect both source and generated assembly under at least `-O3` and `-Os`; if a divide is needed only for public data, document why the operand is public at the call site.

For ML-DSA-bug-style risks, assume that **small changes that interoperate are still dangerous**. Do not accept an optimization merely because signatures verify, KAT hashes match on one seed, or the change is local. Samplers, challenge generation, rejection decisions, norm checks, rANS parsing, zero padding, and transcript hashes must have negative tests and differential tests that exercise boundary cases, malformed inputs, and cross-backend equivalence. A performance patch that changes nonce scheduling, rejection-loop consumption, output-indexed signs, `z` splitting, hint ranges, or parsing failure order needs a reviewer to state the preserved invariant explicitly. The only retry/timing variations intentionally outside the risk list are the deployment-classified keygen loop, the near-$2^{-100}$ signing/rANS retry events fixed by our parameters, and Gaussian sampler variance covered by isochronous sampling.

The current audit of Lithium-Code found and tightened several paper-inspired hazards. `highbits()` now uses an explicit shift for `LITHIUM_TAU`, so `ref/rounding.c` does not rely on the compiler turning `/ LITHIUM_TAU` into a shift. Verification reconstructs `z2` through `reduce2q(diff)` instead of `% LITHIUM_DQ`, avoiding a compiler-selected divide in another hot arithmetic path. Signing packs `w0bits` branchlessly, scalar `poly_sqnorm()` centers coefficients without a secret-dependent branch, sampler sign-bit indexing uses shifts and masks rather than `/ 8` and `% 8`, and the AVX-512 Gaussian finalizer zeroes its sign-bit padding before wide loads. These fixes are deliberately small and KAT-neutral.

Performance work should keep the following checklist close to the code review:

- **Assembly is evidence.** For every secret-handling object, run a forbidden-instruction scan on the built object, not just a source grep; include optional exp variants and size-optimized builds when the source contains `/`, `%`, floating point, table lookups, or compiler-sensitive comparisons.
- **Distribution comes before cycles.** Any sampler, IRS, rANS, or challenge optimization must be checked against the scalar reference on random and crafted edge cases, including refill boundaries, zero folding, tail bits, noncanonical encodings, and overflow paths.
- **Retry timing is classified.** Keygen rejection is a one-time provisioning cost in common deployments such as certificate issuance and TLS 1.3 server authentication, where the resulting key is long-lived; signing retries from the norm gate and rANS overflow are parameterized to occur with probability near $2^{-100}$; Gaussian sampler runtime variance follows the isochronous Gaussian-sampling discipline and does not expose attacker-usable information. Benchmarks should still report tail ratios so accidental engineering regressions are visible.
- **Generated constants are code.** Keep constants in `@@AUTOGEN@@` regions reproducible, diff them after regeneration, and require the generator to enforce accuracy or probability bounds instead of trusting copied tables.
- **Backends must be byte-for-byte comparable.** Scalar, AVX2, AVX-512, and optional exp backends need shared KATs or explicitly separate KATs; if a fast path is not byte-identical, the README/design note must say which transcript changes and why it remains valid.
- **Tests must fail on the bug class.** KATs catch drift, but they do not by themselves catch side channels, parser leniency, or rare boundary failures; pair them with negative parsing tests, sanitizer/fuzz targets where possible, and static CT scans.


## Verification Gates

Lithium-Code uses named verification gates instead of priority labels. The root `./run_tests.sh` command runs the constant-time scan, the security audit matrix, and the integration audit matrix in order; the retained reports live under `test/` so reviewers can inspect the last local evidence without rerunning every check.

### Machine-code constant-time scan

```bash
./ct_scan_matrix.sh
```

`./ct_scan_matrix.sh` writes `test/ct_scan_matrix.txt` and covers `ref`, `avx2`, and `avx512`; modes `120 180 260`; `gcc` and `clang`; `-O3` and `-Os`; plus the optional `AVX2EXP=1` and `AVX512IFMAEXP=1` exp variants. Backend-local quick checks remain available as `make -C ref ct-scan`, `make -C avx2 ct-scan`, and `make -C avx512 ct-scan`, but the root matrix is the command that proves full coverage. The scanner is intentionally strict on secret-handling code, with a narrow documented allowlist for public/fixed arithmetic: SIMD SHAKE gather instructions in `fips202x4.c` and `fips202x8.c` are limited to public lane pointer/offset handling, and GCC `-Os` divisions in `fips202.c`, `fips202x4.c`, `fips202x8.c`, and `sampler.c` are limited to public one-shot SHAKE output-length block counts or public Gaussian initial-refill block counts. The last local run in `test/ct_scan_matrix.txt` passed on 2026-06-06 with GCC 11.4.0 and Clang 14.0.0 and records these entries as `allowlisted_gather` or `allowlisted_public_div`.

TIMECOP/Valgrind support is optional dynamic-taint evidence that cross-checks the static object-code CT scan, not a mandatory gate. One command from the repository root runs the whole workflow:

```bash
sh tools/timecop/run_all.sh
```

`run_all.sh` builds a KyberSlash-style patched Valgrind on first use (under `external/timecop/valgrind-varlat/`, which is gitignored), proves that a tiny secret-division C program is detected (and a public division is not), then runs the signing smoke for the `ref` and `avx2` backends across modes `120 180 260` with signing key material and signing randomness marked secret. The gate is the KyberSlash property only — *fail iff a secret operand reaches a variable-latency instruction* (`--variable-latency-errors=yes`); secret-dependent control flow stays owned by `./ct_scan_matrix.sh` and `dudect-sign`, and the intentional constant-time `cmov` idioms those scanners certify are suppressed via `tools/timecop/lithium-ct.supp`. The `avx512` backend is **SKIPPED**: Valgrind 3.23 cannot decode AVX-512 (it aborts the guest with `SIGILL`), so the avx512 variable-latency property is enforced statically by `ct_scan` + `dudect-sign` instead (`TIMECOP_TRY_AVX512=1` reproduces the SIGILL). The methodology, host details, backend coverage, and last recorded result (2026-06-09: positive control + `ref`/`avx2` x all three modes PASS, variable-latency operands = 0) live in `docs/design-notes/TIMECOP.md`. The KyberSlash `valgrind-varlat` patch targets Valgrind 3.23.0; `build_valgrind_varlat.sh` downloads, filters, patches, and builds it automatically, so no system Valgrind is required.

### Security audit matrix

```bash
./security_audit_matrix.sh
```

`./security_audit_matrix.sh` writes `test/security_audit_matrix.txt` and delegates to three retained evidence reports: `test/fault_injection_matrix.txt`, `test/parser_negative_matrix.txt`, and `test/trace_kat_matrix.txt`. The root `./run_tests.sh` includes this script after `./ct_scan_matrix.sh`, so the default regression gate proves both the machine-code constant-time envelope and the parser/fault/trace coverage. TIMECOP remains optional and is intentionally not a default hard gate; it is available opt-in as `TIMECOP=1 ./run_tests.sh`, which appends `tools/timecop/run_all.sh` after the standard gates.

Single-gate commands are available when debugging a failure:

```bash
make -C ref MODES="120 180 260" TESTS=security_fault_injection check
make -C avx2 MODES="120 180 260" TESTS=parser_negative check
make -C avx512 MODES="120 180 260" TESTS=trace_kat check
```

The fault-injection gate covers adjacent coefficient repetition, repeated nonce, cleared signing secret, wrong hash input length, and wrong rejection counter behavior. A PASS means the detector reports `caught` for each injected bug class; a broken build or silently accepted faulty transcript is not counted as success.

The parser-negative gate signs a deterministic valid message, mutates the signature into signature-length edges, rANS reserve overflow, nonzero padding, noncanonical rANS stream, out-of-range hint/rANS corruption, out-of-bound `z`, and malformed challenge packing, and requires `ref`, `avx2`, and `avx512` to reject every case.

The trace KAT records a deterministic transcript hash plus randombytes call/byte counts, sign attempt count, and rANS overflow retry count. The matrix requires mode `120`, `180`, and `260` trace hashes to match across `ref`, `avx2`, and `avx512`; the last local run in `test/trace_kat_matrix.txt` passed with hashes `6853015861215245044`, `9614476440887310438`, and `15007787863907802895` respectively.

The secret/public classification baseline is `docs/design-notes/SECRET_PUBLIC_AUDIT.md`. Future public-data division/modulo allowlists, lookup/branch/gather allowlists, parser early exits, and declassification claims must cite that note. TIMECOP status and the decision to keep it external live in `docs/design-notes/TIMECOP.md`.

Last local security-audit result: `test/security_audit_matrix.txt` reports `OVERALL: PASS` on 2026-06-06 for backends `ref avx2 avx512` and modes `120 180 260`.

### Integration audit matrix

```bash
./integration_audit_matrix.sh
```

`./integration_audit_matrix.sh` writes `test/integration_audit_matrix.txt`, `test/diff_fuzz_matrix.txt`, and `test/perf_review_checklist_status.txt`. It builds and runs `diff_fuzz` for `ref`, `avx2`, and `avx512` across modes `120 180 260` with default `DIFF_FUZZ_N=32`, then compares normalized `CASE` lines and `transcript_hash` lines across backends. The root `./run_tests.sh` includes this script after `./security_audit_matrix.sh`, so the default regression gate proves constant-time coverage, parser/fault/trace coverage, differential fuzzing, and checklist coverage.

Single-gate commands are available when debugging a failure:

```bash
make -C ref MODES="120 180 260" TESTS=diff_fuzz check
make -C avx2 MODES="120 180 260" TESTS=diff_fuzz check
make -C avx512 MODES="120 180 260" TESTS=diff_fuzz check
```

The fuzz harness is deterministic and bounded, not a replacement for long-running fuzzing. Each case generates a keypair, a message, a signature, a tampered-message rejection case, and a malformed-signature rejection case; it records signature length, accept/reject bitmap, transcript hash, parsed rANS length, zero-padding status, challenge weight, `unpack_sig` return code, and successful `c`/`h`/`z1` hashes. The compact unpack-state summary uses the public `unpack_sig` helper and the signature bytes themselves; it does not add a public API.

Performance-sensitive patches must use `docs/design-notes/PERF_PATCH_REVIEW.md`; `.github/pull_request_template.md` asks authors to summarize distribution changes, randomness consumption changes, parser rejection-condition changes, machine-code CT-envelope changes, byte-for-byte KAT status, and new public-data-only allowlists. `./integration_audit_matrix.sh` checks that README links the checklist and records the result in `test/perf_review_checklist_status.txt`.

Component-level dudect is optional evidence, not a hard gate, because local timing noise and microarchitecture effects can produce WARN/FAIL observations that need interpretation rather than automatic blocking:

```bash
./run_dudect_components.sh
```

`./run_dudect_components.sh` writes `test/dudect_components_matrix.txt` and runs `dudect_components` for `poly_sqnorm`, `poly_ntt`/`poly_invntt_tomont` roundtrip, and a rounding/highbits/lowbits/lsb batch. It deliberately excludes Gaussian sampling, IRS full loops, rANS overflow retry, keygen retry, and the full signing retry loop. The last local optional component run on 2026-06-06 is recorded as `WARN` on this host; interpret it as component smoke evidence only. This is separate from the native-Ubuntu end-to-end signing dudect run retained in `test/dudect_all.txt`, which passed for every backend and mode with `DUDECT_N=200000`.

Last local integration-audit result: `test/integration_audit_matrix.txt` reports `OVERALL: PASS` on 2026-06-06 for backends `ref avx2 avx512`, modes `120 180 260`, and `DIFF_FUZZ_N=32`.

## `ref/` is a pure scalar reference

**`ref/` is a guaranteed pure scalar reference implementation: it contains NO SIMD/AVX2 code and no AVX2-awareness.** There are no `#ifdef __AVX2__` blocks and no `LITHIUM_USE_AVX2_*` routing macros in `ref/` — every `ref/` source compiles to portable scalar C and behaves identically with or without `-mavx2`. All vectorization lives in [`avx2/`](avx2/).

The `avx2/` build shares `ref/` as the single source of truth via **symlinks**, and maintains its **own real copies** only of the files that genuinely differ — `irs.c`, `sign.c`, `sampler.c`, `polyvec.c` plus `fips202x4.c`, `f1600x4.S`, `approx_exp_avx2.c`, and the **vendored per-mode NTT** in `<qset>/` (`ntt.S` + `ntt_consts.c`, copied from `tools/ntt` so the build is self-contained). `sampler.c` and `polyvec.c` carry all the AVX2 fast paths directly (consolidated — `sampler.c`: the 8-way $\sigma_2$/noise CDT + keccak4x 4-way backends; `polyvec.c`: the SIMD-NTT routing + shake128x4 uniform), gated by `LITHIUM_USE_AVX2_*`; the float-exp constants live in an `@@AUTOGEN` region of `approx_exp_avx2.h`. Everything else is a symlink back to `ref/`. The one opt-in exception that touches shared code is the float ApproxExp dispatch (`#ifdef LITHIUM_AVX2_EXP` in the shared `approx_exp.c`), which is dormant unless you build `make AVX2EXP=1`.

## AVX2 implementation

An AVX2-optimized implementation lives in [`avx2/`](avx2/) — AVX2 NTT, 8-way BaseSampler, 4-way SHAKE, vectorized `a_q60`/`y`-unpack, and an opt-in float ApproxExp (see Performance above). The default build is **byte-identical to the `ref/` KAT** (same three FNV hashes) and constant-time (no division/sqrt on secret data). It is a **pure AVX2 path**: `-mno-avx512f` forbids any AVX-512 that `-march=native` may otherwise enable on AVX-512-capable systems, enforced by `make check-no-avx512`. See [`avx2/DESIGN.md`](avx2/DESIGN.md) for the full design, KAT/CT/reproducible-constant policy, the §16 port-model audit, and benchmarks.

```bash
cd avx2 && make check        # 39 PASS (13 tests × 3 modes), KAT == ref
make check-kat               # assert the three recorded hashes
make check-consts            # every constant regenerates byte-for-byte + accuracy gate
make ct-scan                 # no div/sqrt in secret-handling code
make check-refill            # force the 4-way SHAKE refill path; KAT unchanged
make TESTS=speed_sign        # end-to-end cycle benchmark
make speed_txt               # -> test/speed.txt   (machine info + speed, both exp variants)
make profile                 # per-component time + randomness breakdown of sign
make profile_txt             # -> test/profiling.txt (machine info + profile, both exp variants)
make dudect-sign             # end-to-end fixed-secret vs random-secret signing t-test
make AVX2EXP=1 check         # opt-in float exp variant (own KAT)
```

`make profile` (here or in `ref/` for the scalar baseline) gives a per-component cycle + SHAKE-byte breakdown of `crypto_sign` (sample_y{shake/BaseSampler/ ApproxExp}, NTT, sign_core, pack, …) via zero-overhead macro-gated probes. See [`avx2/DESIGN.md`](avx2/DESIGN.md) §14.

## AVX-512 implementation

An AVX-512 implementation lives in [`avx512/`](avx512/), built the same way (symlinks to `ref/` + its own real copies, with the SIMD NTT/sampler consolidated into `polyvec.c`/`sampler.c` and the per-mode NTT vendored in `<qset>/`): **8-way Keccak** (`shake128/256x8`, XKCP-derived `keccakf1600x8`), the **AVX-512 NTT** (`ntt_avx512.S`, vendored from `tools/ntt`), a **16-lane** $\sigma_2$/noise CDT, `vpcompress`-based Gaussian finalize + uniform-`A0` rejection, `vpermb` 7-bit `y`-unpack, the 2-way interleaved rANS, and an opt-in 8-lane **IFMA integer ApproxExp**.  The default build is **byte-identical to the `ref/` KAT** (same three FNV hashes) and constant-time.  The design notes report measurements from one i7-11700K Rocket Lake native-Ubuntu test system with Turbo Boost and Hyper-Threading disabled as reproducibility context; see [`avx512/DESIGN.md`](avx512/DESIGN.md), where §16 covers the observed port model and the optimization campaign.

```bash
cd avx512 && make check       # 42 PASS (14 tests × 3 modes), KAT == ref
make check-kat                # assert the three recorded hashes (== ref == avx2)
make check-kat-avx512ifma     # lock the opt-in IFMA-exp KAT
make check-consts             # every constant regenerates byte-for-byte + accuracy gate
make ct-scan                  # no div/sqrt/gather in secret-handling code
make check-refill             # force the 8-way SHAKE refill path; KAT unchanged
make TESTS=speed_sign         # end-to-end cycle benchmark
make speed_txt                # -> test/speed.txt   (machine info + speed, both exp variants)
make profile                  # per-component time + randomness breakdown of sign
make profile_txt              # -> test/profiling.txt (machine info + profile, both exp variants)
make dudect-sign              # end-to-end fixed-secret vs random-secret signing t-test
make AVX512IFMAEXP=1 check     # opt-in IFMA integer exp variant (own KAT)
```

## Security levels

| mode  | n   | q     | ($\ell$, m) | κ  | τ    | NTT          | pk   | sk   | sig  |
|-------|-----|-------|--------|----|------|--------------|------|------|------|
| Lithium-120 | 256 | 57601 | (2, 2) | 30 | 512  | incomplete (deg-2) | 1056 | 1984 | 1187 |
| Lithium-180 | 256 | 30977 | (3, 3) | 58 | 512  | incomplete (deg-2) | 1472 | 2720 | 1618 |
| Lithium-260 | 512 | 64513 | (2, 2) | 58 | 1024 | full         | 2080 | 3808 | 2348 |

Sizes are in bytes; `pk` matches `param.md` exactly. The signature decomposes as `challenge (⌈n/8⌉) + rANS block (2 + R) + z_lo bit-pack`, with the rANS reserve `R = 481/688/746` sized for a strict $2^{-100}$ overflow probability (`tools/SigSize.py`).

## Layout

```
Lithium-Code/
├── ref/          C99 PURE-SCALAR reference implementation (per-level namespaces;
│   │             no SIMD/AVX2 — see "`ref/` is a pure scalar reference")
│   ├── *.c *.h   primitives + scheme (see "Source map" below)
│   │             Generated constants live in @@AUTOGEN@@ regions inside the
│   │             hand-written headers (params.h, approx_exp.h, ntt.h, irs_rans.h)
│   │             and sampler.c — regenerated in place by `make tables`; never
│   │             hand-edit a region between its @@AUTOGEN markers.
│   ├── namespace.h  central symbol renames (toggle with -DDISABLE_NAMESPACE=1)
│   ├── test/     unit tests (one gate per module) + prof.{c,h}
│   └── vector/   deterministic KAT generator + recorded hashes (KAT.txt)
├── avx2/         AVX2 build: symlinks to ref/ + real copies of the vectorized
│                 files (see "`ref/` is a pure scalar reference")
├── tools/        Python generators + bit-exact golden references
├── run_tests.sh  fail-fast end-to-end wrapper around the Makefile and audit gates
├── security_audit_matrix.sh     fault-injection, parser-negative, and trace-KAT audit
└── integration_audit_matrix.sh  differential-fuzz and performance-review checklist audit
```

## Build & test

The core scalar build requires a C99 compiler, GNU make, Python 3.8+ with only the standard library, and the usual binutils tools used by the Makefiles. The full verification gate additionally expects `gcc`, `clang`, `objdump`, `git`, and an x86-64 host with the SIMD features described in Quick Artifact Evaluation; unsupported SIMD backends should be skipped rather than emulated. Optional TIMECOP evidence downloads and builds Valgrind and therefore also needs `curl`, `tar`, a POSIX shell toolchain, and network access. Optional third-party benchmark reproduction uses `git` and network access to fetch the upstream projects listed under `third_party/`.

```bash
# Regenerate every generated header from tools/, then build and run the ref gates.
cd Lithium-Code/ref
make tables                  # regenerate every generated header from tools/
make check                   # build + run all gates, all three modes (default)
make MODES=120 check         # build + run the gates for a single mode
make check-rans-failrate     # statistical smoke test for rANS pack overflow
make MODES=120 DUDECT_N=200000 dudect-sign  # end-to-end signing timing smoke test

# Or run the full ref + avx2 + avx512 gate from the repository root.
cd Lithium-Code
./run_tests.sh
./integration_audit_matrix.sh      # integration audit only
./run_dudect_components.sh         # optional component timing smoke
```

Pick the security level at compile time with `LITHIUM_MODE` $\in$ `{120, 180, 260}` (default 120), e.g. `-DLITHIUM_MODE=180`.

## End-to-end signing dudect

The `ref/`, `avx2/`, and `avx512/` builds provide `make dudect-sign`, a small dudect-style Welch t-test harness for `crypto_sign_signature()`. The harness measures class 0 as a fixed secret key and class 1 as a key drawn from a pre-generated random secret-key pool, overrides `randombytes()` so measured signing consumes per-trace signing randomness prepared outside the timed region, and copies the selected key to a fixed buffer before timing so the timing class is not just a key-buffer address class. This test is intended as a long-term-secret timing smoke test; Gaussian rejection-sampling variance can still make individual runs noisy, so use `DUDECT_N` in the hundreds of thousands or millions on a pinned, low-noise machine before interpreting `|t|`. Native-Ubuntu retained evidence is available in `test/dudect_all.txt`: on 2026-06-05, Linux `6.11.0-26-generic`, `DUDECT_N=200000`, all `ref`, `avx2`, and `avx512` modes `120 180 260` stayed below the usual `4.5` threshold, with maximum observed `|t|=1.23` (`avx512` mode 180).

```bash
cd Lithium-Code/ref
make MODES='120 180 260' DUDECT_N=200000 dudect-sign

cd ../avx2
make MODES=120 DUDECT_N=200000 dudect-sign
make AVX2EXP=1 MODES=120 DUDECT_N=200000 dudect-sign

cd ../avx512
make MODES=120 DUDECT_N=200000 dudect-sign
make AVX512IFMAEXP=1 MODES=120 DUDECT_N=200000 dudect-sign
```

## API

NIST/PQClean-style signature interface (`ref/api.h`):

```c
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen, const uint8_t *sk);
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen, const uint8_t *pk);
/* plus the combined crypto_sign / crypto_sign_open wrappers */
```

Byte lengths are `CRYPTO_PUBLICKEYBYTES`, `CRYPTO_SECRETKEYBYTES`, `CRYPTO_BYTES`.

## Python-reproduction principle

Every design-level constant the C code uses is **generated** by a script in `tools/` **and** that same script doubles as a **bit-exact / numerically-exact golden reference** the C is checked against. No magic number is hand-written into the C without a script that reproduces and validates it. See [tools/README.md](tools/README.md) for the full table; in brief:

Each generator now rewrites a delimited `@@AUTOGEN:<tag>@@` region **inside** a hand-written file (via `tools/autogen.py`), instead of owning a whole header; `make tables` regenerates every region in place and the avx2 `make check-consts` gate re-runs them and diffs (drift $\Rightarrow$ failure).

| generator            | patches region in            | what it pins down |
|----------------------|------------------------------|-------------------|
| `gen_zetas.py`       | `ntt.h` (`@@ntt@@`)          | NTT twiddles + Montgomery constants |
| `gen_bounds.py`      | `params.h` (`@@bounds@@`)    | exact integer squared norm thresholds |
| `gen_approx_exp.py`  | `approx_exp.h` (`@@approx_exp@@`) | sampler minimax `exp(-a)` (B=8, degree-6, target $2^{-43}$) plus wide/small Taylor helpers |
| `gen_rcdt.py`        | `sampler.c` (`@@rcdt@@`, ref+avx2) | 93-bit RCDT sampler tables |
| `gen_rans_tables.py` | `irs_rans.h` (`@@rans_tables@@`) | the two rANS frequency tables |
| `SigSize.py`         | `irs_rans.h` (`@@rans_meta@@`) | `LITHIUM_RANS_RESERVED_BYTES` for $2^{-100}$ |
| `irs_ref.py`         | `irs_rans.h` (`@@irs_consts@@`) | repetition rate M, `LITHIUM_REJ_T = 33`, and the conditional-output error gate (`max_rel < $2^{-40}$`) |

Pure references / helpers (no emitted region): `ntt_ref.py`, `reduce_ref.py`, `rans.py`, `check_rans.py`, `autogen.py` (the region patch/read helper).

## Design highlights

- **NTT splits by level.** Full negacyclic NTT (`base_deg=1`, $q \equiv 1 \bmod 2n$) for Lithium-260; incomplete Kyber-style NTT (`base_deg=2`, deg-2 basemul, $q \equiv 1 \bmod n$) for Lithium-120/180. **Unsigned uint16 Montgomery reduction, R = 2^16** (canonical `[0,q)`; every NTT add/sub is reduced via `addm16`/`subm16` since `2q` can exceed `2^16`). The signed scheme domain (sampling, rejection, mod-2q, norm, packing) stays `int32`; conversion happens at the NTT boundary (`poly`⇄`poly16`).
- **Polynomial inversion `f^{-1}`** in the NTT residue domain (Fermat per-coeff for the full NTT; $2\times 2$ quadratic-ring inverse for the incomplete NTT); resample `f0` if non-invertible.
- **Iterative rejection sampling** via the exp-recurrence three-region method, series truncated to **T = 33** terms (`LITHIUM_REJ_T`). Our target-B parameters have $\alpha \approx 4$, where the alternating series converges far slower than the paper's $\alpha \approx 1.16$; the current `tools/irs_ref.py` check validates the conditional-output `max_rel` below $2^{-40}$ in all three modes (observed about $2^{-45}$). The keygen norm window `B_k' ≤ ‖s‖ < B_k` bounds $\gamma$ from below; see `docs/design-notes/Trunc.md`.
- **No slot-0 compression.** Every `z` component is split uniformly into a t-bit bit-packed tail + an rANS-coded head. Two rANS frequency tables (true rounded z-head; true MakeHint bucket crossing, including $|h| \ge 2$ tails), merged into a single rANS stream, `prob_bits = 10`, `L = 2^23`.
- **Binary challenge** (0/1, weight κ), transmitted as an n-bit string.
- **$2^{-100}$ rANS reserve** (`tools/SigSize.py`), with the zero-pad authenticated by the verifier (no signature malleability).

## Source map (`ref/`)

- **primitives:** `reduce.c`, `ntt.c` (+ Montgomery/twiddle constants in `ntt.h`), `rounding.c` (mod-2q HighBits/LowBits/MakeHint/UseHint), `sampler.c` (+ the RCDT tables, merged in), `approx_exp.c` (+ constants in `approx_exp.h`), `irs.c` (rejection step), `irs_rans.c` (rANS + the rANS tables / irs constants)
- **scheme:** `polyvec.c` (ring arithmetic incl. `poly_invert` + NTWE matrix/ vector ops), `packing.c` (raw + compact rANS signature), `sign.c` (keygen / sign / verify)
- **symmetric / RNG:** `fips202.c`, `symmetric-shake.c`, `randombytes.c`
- **tests** (`ref/test/`): `test_ntt`, `t_poly`, `t_exp`, `test_sampler`, `test_irs`, `test_mod2q`, `test_rans`, `test_keys`, `test_sign`, `test_kat`, `test_rans_failrate`; plus `speed*` benchmarks and `prof.{c,h}` (macro-gated profiler). KAT hashes recorded in `ref/vector/KAT.txt`.

## References

- D. J. Bernstein, K. Bhargavan, S. Bhasin, A. Chattopadhyay, T. K. Chia, M. J. Kannwischer, F. Kiefer, T. B. Paiva, P. Ravi, and G. Tamvada. *KyberSlash: Exploiting secret-dependent division timings in Kyber implementations.* January 2025.
- D. J. Bernstein. *Exploiting ML-DSA bugs.* Permanent ID ddd73b60e0c76be3f8236a304feba7d5805f9f60, June 1, 2026.
- J. Howe, T. Prest, T. Ricosset, and M. Rossi. *Isochronous Gaussian Sampling: From Inception to Implementation.* PQCrypto 2020, pages 53--71, DOI `10.1007/978-3-030-44223-1_4`.
- J. Gärtner. *Compact Lattice Signatures via Iterative Rejection Sampling.* Algorithm/specification notes: `docs/design-notes/IRS.md`.
- Parameter notes: `docs/design-notes/param.md`.
- rANS layer: `docs/design-notes/IRS_NIST_rANS.md`. Truncation depth: `docs/design-notes/Trunc.md`.
