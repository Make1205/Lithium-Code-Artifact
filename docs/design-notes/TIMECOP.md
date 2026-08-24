# TIMECOP / Valgrind Variable-Latency Notes

This note records the TIMECOP feasibility setup for Lithium-Code and its last
recorded result. Patched Valgrind is **dynamic** confirmation of the KyberSlash
variable-latency property; it complements, and does not replace, the
**static** object-code constant-time scan (`ct_scan_matrix.sh`) and the
statistical `dudect-sign` harness.

As of 2026-06-09 the full workflow builds and passes on native Ubuntu for the
`ref` and `avx2` backends (see [Backend Coverage](#backend-coverage),
[Environment](#environment), and [Result Report](#result-report)). The `avx512`
backend cannot be run under Valgrind and is covered statically instead.

## Tools

- One-shot driver (build-if-needed + positive control + smoke for each runnable
  backend × mode): `tools/timecop/run_all.sh`
- Build patched Valgrind: `tools/timecop/build_valgrind_varlat.sh`
- Positive-control secret-division detector:
  `tools/timecop/prove_varlat.sh --valgrind external/timecop/valgrind-varlat/install/bin/valgrind`
- Per-backend signing smoke (`--backend ref|avx2|avx512`, `--mode 120|180|260`):
  `tools/timecop/run_sign_smoke.sh --valgrind … --backend avx2 --mode 180`
- Constant-time cmov suppression: `tools/timecop/lithium-ct.supp`

From the repository root, `sh tools/timecop/run_all.sh` is the only command
needed; it builds the patched Valgrind on first use and then runs every probe.
`TIMECOP_BACKENDS`, `TIMECOP_MODES`, and `TIMECOP_TRY_AVX512=1` tune what runs.

## Methodology

The patched Memcheck (KyberSlash `valgrind-varlat` patch) adds
`--variable-latency-errors=yes`, which flags when a secret/undefined operand
reaches a **variable-latency instruction** (`div`/`mod`/...). This is the exact
class KyberSlash exploits.

- **Scope of the smoke = the variable-latency property only.** The gate is:
  *fail iff Memcheck reports `Variable-latency instruction operand ...`.*
  Secret-dependent **control flow** is owned separately by the object-code CT
  scan (`tools/ct_scan.py` disassembles secret-handling objects and rejects
  `idiv|div|sdiv|udiv|vdiv` and secret gather/scatter) and by `dudect-sign`.
- **`--undef-value-errors` must stay at its default (`yes`).** The
  variable-latency check is implemented on top of Memcheck's undefined-value
  instrumentation; `--undef-value-errors=no` silently disables it (an earlier
  version of these scripts passed that flag and therefore tested nothing).
- **`Memcheck:Cond` is suppressed; `Memcheck:Value` is not.** Lithium's
  sampler / `approx_exp` / IRS / packing use branchless constant-time idioms
  (full-table scans, masked selects, `ct_lt`/`ct_sel`/`ct_is_zero`) that GCC
  `-O2`/`-O3` lowers to `cmov`/`setcc`. Memcheck reports each secret-conditioned
  `cmov` as `Conditional jump or move depends on uninitialised value` — a
  known false positive on constant-time code — so those are suppressed via
  `lithium-ct.supp`. The patched Valgrind records a variable-latency hit as an
  `Err_Value` (the same kind as a plain use-of-uninitialised value), so a
  blanket `Memcheck:Value` suppression would blind the gate; we therefore keep
  `Memcheck:Value` un-suppressed and let the message-text grep isolate the
  variable-latency signal. Each smoke run also saves the **unsuppressed**
  stderr (`<backend>/out/timecop_smoke_<backend>_<mode>.full.err`) as an audit
  artifact.

## Backend Coverage

The harness `ref/test/timecop_smoke.c` is backend-agnostic and is compiled
against each backend's signing sources (the `avx2`/`avx512` builds carry
`-I../ref/test`, so no per-backend copy is needed). Each smoke is built with the
backend's real production flags (mirroring its Makefile `test_sign` list, minus
`randombytes.c` since the harness defines `randombytes()`).

| Backend | Build flags                                   | Under Valgrind 3.23 | Variable-latency owner here          |
| ------- | --------------------------------------------- | ------------------- | ------------------------------------ |
| `ref`   | `-O2` scalar                                  | runs                | this smoke (dynamic)                 |
| `avx2`  | `-O3 -mavx2 … -mno-avx512f` (pure AVX2)        | runs                | this smoke (dynamic)                 |
| `avx512`| `-O3 -mavx512f/bw/vl/dq/ifma -march=native`   | **SIGILL** — cannot | static only (`ct_scan` + `dudect`)   |

**Why `avx512` is skipped.** Valgrind 3.23's AVX-512 support is incomplete; it
cannot decode the EVEX-encoded instructions the backend emits (the built mode-120
binary contains ~7.6k `zmm`/k-mask instructions, and `-march=native` even
auto-vectorizes the "scalar" Keccak to AVX-512). The native binary runs fine, but
under Valgrind the guest aborts with `SIGILL — Illegal opcode` in
`shake256_init` before signing starts. This is a Valgrind limitation, not a
finding about the code, and it is not fixed by a newer Valgrind (AVX512IFMA in
particular remains unsupported) — so `run_sign_smoke.sh --backend avx512` reports
**SKIP** (exit 3), not FAIL. The avx512 no-variable-latency-division property is
still enforced statically: `tools/ct_scan.py` disassembles the avx512 objects
across `-O3`/`-Os` and `gcc`/`clang`, and `make -C avx512 dudect-sign` measures
the avx512 signing path statistically.

## Smoke Taint Scope

`ref/test/timecop_smoke.c` marks packed signing secret-key material and signing
randomness as undefined/secret. Sampler and IRS inputs are then tainted by
propagation through the normal signing data flow: the packed `key` and the
fresh signing randomness feed `rhoprime`, `rhoprime` feeds Gaussian sampling,
and `key || mu || seedc` feeds the IRS RNG. The taint offsets are computed from
the parameter macros (`SEEDBYTES`, `TRBYTES`, `SKSEEDBYTES`, `K`, `M`,
`POLYNOISE_PACKEDBYTES`, `POLYPK_PACKEDBYTES`), so they track `pack_sk`'s layout
in every mode. Public `rho`, `tr`, cached `bn`, public key, message, verifier
inputs, and the produced signature (after signing) are kept defined/declassified.

## Patch Sources

- KyberSlash artifacts page: `https://kyberslash.cr.yp.to/papers.html`
- Preferred patch: `https://kyberslash.cr.yp.to/valgrind-varlat-patch-20240808.txt`
  (native target: Valgrind 3.23.0)
- Fallback patch: `https://kyberslash.cr.yp.to/valgrind-varlat-patch-20250805.txt`

The patch is cut from the Valgrind git tree, so its first hunk edits
`.gitignore`, which the released source tarballs do not ship.
`build_valgrind_varlat.sh` filters out hunks that modify files absent from the
extracted tree (keeping new-file creations), applies the remainder with fuzz,
and verifies success by capability (`--variable-latency-errors` present and a
working build) rather than by `patch`'s raw exit code.

## Current Decision

TIMECOP is an **optional, integrated** audit workflow:

- It is **not** part of the default hard gate (`./run_tests.sh`,
  `./security_audit_matrix.sh`), because it requires a one-time local
  patched-Valgrind build (download + ~3 min compile).
- It is wired in as an **opt-in**: `TIMECOP=1 ./run_tests.sh` appends
  `tools/timecop/run_all.sh` (ref + avx2 across all modes) after the standard
  gates.
- The static object-code CT scan remains the mandatory owner of the
  no-variable-latency-division and no-secret-gather properties for **all**
  backends; the TIMECOP smoke is the dynamic cross-check for the backends
  Valgrind can run (`ref`, `avx2`).

## Environment

Recorded on the machine where the workflow was last validated (2026-06-09).

| Component        | Value                                                             |
| ---------------- | ----------------------------------------------------------------- |
| OS               | Ubuntu 24.04.2 LTS                                                 |
| Kernel           | Linux 6.11.0-26-generic, x86_64                                    |
| CPU              | 11th Gen Intel Core i7-11700K @ 3.60GHz (AVX2 + AVX-512; `nproc`=8)|
| RAM              | 31 GiB                                                             |
| Compiler         | GCC 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1)                        |
| Binutils / ld    | GNU ld 2.42                                                        |
| make             | GNU Make 4.3                                                       |
| patch            | GNU patch 2.7.6                                                    |
| autoconf/automake| 2.71 / 1.16.5                                                      |
| glibc            | 2.39 (`libc6-dbg` 2.39-0ubuntu8.7 — required by Valgrind)          |
| Built Valgrind   | `valgrind-3.23.0`, reports `--variable-latency-errors=no\|yes`     |

**Build dependencies installed** (`apt-get install`): `build-essential`,
`automake`, `autoconf`, `libtool`, `libc6-dbg`, `pkg-config`.

**Sources and checksums** (downloaded to
`external/timecop/valgrind-varlat/downloads/`):

| File                                  | SHA-256                                                            |
| ------------------------------------- | ----------------------------------------------------------------- |
| `valgrind-3.23.0.tar.bz2`             | `c5c34a3380457b9b75606df890102e7df2c702b9420c2ebef9540f8b5d56264d` |
| `valgrind-varlat-patch-20240808.txt`  | `b4fe3b9f1badcbec078e2dda08718ccc1571d3cf4f79c4b80ef6adb157d3999e` |

The build itself (`build_valgrind_varlat.sh`): downloads the tarball + patch,
filters the patch (drops the `.gitignore` hunk absent from the tarball), applies
it with fuzz, then `./configure --prefix=.../install`, `make -j8`, `make
install`; full log at `external/timecop/valgrind-varlat/build.log`. All build,
download, and run artifacts live under `external/` and `<backend>/out/`, both
gitignored, so nothing built here is committed.

## Result Report

**2026-06-09 — PASS (native Ubuntu).** Driver: `sh tools/timecop/run_all.sh`
(~23 s end-to-end once Valgrind is built: positive control + ref/avx2 × 3 modes).

### 1. Patched Valgrind build

```
$ external/timecop/valgrind-varlat/install/bin/valgrind --version
valgrind-3.23.0
$ ... --help | grep variable-latency
    --variable-latency-errors=no|yes check for variable latency errors [no]
```

### 2. Positive control — `prove_varlat.sh` → PASS

```
# secret_div:  storage = 100 / (13 | secret);   (secret = uninitialised malloc byte)
==…== Variable-latency instruction operand of size 4 is secret/uninitialised
==…==    at 0x109094: main (in …/external/timecop/probes/secret_div)

# public_div:  storage = 100 / 13;              (constant divisor)
<no output — no variable-latency report>
```

secret_div: **1** variable-latency report → detector works. public_div: **0**.

### 3. Lithium signing smoke — `run_sign_smoke.sh` → PASS (ref + avx2)

One keygen + sign + verify per backend × mode under the patched Memcheck, with
signing secret-key material and signing randomness marked secret. `varlat` is
the gate; `Cond`/`Value` are the unsuppressed-run counts (out of scope, see
Methodology). `siglen` equals each mode's `CRYPTO_BYTES`, confirming distinct
parameter sets are really exercised.

| Backend | Mode | `siglen` | `varlat` (gate, must be 0) | `Cond` | `Value` | Result |
| ------- | ---- | -------- | -------------------------- | ------ | ------- | ------ |
| ref     | 120  | 1187     | **0**                      | 14     | 9       | PASS   |
| ref     | 180  | 1618     | **0**                      | 14     | 9       | PASS   |
| ref     | 260  | 2348     | **0**                      | 13     | 9       | PASS   |
| avx2    | 120  | 1187     | **0**                      | 12     | 11      | PASS   |
| avx2    | 180  | 1618     | **0**                      | 13     | 11      | PASS   |
| avx2    | 260  | 2348     | **0**                      | 11     | 9       | PASS   |

`varlat = 0` in every runnable backend × mode → no secret operand reaches a
variable-latency instruction in the scalar or AVX2 signing path.

### 4. avx512 — SKIP (Valgrind cannot execute AVX-512)

```
$ run_sign_smoke.sh --backend avx512 --mode 120     # exit 3 (SKIP, not FAIL)
timecop avx512 mode-120 smoke: SKIP — Valgrind cannot execute this backend
  (Illegal opcode at address 0x11EBED); covered statically by ct_scan + dudect-sign

# avx512/out/timecop_smoke_avx512_120.full.err:
==…== Process terminating with default action of signal 4 (SIGILL)
==…==  Illegal opcode at address 0x11EBED
==…==    at 0x11EBED: lithium_fips202_ref_shake256_init (fips202.c:617)
```

The native avx512 binary runs fine; only Valgrind cannot decode its AVX-512
instructions. `TIMECOP_TRY_AVX512=1 sh tools/timecop/run_all.sh` reproduces this.

### Artifacts (gitignored, regenerated on each run)

- `external/timecop/valgrind-varlat/{build.log,BUILD_INFO,install/}`
- `external/timecop/probes/{secret_div,public_div}.{out,err}`
- `<backend>/out/timecop_smoke_<backend>_<mode>.{out,err,full.err,full.out}` —
  `.err` is the suppressed gate run; `.full.err` is the unsuppressed audit run.

Subsequent `run_*` invocations append their own command-log entries below this
section.
