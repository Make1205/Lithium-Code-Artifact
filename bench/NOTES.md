# Benchmark notes

Initial worktree: clean. `bench` was created from the recorded initial commit.
Only files under `bench/` are added. Algorithm source, parameters and original
build files remain unchanged. The runner extracts the selected commit with local
`git archive`, so an older upstream revision can be measured without resetting
the working branch or rewriting files. No network is used at run time.

Initial commit: `81389d1971ef64f6469167ab10855e3594ea7729`.
Tested commit: `81389d1971ef64f6469167ab10855e3594ea7729`.

## Timing source and measurement

The four unmodified timing files come from pq-crystals/dilithium
`bba9534a3c5cbef66481cd6518d786e23ae601b5`, `ref/test/` (there is no root
`test/` directory at that revision). They match the Antrag bench copies byte
for byte. RDTSC is used; USE_RDPMC is not defined. The unchanged print_results
function receives iterations+1 timestamps for exactly iterations operations,
subtracts its own counter overhead, and computes the median. Thousand cycles
means RDTSC ticks divided by 1000. Loop/call overhead is included consistently.

Historical setup (superseded by the calibration audit below): system GCC
-O3 -march=native -std=gnu11; ref disabled automatic loop/SLP vectorization; algorithm SIMD intrinsics are selected only in the
upstream SIMD implementation. Exact flags are recorded per variant in
config.json, metadata.json and build logs. Messages have 32 bytes. Signing and
verification receive serialized byte keys and include all decoding/expansion
performed by the upstream byte API. KeyGen returns serialized keys.
Expanded-key Sign is reported only when a supported API is available.

QUICK=1: 11 timing operations, upstream KAT/self-tests, 100 generated signatures.
Default FULL: 10001 timing operations, upstream KAT/self-tests, 10000 generated
signatures. Only QUICK is executed in cloud. Full-scale runs are for the user's
local machine. No test assertions are disabled. The correctness and length
statistics reuse the same signature pool, using a fixed-seed key per build/set.
All signatures self-verify, cross-verify when another build exists, and are
checked for rejection after changing the message. Equality expectations are
configured explicitly; floating-point implementations can use cross verification
when byte equality is inappropriate. No result is declared complete unless all
required commands finish with zero status. The CSV contains medians and sizes;
raw timer output and correctness evidence are retained.

A bench-only xorshift64* byte stream (seed 1234567, fixed update and byte order in
driver.c) supplies reproducible test entropy where the API accepts an external
provider. This is for correctness vectors, not production randomness. For NIST
APIs the timing provider is the original randombytes implementation, renamed only
in its provider translation unit and routed through a bench wrapper; the algorithm
call sites remain unchanged. Callback APIs use the repository SHAKE generator
with fixed test seeds and OS-seeded timing state as documented below.

Missing system-library headers/libraries are compile/link checked at startup,
before algorithm builds. The runner never installs packages. Bash, Python 3.8+,
Git and system GCC are required; additional apt dependencies are listed below.

## Scheme-specific configuration and format

Default integer-exp ref and AVX2. Signature is fixed/padded; actual length = fixed length minus reserved rANS block size plus its encoded two-byte little-endian rANS length. Ratios use user-provided COPA calibration medians; QUICK ratios are format checks only. No expanded-key signing API is exposed. Timing uses the upstream native randombytes provider; fixed-seed correctness uses a bench-only deterministic byte provider.

## Native-source compilation changes

No algorithm source patches. All source copies come from the tested commit.
Additional source selection, preprocessor options and include paths appear in
config.json. Configuration-only exceptions, if any, are documented below.

## Timing file checksums

* cpucycles.c: `31cc4c16b1e896b2cf7cc154aa21430329c0047cf9898a1865f9bfe9904f25e4`
* cpucycles.h: `7f1e2b68cf3e3c76076364af3c5350de61442450c477aec809da8de10e131838`
* speed_print.c: `533fc7bd12204078a5390cfb58c05d550b214850876609f8e98f7377b9d4f1b3`
* speed_print.h: `3946b0226c63125664eb443a20d43563e6d1ba7fa3be5c1b952cd29d1477b2ef`

Additional apt dependencies: none beyond build-essential, python3 and git.

The requested upstream commit equals the initial commit, so no revision fallback
was needed. Default integer exp is measured, not the optional AVX2Exp/AVX-512
variants. AVX2 needs BMI2 and POPCNT. It uses USE_AVX2_NTT, USE_AVX2_SAMPLER and
USE_AVX2_SHAKE4X and -mno-avx512f to retain the upstream pure-AVX2 configuration.
The native test/test_kat.c hashes are checked against the three recorded values
in avx2/Makefile, not just the test exit status. The default KAT is ten cases.

The external driver preincludes the backend's config.h before api.h to select
the correct namespace through symlinked public headers. Test header directories
are included because the upstream implementation includes prof.h. These are
compiler options only. No source or KAT fixture was patched.

Paper q/pk/fixed-sig values match all three sets. Actual unpadded size comes
from the two-byte rANS length within the signature. It removes only unused
reserved bytes in that internal block, preserving the encoded tail after it.
Calibration ratios are current median / user-provided median; QUICK ratios
are printed solely to validate the format, not to validate machine performance.

Before any algorithm build, the runner also compile/links a standard C
header/runtime/libm probe. Failure reports apt package build-essential
(including libc6-dev). Exact extra per-build flags are printed near the
SUMMARY heading as well as retained in metadata and logs. This final
preflight/report change leaves the measured driver and QUICK values unchanged.


## COPA local_check audit and calibration gate (2026-11-10)

COPA was read with `git show` at
`3c676bd2acba2d1d8d2f9f564f097f0e5c4588ac`; its checkout was not modified.
LITHIUM is still the unmodified upstream
`81389d1971ef64f6469167ab10855e3594ea7729`.
Evidence, line-numbered excerpts and source diffs are in `calibration/`.
Git symlinks were resolved before comparing source contents.

| Item | COPA local_check at 3c676bd | Previous LITHIUM bench at 36fe62a | This revision |
|---|---|---|---|
| Compiler and ref options | System GCC; effective `-O3 -march=native -std=c99 -Wall -Wextra`; no vectorization suppression | `-O3 -march=native -std=gnu11 -fno-tree-vectorize -fno-tree-slp-vectorize` | `-O3 -march=native -std=c99`; removed both extra suppression options |
| Native Makefiles | ref defaults include warning/stack-protector options, but local_check overrides CFLAGS; AVX2 uses common/copa-backend.mk and copa-benchmark.mk | Builds selected upstream translation units directly | Same direct source build, exact commands retained in build logs; no stack-protector claim inferred from unused defaults |
| AVX2 selectors | `-mavx2 -mno-avx512f -mbmi2 -mpopcnt`, USE_AVX2_NTT/SAMPLER/SHAKE4X | Same selectors | Unchanged; no optional AVX2Exp/AVX-512 path |
| Sources | ref ntt/reduce/polyvec/rounding/sampler/approx_exp/irs/packing/irs_rans/sign/symmetric-shake/fips202/randombytes; AVX2 also fips202x4, f1600x4 assembly and modulus-specific NTT assembly/constants | Corresponding upstream sources, no COPA patches | Unchanged upstream sources; source identity differences detailed below |
| Bench macros | BENCH_LITHIUM, BENCH_COUNT_ATTEMPTS, LITHIUM_MODE and linker --wrap=randombytes | LITHIUM_MODE plus config.h preinclude for namespaces | Same namespace/mode setup; no BENCH_COUNT_ATTEMPTS source instrumentation |
| Sign interval | Byte-sk API: unpack/decode sk, message hashing, matrix expansion, rejection sampling and encoding; instrumented attempt counter and two internal serialized timestamps | Same byte-sk API work, without COPA instrumentation | Same upstream byte-sk API; only API call inside interval |
| Verify interval | Byte-pk API: unpack pk/sig, hash pk and message, expand A0 and verification arithmetic | Same byte-pk API work | Unchanged; no predecoded pk or matrix cache supplied by driver |
| Message | 59 zero bytes | 32 zero bytes | 32 zero bytes, as required |
| Entropy for timing | Script defaults BENCH_SEED to copa-local-v3-fixed-corpus; wrapper SHAKE256(seed bytes concatenated with LE64 call counter), starting at zero | Native OS randombytes | Native OS by default; CALIBRATION=1 sets the same default seed and uses the identical SHAKE wrapper construction, entirely in bench/rng.c |
| Sample population | Fresh KeyGen, Sign and Verify triple each round | Separate operation loops, Sign on one last generated key, Verify on one last signature | Fresh byte-key/signature triple each round |
| Warmup | 20 complete triples discarded | One initial keygen/sign, no 20-triple warmup | 20 complete triples discarded |
| Counter | LFENCE; RDTSCP; LFENCE around calls; asserts inside interval; no counter overhead subtraction | Dilithium RDTSC files, adjacent timestamps in operation loops; print_results subtracts counter overhead | Unmodified Dilithium timing files retained; explicit start/end per API, duration arrays passed as cumulative timestamps to unchanged print_results; assertions outside interval |
| Median | Sort each operation duration vector, index N/2 (10001 samples) | Unchanged Dilithium print_results median and overhead subtraction | Same unchanged print_results; odd N=10001 FULL or 11 QUICK |

The extra ref vectorization suppression was an actual compiler mismatch and is
removed. The timed byte-key parsing, ExpandA and required message hashing were
already included in both drivers; omitting them does not explain the ratios.
The single-key versus fresh-key corpus and OS versus seeded entropy can change
signing rejection distributions. These are identified confounders, not measured
attributions of the observed 1.561 ratio. Both 32+64 and 59+64 bytes fit one
SHAKE256 rate block; message length alone is not an established explanation for
that large difference. AVX2 KeyGen near 1 does not establish timer equivalence.

COPA also changes executable code compared with original LITHIUM: both sign.c
files add attempt/timing instrumentation; avx2/fips202x4.c changes partial-block
loads to zero-initialized memcpy tails and changes output stores. The preserved
diffs show these changes. The SHA256 of COPA's fips202x4.c is
`6786911e551533f68ddb3c563959d0d3a985fef06b7b1ea2853ef3f494c9a306`;
upstream is `bedf5f7a5d90203a21a6ae46a76fb3f0ef043e8b077252deacbfaf123014bd60`.
They are not identical benchmark source trees, even if cryptographic outputs
match. We do not import these patches or modify any algorithm code to chase the
historical numbers. No claim is made that COPA's cryptographic results are wrong.

### User-supplied FULL reference, exact cycles

CPU i7-10700K, GCC 11.4.0, pinned CPU 0, samples 10001, warmup 20,
59-byte messages, byte keys. BENCH_SEED was unset by the user; local_check sets
the default above. These exact numbers, not rounded thousand-cycle values, are
used in config.json as the ratio denominators.

| Set | ref KeyGen/Sign/Verify (cycles) | AVX2 KeyGen/Sign/Verify (cycles) |
|---|---|---|
| 120 | 731859/746648/170261 | 283643/348866/89493 |
| 180 | 1021189/1185004/250273 | 386785/520265/120965 |
| 260 | 1426175/1522873/310186 | 557425/645780/143338 |

### Acceptance status and local command

The unified rules are in [PROTOCOL.md](PROTOCOL.md). The revised driver aligns
flags, fresh-key population, warmup and optional seeded entropy while keeping
32-byte messages, original algorithm sources and the required unmodified timing
files. It intentionally does not masquerade as the historical COPA workload.

Run `CALIBRATION=1 bash bench/run.sh` on the stated local CPU/GCC. FULL runs all
normal KAT/correctness/length/timing work and additionally writes calibration.json.
It returns exit 2 unless all 18 medians are within 3% and CPU/GCC/affinity/seed
match the supplied baseline. A passing gate means numerical agreement within
3% for these two disclosed workloads; it does not prove source/workload identity.
QUICK always says NOT_VALIDATED_QUICK and cannot accept the numerical target.
The cloud AMD EPYC/GCC14 run passed six upstream KATs and 100-signature equality,
self-verification and cross-verification for each parameter set. It is not a
same-machine FULL calibration. The run's metadata retains benchmark hashes and
its dirty pre-commit status so its provenance is reviewable.

The user's instruction says to stop if this calibration cannot be reproduced.
The numerical condition is currently unvalidated, and source/message/timer
mismatches prevent declaring it reproduced by inspection. Other repository
changes (including VecFalcon and c-fn-dsa) are therefore deferred at this gate;
no all-repository rerun is claimed complete. The local gate provides the next
concrete measurement without changing COPA or any algorithm source.
