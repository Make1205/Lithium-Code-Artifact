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

System GCC: -O3 -march=native -std=gnu11. Ref additionally disables automatic
loop/SLP vectorization; algorithm SIMD intrinsics are selected only in the
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
