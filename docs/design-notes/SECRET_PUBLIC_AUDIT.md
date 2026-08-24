# Secret/Public Audit Baseline

This note is the secret/public classification baseline for Lithium-Code. Any future division, modulus, lookup, branch, gather allowlist, rejection timing decision, or declassification in secret-handling code must cite one of the entries below or extend this note in the same style before the code change lands.

## Classification Labels

- `secret`: Derived from signing secret material, signing randomness, Gaussian/noise samples, IRS randomness, or unreleased intermediate values whose distribution or value must not control variable-latency operations.
- `public`: Known to the verifier, attacker-controlled, deterministic from public inputs, or already serialized in public artifacts.
- `public length`: A byte count, mode constant, loop bound, fixed retry budget, or one-shot buffer size that may drive branches or divisions because it is independent of secret values.
- `declassified`: A value derived from secret data but intentionally released by the scheme through the public key, signature, or accept/reject output after the documented check.

## Keygen

- `seedbuf`, expanded key seeds before packing, `rpseed`, `f0`, `f`, `s0`, `e`, the assembled secret vector `s`, keygen norm `nrm`, and invertibility retry state are `secret` until the public and secret keys are packed.
- `rho`, expanded `A0`, loop bounds, mode constants, polynomial dimensions, and pack sizes are `public` or `public length`.
- `b` and canonical `bn` are `declassified` when packed into `pk` and copied into `sk`; after that point `bn` may be treated as public cached key material.
- Keygen rejection timing is deployment-classified as a provisioning-time retry, not a per-message side channel, but operations inside each attempt still follow `secret` arithmetic rules.

## Sign

- Unpacked `key`, `s`, signing randomness `rnd`, `rhoprime`, Gaussian `y`, `z` before signature release, IRS RNG state, norm value, rejection decisions, and rANS overflow retry state are `secret`.
- Message bytes and message length are attacker-controlled `public`; `tr`, `mu`, and challenge hash input length are public once bound to the public key and message transcript.
- `rho`, cached public `bn`, expanded `A0`, imported `bh`, mode constants, `LITHIUM_BYTES`, and rANS reserve lengths are `public` or `public length`.
- The challenge `c`, hints `h`, `z1`/`zlow`, rANS stream length, and zero padding are `declassified` only through `pack_sig`; before successful packing they remain secret-derived and must not drive variable-latency operations except through documented fixed-format writers.

## Verify

- Public key, message, signature bytes, signature length, unpacked `c`, `h`, `z1`, reconstructed `w0bits`, and `A0` are `public` because verification operates on adversary-provided or public-key-derived data.
- Verify accept/reject codes and early exits are public parser behavior; they must be deterministic and covered by negative corpus cases, but they are not secret-dependent branches.
- Reconstructed `z2` and norm `nrm` in verify are `public` because they are functions of public signature and public key data.

## Sampler

- Gaussian and noise sample values, sign bits, CDT decisions, approximate-exp acceptance values, refill state that depends on accepted samples, and generated secret-key or signing vectors are `secret`.
- SHAKE output block counts derived from fixed sampler batch sizes, initial refill sizes, lane counts, and mode constants are `public length`.
- RCDT tables, approximate-exp constants, and fixed distribution tables are `public`; indexes into them are allowed only through full scans, fixed-lane SIMD selection, or documented constant-time table selection.

## IRS

- IRS input secret vector, challenge-expanded secret contribution, random bits from the IRS RNG, rejection choices, and adjusted `z` values are `secret` until signature packing declassifies the accepted output.
- Challenge polynomial support size `kappa`, polynomial dimensions, and loop trip counts are `public length`.
- Any future shortcut in IRS rejection or coefficient adjustment must cite this section and prove that it does not expose secret-dependent variable latency.

## Challenge

- Sign-side `buckets`, `w0bits`, `mu`, `seedc`, and `c` are secret-derived until the final signature declassifies `c`; hash absorption lengths and mode constants are `public length`.
- Verify-side challenge recomputation inputs are `public`; challenge mismatch is public reject behavior.
- Challenge packing is fixed-size and public after serialization; malformed challenge bytes are parser-negative corpus inputs.

## rANS

- Sign-side z-head and hint symbols, rANS encoder states, stream bytes before publication, rANS used length, and overflow retry decisions are `secret` until the signature is successfully emitted.
- rANS frequency tables, CDF tables, reciprocal constants, interleaving count, reserve size, and final fixed signature size are `public` or `public length`.
- Verify-side rANS bytes, decoded symbols, noncanonical-state checks, reserve overflow checks, and padding checks are `public` parser behavior and must stay in the parser-negative corpus.

## Packing

- Key packing treats secret-key seed and secret polynomials as `secret`; public key bytes and cached `bn` are `declassified` public-key material.
- Sign-side signature packing handles secret-derived `c`, `h`, and `z1` until the complete fixed-size signature is emitted; zero padding is `declassified` and authenticated by verify.
- Verify-side unpacking treats all signature bytes as attacker-controlled `public`; branches on malformed lengths, padding, rANS decode failure, and signature size are public parser behavior.

## NTT Boundary

- Secret-key and signing vectors entering NTT or inverse NTT routines are `secret`; public matrix expansion `A0`, public-key `bn`, and verifier-side imported NTT data are `public`.
- Backend-native NTT slot order is an internal representation detail, not a declassification event; only `pack_pk`, `pack_sk`, and `pack_sig` define serialized declassification boundaries.
- SIMD gathers or shuffles may be allowlisted only when their addresses depend on fixed permutations, public lane layout, public SHAKE lane pointers, or mode constants; secret-indexed gather remains forbidden.

## Current Allowlist Citations

- `ct_scan_matrix.sh` allowlisted public divisions in SHAKE block-count code and public Gaussian initial-refill block counts cite the `Sampler` and `Packing` public-length entries.
- SIMD SHAKE gather allowlists cite the `NTT Boundary` and `Sampler` entries because lane pointer/offset selection is public and fixed by the backend schedule.
- Verify parser early exits cite the `Verify`, `rANS`, and `Packing` public parser behavior entries.
