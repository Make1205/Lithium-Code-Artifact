# Public-key unpack (`unpack_pk_bn`)

`unpack_pk` (and the identical b-hat loop in `unpack_sk`) used to read each packed coefficient into a 32-bit `tmp[N]` via the generic `bitunpack`, then call `freeze()` (a Barrett `reduce32` plus `caddq`, about 20 instructions) to land it in $[0, q)$. This is on the verifier hot path: the "unpack pk/sig" bucket is the largest non-NTT cost in verify. Both the 32-bit staging and the full `freeze()` are unnecessary because the data is fundamentally 16-bit and already nearly in range.

The work is now a single backend hook, `unpack_pk_bn(poly16 bn[M], const uint8_t *packed)`, declared in `polyvec.h`, namespaced in `namespace.h`, and implemented per backend in the forked `polyvec.c` (mirroring how `poly_ntt_import` is forked). `packing.c` stays shared and simply calls the hook from both `unpack_pk` and `unpack_sk`.

## `freeze()` is one conditional subtract here

A packed `QBITS`-bit field is some value $c \in [0, 2^{\text{QBITS}})$. In every mode $2^{\text{QBITS}} < 2q$ (mode 120: $65536 < 115202$; mode 180: $32768 < 61954$; mode 260: $65536 < 129026$), so $c \in [0, 2q)$ and a single conditional subtract maps it canonically into $[0, q)$:

$$c \mathrel{-}= q \mathbin{\&} {-}(c \ge q).$$

This is bit-exact with `freeze()` over the whole field domain, not just for honest keys: for $c \in [0, q)$ both are the identity, and for $c \in [q, 2^{\text{QBITS}}) \subset [q, 2q)$ both yield $c - q$. It was confirmed by brute force over all of $[0, 2^{\text{QBITS}})$ for every mode (zero mismatches). Because the result is identical for every input, the change is KAT-neutral: the recorded hashes are unchanged.

## Direct 16-bit unpacking

For `QBITS == 16` (modes 120 and 260) each coefficient occupies exactly two bytes with no cross-byte straddle, so unpacking is a plain little-endian `uint16` load followed by the conditional subtract. The SIMD forks load 16 (AVX2) or 32 (AVX-512) coefficients at a time and do the conditional subtract with `min_epu16(c, c - q)` (AVX2) or `cmpge_epu16_mask` plus `mask_sub` (AVX-512).

For `QBITS == 15` (mode 180) the fields straddle bytes. The scalar path is a fused least-significant-bit-first 15-bit reader straight to `uint16` (a `uint32` accumulator suffices, since $7 + 15 < 32$; no `uint64`, no `tmp[N]`, no `freeze`). The SIMD path extracts 8 coefficients per iteration: a `pshufb` spreads each field's up-to-three source bytes into a 32-bit lane, a per-lane `srlv` aligns it, a 15-bit mask isolates it, then the conditional subtract and a narrow back to `uint16`. The last group is handled by the scalar tail so the 16-byte broadcast load never reads past the final polynomial.

This removes the `tmp[N]` staging and replaces `freeze()` with one subtract on the scalar reference (which, being in the shared `polyvec.c` only for `ref`, also speeds the portable build), and adds a vectorized load for the AVX backends. Combined with the rANS format change in `IRS_NIST_rANS.md`, the verify "unpack pk/sig" bucket on the i7-11700K dropped from about $23600$ to about $10500$ cycles.
