/*
 * packing.c - serialization for IRS (RAW validation packers + compact
 * on-wire).
 */
#include "packing.h"

#include <stdint.h>
#include <string.h>

#include "irs_rans.h"
#include "polyvec.h"
#include "reduce.h"

/* ---- generic little-endian bit packers ---- */
static void bitpack(uint8_t *out, const uint32_t *v, int n, int bits)
{
    /* LSB-first bit-accumulator writer: the byte-emit schedule depends
     * only on the loop index (not the data), so there are NO
     * data-dependent branches (the old per-bit `if (x & (1<<b))`
     * mispredicted ~50% of the time) and no memset. Byte-identical output
     * -> KAT-neutral.  bits <= 32 keeps acc (<= 7+32 bits) within 64; the
     * final partial byte is zero-padded in its high bits. */
    const uint32_t mask = (bits >= 32) ? 0xffffffffu : ((1u << bits) - 1u);
    uint64_t acc = 0;
    int accbits = 0;
    size_t outpos = 0;
    for (int i = 0; i < n; i++) {
        acc |= (uint64_t)(v[i] & mask) << accbits;
        accbits += bits;
        while (accbits >= 8) {
            out[outpos++] = (uint8_t)acc;
            acc >>= 8;
            accbits -= 8;
        }
    }
    if (accbits > 0)
        out[outpos] =
            (uint8_t)acc; /* trailing partial byte (zero-padded) */
}
static void bitunpack(uint32_t *v, const uint8_t *in, int n, int bits)
{
    /* LSB-first bit-accumulator reader: byte-read schedule depends only on
     * the loop index, so no data-dependent branches (vs the old per-bit
     * `if`).  Reads exactly ceil(n*bits/8) bytes -> byte-identical inverse
     * of bitpack (KAT-neutral). */
    const uint32_t mask = (bits >= 32) ? 0xffffffffu : ((1u << bits) - 1u);
    uint64_t acc = 0;
    int accbits = 0;
    size_t inpos = 0;
    for (int i = 0; i < n; i++) {
        while (accbits < bits) {
            acc |= (uint64_t)in[inpos++] << accbits;
            accbits += 8;
        }
        v[i] = (uint32_t)acc & mask;
        acc >>= bits;
        accbits -= bits;
    }
}

void pack_pk(uint8_t *pk, const uint8_t rho[SEEDBYTES], const poly16 bn[M])
{
    memcpy(pk, rho, SEEDBYTES);
    uint8_t *p = pk + SEEDBYTES;
    uint32_t tmp[N];
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < N; k++)
            tmp[k] = (uint32_t)bn[i].coeffs[k];
        bitpack(p, tmp, N, QBITS);
        p += POLYPK_PACKEDBYTES;
    }
}
void unpack_pk(uint8_t rho[SEEDBYTES], poly16 bn[M], const uint8_t *pk)
{
    memcpy(rho, pk, SEEDBYTES);
    unpack_pk_bn(bn, pk + SEEDBYTES);
}

static void pack_noise(uint8_t *out, const poly *p)
{
    uint32_t tmp[N];
    uint32_t mask = (1u << NOISE_BITS) - 1u;
    for (int k = 0; k < N; k++)
        tmp[k] = ((uint32_t)p->coeffs[k]) & mask;
    bitpack(out, tmp, N, NOISE_BITS);
}
static void unpack_noise(poly *p, const uint8_t *in)
{
    uint32_t tmp[N];
    bitunpack(tmp, in, N, NOISE_BITS);
    const int32_t sgnbit = 1 << (NOISE_BITS - 1);
    const int32_t range = 1 << NOISE_BITS;
    for (int k = 0; k < N; k++) {
        int32_t x = (int32_t)tmp[k];
        int32_t m = -((x & sgnbit) >> (NOISE_BITS - 1));
        p->coeffs[k] = x - (range & m);
    }
}

void pack_sk(uint8_t *sk, const uint8_t rho[SEEDBYTES],
             const uint8_t tr[TRBYTES], const uint8_t key[SKSEEDBYTES],
             const poly *f0, const poly s0[ELL], const poly e[M],
             const poly16 bn[M])
{
    uint8_t *p = sk;
    memcpy(p, rho, SEEDBYTES);
    p += SEEDBYTES;
    memcpy(p, tr, TRBYTES);
    p += TRBYTES;
    memcpy(p, key, SKSEEDBYTES);
    p += SKSEEDBYTES;
    pack_noise(p, f0);
    p += POLYNOISE_PACKEDBYTES;
    for (int j = 0; j < ELL; j++) {
        pack_noise(p, &s0[j]);
        p += POLYNOISE_PACKEDBYTES;
    }
    for (int i = 0; i < M; i++) {
        pack_noise(p, &e[i]);
        p += POLYNOISE_PACKEDBYTES;
    }
    uint32_t tmp[N];
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < N; k++)
            tmp[k] = (uint32_t)bn[i].coeffs[k];
        bitpack(p, tmp, N, QBITS);
        p += POLYPK_PACKEDBYTES;
    }
}
void unpack_sk(uint8_t rho[SEEDBYTES], uint8_t tr[TRBYTES],
               uint8_t key[SKSEEDBYTES], poly *f0, poly s0[ELL], poly e[M],
               poly16 bn[M], const uint8_t *sk)
{
    const uint8_t *p = sk;
    memcpy(rho, p, SEEDBYTES);
    p += SEEDBYTES;
    memcpy(tr, p, TRBYTES);
    p += TRBYTES;
    memcpy(key, p, SKSEEDBYTES);
    p += SKSEEDBYTES;
    unpack_noise(f0, p);
    p += POLYNOISE_PACKEDBYTES;
    for (int j = 0; j < ELL; j++) {
        unpack_noise(&s0[j], p);
        p += POLYNOISE_PACKEDBYTES;
    }
    for (int i = 0; i < M; i++) {
        unpack_noise(&e[i], p);
        p += POLYNOISE_PACKEDBYTES;
    }
    unpack_pk_bn(bn, p);
}

/* ---- RAW signature (validation) ---- */
void pack_sig_raw(uint8_t *sig, const poly *c, const poly h[M],
                  const poly z1[ELL + 1])
{
    uint8_t *p = sig;
    memset(p, 0, CHALLENGE_PACKEDBYTES);
    for (int k = 0; k < N; k++)
        if (c->coeffs[k])
            p[k >> 3] |= (uint8_t)(1u << (k & 7));
    p += CHALLENGE_PACKEDBYTES;
    uint32_t tmp[N];
    uint32_t hmask = (1u << W1_BITS) - 1u;
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < N; k++)
            tmp[k] = ((uint32_t)h[i].coeffs[k]) & hmask;
        bitpack(p, tmp, N, W1_BITS);
        p += HINT_PACKEDBYTES;
    }
    for (int i = 0; i < ELL + 1; i++) {
        for (int k = 0; k < N; k++) {
            int16_t v = (int16_t)z1[i].coeffs[k];
            p[2 * k] = (uint8_t)(v & 0xff);
            p[2 * k + 1] = (uint8_t)((v >> 8) & 0xff);
        }
        p += Z1_PACKEDBYTES;
    }
}
void unpack_sig_raw(poly *c, poly h[M], poly z1[ELL + 1],
                    const uint8_t *sig)
{
    const uint8_t *p = sig;
    for (int k = 0; k < N; k++)
        c->coeffs[k] = (p[k >> 3] >> (k & 7)) & 1;
    p += CHALLENGE_PACKEDBYTES;
    uint32_t tmp[N];
    int32_t sgn = 1 << (W1_BITS - 1);
    for (int i = 0; i < M; i++) {
        bitunpack(tmp, p, N, W1_BITS);
        for (int k = 0; k < N; k++) {
            int32_t x = (int32_t)tmp[k];
            if (x & sgn)
                x -= (1 << W1_BITS);
            h[i].coeffs[k] = x;
        }
        p += HINT_PACKEDBYTES;
    }
    for (int i = 0; i < ELL + 1; i++) {
        for (int k = 0; k < N; k++)
            z1[i].coeffs[k] = (int16_t)((uint16_t)p[2 * k] |
                                        ((uint16_t)p[2 * k + 1] << 8));
        p += Z1_PACKEDBYTES;
    }
}

/* ---- compact (on-wire) signature: rANS(z-head ++ hint) + raw z-lo ----
 */
#define NZH ((ELL + 1) * N)
#define NHN (M * N)
#define HALF_KSPLIT (1 << (T - 1))
#define PACK_RANS_PB RANS_PROB_BITS
#define PACK_RANS_PSCALE (1u << PACK_RANS_PB)
#define PACK_RANS_STREAM_MASK RANS_INTERLEAVED_MASK

static inline unsigned pack_rans_stream_index(size_t t)
{
    return (unsigned)(t & PACK_RANS_STREAM_MASK);
}

static inline int pack_rans_state_in_range(uint32_t x)
{
    return x >= RANS_L && x < (RANS_L << 8);
}

static inline void pack_sig_split_z(int32_t z, int32_t *head, int32_t *low)
{
    *head = (z + HALF_KSPLIT) >> T;
    *low = z - (*head * KSPLIT);
}

static inline void pack_rans_put(uint32_t *x, uint8_t *out, size_t *pos,
                                 uint32_t freq, uint32_t rcp, uint32_t rsh,
                                 uint32_t bias)
{
    uint32_t x_max = ((RANS_L >> PACK_RANS_PB) << 8) * freq;
    while (*x >= x_max) {
        out[--(*pos)] = (uint8_t)(*x & 0xff);
        *x >>= 8;
    }
    uint32_t q = (uint32_t)(((uint64_t)*x * rcp) >> 32) >> rsh;
    *x = *x + bias + q * (PACK_RANS_PSCALE - freq);
}

static inline int pack_sig_encode_hint(uint32_t *x, uint8_t *out,
                                       size_t *pos, int32_t sym)
{
    int slot = (int)(sym - RANS_HINT_LO);
    if (slot < 0 || slot >= RANS_HINT_N)
        return -2;
    pack_rans_put(x, out, pos, RANS_HINT_FREQ[slot], RANS_HINT_RCP[slot],
                  RANS_HINT_RSH[slot], RANS_HINT_BIAS[slot]);
    return 0;
}

static inline int pack_sig_encode_zhead(uint32_t *x, uint8_t *out,
                                        size_t *pos, int32_t sym)
{
    int slot = (int)(sym - RANS_ZH_LO);
    if (slot < 0 || slot >= RANS_ZH_N)
        return -2;
    pack_rans_put(x, out, pos, RANS_ZH_FREQ[slot], RANS_ZH_RCP[slot],
                  RANS_ZH_RSH[slot], RANS_ZH_BIAS[slot]);
    return 0;
}

static inline void unpack_sig_decode_zhead(uint32_t *x, poly *z1i, int k)
{
    uint32_t val = *x & (PACK_RANS_PSCALE - 1);
    unsigned slot = RANS_ZH_SLOT[val];
    z1i->coeffs[k] = RANS_ZH_LO + (int)slot;
    *x = RANS_ZH_FREQ[slot] * (*x >> PACK_RANS_PB) + val -
         RANS_ZH_CDF[slot];
}

static inline void unpack_sig_decode_hint(uint32_t *x, poly *hi, int k)
{
    uint32_t val = *x & (PACK_RANS_PSCALE - 1);
    unsigned slot = RANS_HINT_SLOT[val];
    hi->coeffs[k] = RANS_HINT_LO + (int)slot;
    *x = RANS_HINT_FREQ[slot] * (*x >> PACK_RANS_PB) + val -
         RANS_HINT_CDF[slot];
}

static inline int unpack_sig_renorm(uint32_t *x, const uint8_t *in,
                                    size_t in_len, size_t *bp)
{
    while (*x < RANS_L) {
        if (*bp >= in_len)
            return -1;
        *x = (*x << 8) | in[(*bp)++];
    }
    return 0;
}

static void unpack_sig_clear_outputs(poly h[M], poly z1[ELL + 1])
{
    memset(h, 0, (size_t)M * sizeof h[0]);
    memset(z1, 0, (size_t)(ELL + 1) * sizeof z1[0]);
}

static int unpack_sig_rans_decode(poly h[M], poly z1[ELL + 1],
                                  const uint8_t *in, size_t in_len)
{
    if (in_len < (size_t)(4 * RANS_INTERLEAVED_STREAMS))
        return -1;
    uint32_t x[RANS_INTERLEAVED_STREAMS];
    size_t bp = 0;
    for (int s = RANS_INTERLEAVED_STREAMS - 1; s >= 0; s--) {
        uint32_t v = 0;
        for (int b = 0; b < 4; b++)
            v = (v << 8) | in[bp++];
        x[s] = v;
        if (!pack_rans_state_in_range(x[s]))
            return -1;
    }

    size_t t = 0;
    for (int i = 0; i < ELL + 1; i++) {
        for (int k = 0; k < N; k++) {
            unsigned s = pack_rans_stream_index(t++);
            unpack_sig_decode_zhead(&x[s], &z1[i], k);
            if (unpack_sig_renorm(&x[s], in, in_len, &bp) != 0)
                return -1;
        }
    }
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < N; k++) {
            unsigned s = pack_rans_stream_index(t++);
            unpack_sig_decode_hint(&x[s], &h[i], k);
            if (unpack_sig_renorm(&x[s], in, in_len, &bp) != 0)
                return -1;
        }
    }
    if (bp != in_len)
        return -1;
    for (int s = 0; s < RANS_INTERLEAVED_STREAMS; s++)
        if (x[s] != RANS_L)
            return -1;
    return 0;
}

static int pack_sig_rans_encode(uint8_t *out, size_t *out_len,
                                size_t max_bytes, const poly h[M],
                                const poly z1[ELL + 1])
{
    uint32_t x[RANS_INTERLEAVED_STREAMS];
    for (int s = 0; s < RANS_INTERLEAVED_STREAMS; s++)
        x[s] = RANS_L;
    size_t pos = max_bytes;
    size_t ntot = NZH + NHN;
    size_t t = ntot;

    for (int i = M; i-- > 0;) {
        for (int k = N; k-- > 0;) {
            --t;
            unsigned s = pack_rans_stream_index(t);
            if (pos < (size_t)(4 * RANS_INTERLEAVED_STREAMS + 4))
                return -2;
            if (pack_sig_encode_hint(&x[s], out, &pos, h[i].coeffs[k]) !=
                0)
                return -2;
        }
    }
    for (int i = ELL + 1; i-- > 0;) {
        for (int k = N; k-- > 0;) {
            int32_t head, low;
            pack_sig_split_z(z1[i].coeffs[k], &head, &low);
            --t;
            unsigned s = pack_rans_stream_index(t);
            if (pos < (size_t)(4 * RANS_INTERLEAVED_STREAMS + 4))
                return -2;
            if (pack_sig_encode_zhead(&x[s], out, &pos, head) != 0)
                return -2;
        }
    }
    for (int s = 0; s < RANS_INTERLEAVED_STREAMS; s++)
        for (int b = 0; b < 4; b++) {
            if (pos == 0)
                return -2;
            out[--pos] = (uint8_t)(x[s] & 0xff);
            x[s] >>= 8;
        }
    size_t used = max_bytes - pos;
    memmove(out, out + pos, used);
    *out_len = used;
    return 0;
}

static void pack_sig_zlow(uint8_t *out, const poly *z1i)
{
    const uint32_t mask = (1u << T) - 1u;
    uint64_t acc = 0;
    int accbits = 0;
    size_t outpos = 0;
    for (int k = 0; k < N; k++) {
        int32_t head, low;
        pack_sig_split_z(z1i->coeffs[k], &head, &low);
        acc |= (uint64_t)(((uint32_t)low) & mask) << accbits;
        accbits += T;
        while (accbits >= 8) {
            out[outpos++] = (uint8_t)acc;
            acc >>= 8;
            accbits -= 8;
        }
    }
    if (accbits > 0)
        out[outpos] = (uint8_t)acc;
}

static void unpack_sig_zlow(poly *z1i, const uint8_t *in)
{
    const uint32_t mask = (1u << T) - 1u;
    const int32_t sgn = 1 << (T - 1);
    uint64_t acc = 0;
    int accbits = 0;
    size_t inpos = 0;
    for (int k = 0; k < N; k++) {
        while (accbits < T) {
            acc |= (uint64_t)in[inpos++] << accbits;
            accbits += 8;
        }
        int32_t low = (int32_t)((uint32_t)acc & mask);
        acc >>= T;
        accbits -= T;
        if (low & sgn)
            low -= 1 << T;
        z1i->coeffs[k] = z1i->coeffs[k] * KSPLIT + low;
    }
}

int pack_sig(uint8_t sig[CRYPTO_BYTES], const poly *c, const poly h[M],
             const poly z1[ELL + 1])
{
    uint8_t *p = sig;
    memset(p, 0, CHALLENGE_PACKEDBYTES);
    for (int k = 0; k < N; k++)
        if (c->coeffs[k])
            p[k >> 3] |= (uint8_t)(1u << (k & 7));
    p += CHALLENGE_PACKEDBYTES;

    uint8_t *lenp = p;
    p += 2;
    uint8_t *rp = p;
    size_t rlen;
    if (pack_sig_rans_encode(rp, &rlen, RANS_RESERVED_BYTES, h, z1) != 0)
        return -2; /* overflow: signer retries */
    /* pack_sig_rans_encode writes backward then memmoves to the front,
     * leaving stale bytes in [rlen, RESERVED); zero that padding so it is
     * authenticated (the verifier requires it to be zero -> no signature
     * malleability). */
    memset(rp + rlen, 0, RANS_RESERVED_BYTES - rlen);
    lenp[0] = (uint8_t)(rlen & 0xff);
    lenp[1] = (uint8_t)((rlen >> 8) & 0xff);
    p += RANS_RESERVED_BYTES;

    for (int i = 0; i < ELL + 1; i++) {
        pack_sig_zlow(p, &z1[i]);
        p += POLYZ_LO_PACKEDBYTES;
    }
    return 0;
}

int unpack_sig(poly *c, poly h[M], poly z1[ELL + 1],
               const uint8_t sig[CRYPTO_BYTES])
{
    const uint8_t *p = sig;
    for (int k = 0; k < N; k++)
        c->coeffs[k] = (p[k >> 3] >> (k & 7)) & 1;
    p += CHALLENGE_PACKEDBYTES;

    size_t rlen = (size_t)p[0] | ((size_t)p[1] << 8);
    p += 2;
    if (rlen > RANS_RESERVED_BYTES)
        return -1;
    /* The reserve tail beyond the actual rANS stream is zero-padding;
     * require it to be zero so every signature byte is authenticated (no
     * malleability). */
    for (size_t pad = rlen; pad < (size_t)RANS_RESERVED_BYTES; pad++)
        if (p[pad] != 0)
            return -1;
    if (unpack_sig_rans_decode(h, z1, p, rlen) != 0) {
        unpack_sig_clear_outputs(h, z1);
        return -1;
    }
    p += RANS_RESERVED_BYTES;

    for (int i = 0; i < ELL + 1; i++) {
        unpack_sig_zlow(&z1[i], p);
        p += POLYZ_LO_PACKEDBYTES;
    }
    return 0;
}
