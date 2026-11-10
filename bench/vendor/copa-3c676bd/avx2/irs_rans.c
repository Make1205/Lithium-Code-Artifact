/*
 * irs_rans.c - byte-wise rANS encoder/decoder for IRS (two-table,
 * RANS_N-way interleaved stream).
 *
 * Matches tools/rans.py byte-for-byte. 32-bit state, L=2^23, 8-bit renorm,
 * prob_bits=10.  The flat symbol sequence is z-head (table ZH) followed by
 * hint (table HINT); symbol t is coded on its power-of-two interleave
 * state. Encode pushes in reverse t (so decode pulls forward t), all
 * states sharing one backward-written byte stream.  RANS_N independent
 * states break the per-symbol x dependency chain -- each enc_put is a
 * reciprocal-multiply + adds dependent on its state's x, so a single
 * stream is latency-bound on that ~5-cyc multiply.  RANS_N=2 already hides
 * it (isolated encode bench, mode120/1280 symbols: N=1 10080 cyc -> N=2
 * 7304 (1.38x); N>=4 no further gain, just more flush bytes), so N=2 is
 * the size/speed optimum: it costs only 4 extra flush bytes, negligible
 * against RANS_RESERVED_BYTES.
 */
#include "irs_rans.h"

#include <stdint.h>
#include <string.h>

#define PB RANS_PROB_BITS /* 10 */
#define PSCALE (1u << PB)
#define RANS_N RANS_INTERLEAVED_STREAMS
#define RANS_MASK RANS_INTERLEAVED_MASK

static inline unsigned rans_stream_index(size_t t)
{
    return (unsigned)(t & RANS_MASK);
}

static inline int rans_state_in_range(uint32_t x)
{
    return x >= RANS_L && x < (RANS_L << 8);
}

/* one encode step (push): updates *x, writes renorm bytes to out[--*pos].
 * The quotient and remainder are replaced by a multiply-by-reciprocal
 * (q is computed from x and rcp; x = x + bias +
 * q*(PSCALE-freq)), which is bit-exact for the encoder's x range (ryg rANS
 * / Granlund-Montgomery; rcp/rsh/bias are generated + verified in
 * tools/gen_rans_tables.py).  No divide on the (public) signature data;
 * the byte stream -- hence the KAT -- is unchanged. */
static inline void enc_put(uint32_t *x, uint8_t *out, size_t *pos,
                           uint32_t freq, uint32_t rcp, uint32_t rsh,
                           uint32_t bias)
{
    uint32_t x_max = ((RANS_L >> PB) << 8) * freq;
    while (*x >= x_max) {
        out[--(*pos)] = (uint8_t)(*x & 0xff);
        *x >>= 8;
    }
    uint32_t q = (uint32_t)(((uint64_t)*x * rcp) >> 32) >> rsh;
    *x = *x + bias + q * (PSCALE - freq);
}

int rans_encode(uint8_t *out, size_t *out_len, size_t max_bytes,
                const int32_t *z_head, size_t nz, const int32_t *hint,
                size_t nh)
{
    uint32_t x[RANS_N];
    for (int s = 0; s < RANS_N; s++)
        x[s] = RANS_L;
    size_t pos = max_bytes; /* write backwards */
    size_t ntot = nz + nh;

    /* flat S[t]: z_head[t] (t<nz, ZH) then hint[t-nz] (HINT).  Push in
     * REVERSE t; decode mirrors forward t with the same state schedule. */
    for (size_t t = ntot; t-- > 0;) {
        unsigned s = rans_stream_index(t);
        if (pos < (size_t)(4 * RANS_N + 4))
            return -2;
        if (t < nz) {
            int slot = (int)(z_head[t] - RANS_ZH_LO);
            if (slot < 0 || slot >= RANS_ZH_N)
                return -2;
            enc_put(&x[s], out, &pos, RANS_ZH_FREQ[slot],
                    RANS_ZH_RCP[slot], RANS_ZH_RSH[slot],
                    RANS_ZH_BIAS[slot]);
        } else {
            int slot = (int)(hint[t - nz] - RANS_HINT_LO);
            if (slot < 0 || slot >= RANS_HINT_N)
                return -2;
            enc_put(&x[s], out, &pos, RANS_HINT_FREQ[slot],
                    RANS_HINT_RCP[slot], RANS_HINT_RSH[slot],
                    RANS_HINT_BIAS[slot]);
        }
    }
    /* flush all RANS_N states, state 0 first (lands last in the stream, so
     * decode reads state RANS_N-1 first). */
    for (int s = 0; s < RANS_N; s++)
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

int rans_decode(int32_t *z_head, size_t nz, int32_t *hint, size_t nh,
                const uint8_t *in, size_t in_len)
{
    size_t ntot = nz + nh;
    if (in_len < (size_t)(4 * RANS_N))
        return -1;
    uint32_t x[RANS_N];
    size_t bp = 0;
    /* init states in reverse flush order: the stream front is state
     * RANS_N-1, ..., then state 0 (mirror of the encoder's flush). */
    for (int s = RANS_N - 1; s >= 0; s--) {
        uint32_t v = 0;
        for (int b = 0; b < 4; b++)
            v = (v << 8) | in[bp++];
        x[s] = v;
        if (!rans_state_in_range(x[s]))
            return -1;
    }

    for (size_t t = 0; t < ntot; t++) {
        unsigned s = rans_stream_index(t);
        /* val in [0, PSCALE) -> slot via O(1) direct lookup (SLOT table,
         * gen_rans_tables.py); CDF[0]=0, CDF[N]=PSCALE so every val maps
         * to a valid slot (no miss).  Avoids an O(n) linear CDF
         * scan. */
        uint32_t val = x[s] & (PSCALE - 1);
        if (t < nz) {
            unsigned slot = RANS_ZH_SLOT[val];
            z_head[t] = RANS_ZH_LO + (int)slot;
            x[s] = RANS_ZH_FREQ[slot] * (x[s] >> PB) + val -
                   RANS_ZH_CDF[slot];
        } else {
            unsigned slot = RANS_HINT_SLOT[val];
            hint[t - nz] = RANS_HINT_LO + (int)slot;
            x[s] = RANS_HINT_FREQ[slot] * (x[s] >> PB) + val -
                   RANS_HINT_CDF[slot];
        }
        while (x[s] < RANS_L) {
            if (bp >= in_len)
                return -1;
            x[s] = (x[s] << 8) | in[bp++];
        }
    }
    /* Canonical stream check: all bytes must be consumed and all
     * interleaved states must rewind to the encoder's initial state. */
    if (bp != in_len)
        return -1;
    for (int s = 0; s < RANS_N; s++)
        if (x[s] != RANS_L)
            return -1;
    return 0;
}
