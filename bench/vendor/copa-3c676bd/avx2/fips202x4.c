#include "fips202x4.h"

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "fips202.h"
#include "prof.h"

static void keccakx4_absorb_once(__m256i s[25], unsigned int r,
                                 const uint8_t *in0, const uint8_t *in1,
                                 const uint8_t *in2, const uint8_t *in3,
                                 size_t inlen, uint8_t p)
{
    size_t i;
    uint64_t pos = 0;
    __m256i t, idx;

    for (i = 0; i < 25; ++i)
        s[i] = _mm256_setzero_si256();

    idx = _mm256_set_epi64x((long long)in3, (long long)in2, (long long)in1,
                            (long long)in0);
    while (inlen >= r) {
        for (i = 0; i < r / 8; ++i) {
            t = _mm256_i64gather_epi64((long long *)pos, idx, 1);
            s[i] = _mm256_xor_si256(s[i], t);
            pos += 8;
        }
        inlen -= r;

        f1600x4(s, KeccakF_RoundConstants);
    }

    for (i = 0; i < inlen / 8; ++i) {
        t = _mm256_i64gather_epi64((long long *)pos, idx, 1);
        s[i] = _mm256_xor_si256(s[i], t);
        pos += 8;
    }
    inlen -= 8 * i;

    if (inlen) {
        uint64_t tail[4] = {0, 0, 0, 0};
        memcpy(&tail[0], in0 + pos, inlen);
        memcpy(&tail[1], in1 + pos, inlen);
        memcpy(&tail[2], in2 + pos, inlen);
        memcpy(&tail[3], in3 + pos, inlen);
        t = _mm256_loadu_si256((const __m256i *)tail);
        s[i] = _mm256_xor_si256(s[i], t);
    }

    t = _mm256_set1_epi64x((uint64_t)p << 8 * inlen);
    s[i] = _mm256_xor_si256(s[i], t);
    t = _mm256_set1_epi64x(1ULL << 63);
    s[r / 8 - 1] = _mm256_xor_si256(s[r / 8 - 1], t);
}

static void keccakx4_squeezeblocks(uint8_t *out0, uint8_t *out1,
                                   uint8_t *out2, uint8_t *out3,
                                   size_t nblocks, unsigned int r,
                                   __m256i s[25])
{
    unsigned int i;
    __m128d t;
    PROF_SQ(4ULL * (uint64_t)nblocks *
            r); /* 4 lanes; no-op unless PROF_RAND */

    while (nblocks > 0) {
        f1600x4(s, KeccakF_RoundConstants);
        for (i = 0; i < r / 8; ++i) {
            t = _mm_castsi128_pd(_mm256_castsi256_si128(s[i]));
            _mm_storeu_si64(&out0[8 * i], _mm_castpd_si128(t));
            _mm_storeu_si64(&out1[8 * i], _mm_srli_si128(_mm_castpd_si128(t), 8));
            t = _mm_castsi128_pd(_mm256_extracti128_si256(s[i], 1));
            _mm_storeu_si64(&out2[8 * i], _mm_castpd_si128(t));
            _mm_storeu_si64(&out3[8 * i], _mm_srli_si128(_mm_castpd_si128(t), 8));
        }

        out0 += r;
        out1 += r;
        out2 += r;
        out3 += r;
        --nblocks;
    }
}

void shake128x4_absorb_once(keccakx4_state *state, const uint8_t *in0,
                            const uint8_t *in1, const uint8_t *in2,
                            const uint8_t *in3, size_t inlen)
{
    keccakx4_absorb_once(state->s, SHAKE128_RATE, in0, in1, in2, in3,
                         inlen, 0x1F);
}

void shake128x4_squeezeblocks(uint8_t *out0, uint8_t *out1, uint8_t *out2,
                              uint8_t *out3, size_t nblocks,
                              keccakx4_state *state)
{
    keccakx4_squeezeblocks(out0, out1, out2, out3, nblocks, SHAKE128_RATE,
                           state->s);
}

void shake256x4_absorb_once(keccakx4_state *state, const uint8_t *in0,
                            const uint8_t *in1, const uint8_t *in2,
                            const uint8_t *in3, size_t inlen)
{
    keccakx4_absorb_once(state->s, SHAKE256_RATE, in0, in1, in2, in3,
                         inlen, 0x1F);
}

void shake256x4_squeezeblocks(uint8_t *out0, uint8_t *out1, uint8_t *out2,
                              uint8_t *out3, size_t nblocks,
                              keccakx4_state *state)
{
    keccakx4_squeezeblocks(out0, out1, out2, out3, nblocks, SHAKE256_RATE,
                           state->s);
}

void shake128x4(uint8_t *out0, uint8_t *out1, uint8_t *out2, uint8_t *out3,
                size_t outlen, const uint8_t *in0, const uint8_t *in1,
                const uint8_t *in2, const uint8_t *in3, size_t inlen)
{
    unsigned int i;
    size_t nblocks = outlen / SHAKE128_RATE;
    uint8_t t[4][SHAKE128_RATE];
    keccakx4_state state;

    shake128x4_absorb_once(&state, in0, in1, in2, in3, inlen);
    shake128x4_squeezeblocks(out0, out1, out2, out3, nblocks, &state);

    out0 += nblocks * SHAKE128_RATE;
    out1 += nblocks * SHAKE128_RATE;
    out2 += nblocks * SHAKE128_RATE;
    out3 += nblocks * SHAKE128_RATE;
    outlen -= nblocks * SHAKE128_RATE;

    if (outlen) {
        shake128x4_squeezeblocks(t[0], t[1], t[2], t[3], 1, &state);
        for (i = 0; i < outlen; ++i) {
            out0[i] = t[0][i];
            out1[i] = t[1][i];
            out2[i] = t[2][i];
            out3[i] = t[3][i];
        }
    }
}

void shake256x4(uint8_t *out0, uint8_t *out1, uint8_t *out2, uint8_t *out3,
                size_t outlen, const uint8_t *in0, const uint8_t *in1,
                const uint8_t *in2, const uint8_t *in3, size_t inlen)
{
    unsigned int i;
    size_t nblocks = outlen / SHAKE256_RATE;
    uint8_t t[4][SHAKE256_RATE];
    keccakx4_state state;

    shake256x4_absorb_once(&state, in0, in1, in2, in3, inlen);
    shake256x4_squeezeblocks(out0, out1, out2, out3, nblocks, &state);

    out0 += nblocks * SHAKE256_RATE;
    out1 += nblocks * SHAKE256_RATE;
    out2 += nblocks * SHAKE256_RATE;
    out3 += nblocks * SHAKE256_RATE;
    outlen -= nblocks * SHAKE256_RATE;

    if (outlen) {
        shake256x4_squeezeblocks(t[0], t[1], t[2], t[3], 1, &state);
        for (i = 0; i < outlen; ++i) {
            out0[i] = t[0][i];
            out1[i] = t[1][i];
            out2[i] = t[2][i];
            out3[i] = t[3][i];
        }
    }
}

/* Copa's aligned opaque engine: vector lvalues carry may_alias, avoiding
 * both aggregate type punning and a 3200-byte state copy per eight blocks. */
void copa_keccak4_absorb(void *state,const uint8_t *a,const uint8_t *b,const uint8_t *c,const uint8_t *d,size_t n){keccakx4_absorb_once((__m256i *)state,SHAKE256_RATE,a,b,c,d,n,0x1f);}
void copa_keccak4_block(void *state,uint8_t *a,uint8_t *b,uint8_t *c,uint8_t *d){keccakx4_squeezeblocks(a,b,c,d,1,SHAKE256_RATE,(__m256i *)state);}
