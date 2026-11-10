/* Independent scalar comparisons for the SIMD partial products and
 * rotations. */
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "copa_fp_simd.h"
#include "copa_poly_simd.h"
static uint64_t state = UINT64_C(0x625482d4ac37ba19);
static uint64_t next(void)
{
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}
static uint64_t q63(uint64_t a, uint64_t b)
{
    uint64_t lo, hi;
    __asm__("mulq %3" : "=a"(lo), "=d"(hi) : "a"(a), "r"(b) : "cc");
    return (hi << 1 | lo >> 63) + ((lo >> 62) & 1);
}
int main(void)
{
    uint64_t a[8], b[8], out[8];
    const uint64_t edge[] = {0,
                             1,
                             UINT64_MAX,
                             UINT64_C(1) << 63,
                             (UINT64_C(1) << 63) - 1,
                             UINT64_C(0xffffffff),
                             UINT64_C(0x100000000)};
    for (unsigned t = 0; t < 100000; ++t) {
        for (unsigned i = 0; i < 8; ++i) {
            a[i] = t < 49 ? edge[t / 7] : next();
            b[i] = t < 49 ? edge[t % 7] : next();
        }
        copa_backend_mulq63x8(out, a, b);
        for (unsigned i = 0; i < 8; ++i)
            assert(out[i] == q63(a[i], b[i]));
        copa_backend_mulq63x4(out, a, b);
        for (unsigned i = 0; i < 4; ++i)
            assert(out[i] == q63(a[i], b[i]));
    }
    int32_t x[COPA_N], y[COPA_N], s[COPA_N];
    for (unsigned i = 0; i < COPA_N; ++i)
        s[i] = (int32_t)(next() % 401) - 200;
    for (unsigned shift = 0; shift < COPA_N; ++shift)
        for (int minus = 0; minus < 2; ++minus) {
            uint64_t norm = 0;
            for (unsigned i = 0; i < COPA_N; ++i)
                x[i] = y[i] = (int32_t)(next() % 20001) - 10000;
            int sign = minus ? -1 : 1;
            for (unsigned i = 0; i < COPA_N; ++i) {
                unsigned j = (i + shift) % COPA_N;
                y[j] += sign * s[i] * (i + shift < COPA_N ? 1 : -1);
            }
            copa_backend_rotate_add(x, s, shift, -minus);
            assert(!memcmp(x, y, sizeof(x)));
            for (unsigned i = 0; i < COPA_N; ++i)
                norm += (uint64_t)((int64_t)x[i] * x[i]);
            assert(copa_backend_norm2(x, COPA_N) == norm);
        }
    for(unsigned t=0;t<1000;++t){
        unsigned positions[COPA_TAU];int32_t masks[COPA_TAU];
        for(unsigned j=0;j<COPA_N;++j)x[j]=y[j]=(int32_t)(next()%20001)-10000;
        for(unsigned k=0;k<COPA_TAU;++k){
            unsigned shift=positions[k]=(unsigned)(next()%COPA_N);
            int minus=(int)(next()&1);masks[k]=-minus;
            for(unsigned j=0;j<COPA_N;++j){unsigned at=(j+shift)%COPA_N;
                y[at]+=(minus?-1:1)*s[j]*(j+shift<COPA_N?1:-1);}
        }
        copa_backend_sparse_add(x,s,positions,masks,COPA_TAU);
        assert(!memcmp(x,y,sizeof(x)));
    }
    puts(
        "PASS SIMD: 1.2M rounded full-width Q63 products; all 512 signed "
        "negacyclic rotations and norms");
    return 0;
}
