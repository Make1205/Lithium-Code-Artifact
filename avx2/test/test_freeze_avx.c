/* test_freeze_avx.c - differential tests for AVX2 NTT input conversion. */
#include <stdint.h>
#include <stdio.h>

#include "ntt_avx_decls.h"
#include "params.h"
#include "reduce.h"

static uint32_t rng = 1;
static uint32_t xr(void)
{
    rng = 1664525u * rng + 1013904223u;
    return rng;
}

static uint16_t small_ref(int32_t x)
{
    uint32_t ux = (uint32_t)x;
    uint32_t mask = 0u - (ux >> 31);
    return (uint16_t)(ux + (mask & (uint32_t)Q));
}

static int check_case(const int32_t in[N])
{
    int16_t out[N];
    poly_freeze_to_u16_avx(out, in);
    for (int i = 0; i < N; i++) {
        uint16_t want = (uint16_t)freeze(in[i]);
        if ((uint16_t)out[i] != want) {
            printf("FAIL full mode=%d i=%d in=%d got=%u want=%u\n",
                   LITHIUM_MODE, i, in[i], (uint16_t)out[i], want);
            return 1;
        }
    }
    return 0;
}

int main(void)
{
    int32_t in[N];
    int16_t out[N];
    const int32_t corners[] = {0,
                               1,
                               -1,
                               Q,
                               -Q,
                               Q - 1,
                               -(Q - 1),
                               2 * Q,
                               -2 * Q,
                               INT32_MAX,
                               INT32_MIN,
                               INT32_MAX - Q,
                               INT32_MIN + Q};

    for (int rep = 0; rep < 200; rep++) {
        for (int i = 0; i < N; i++)
            in[i] = (int32_t)xr();
        if (check_case(in))
            return 1;
        for (int i = 0; i < N; i++)
            in[i] = (int32_t)(xr() & 0x7fffffffU);
        if (check_case(in))
            return 1;
        for (int i = 0; i < N; i++)
            in[i] = -(int32_t)(xr() & 0x7fffffffU);
        if (check_case(in))
            return 1;
    }
    for (unsigned c = 0; c < sizeof corners / sizeof corners[0]; c++) {
        for (int i = 0; i < N; i++)
            in[i] = corners[(c + (unsigned)i) %
                            (sizeof corners / sizeof corners[0])];
        if (check_case(in))
            return 1;
    }

    for (int i = 0; i < N; i++)
        in[i] = (int32_t)(xr() % (2 * Q - 1)) - (Q - 1);
    poly_freeze_small_to_u16_avx(out, in);
    for (int i = 0; i < N; i++) {
        uint16_t want = small_ref(in[i]);
        if ((uint16_t)out[i] != want) {
            printf("FAIL small mode=%d i=%d in=%d got=%u want=%u\n",
                   LITHIUM_MODE, i, in[i], (uint16_t)out[i], want);
            return 1;
        }
    }

    for (int i = 0; i < N; i++)
        in[i] = (int32_t)(xr() % Q);
    poly_cast_modq_avx(out, in);
    for (int i = 0; i < N; i++) {
        if ((uint16_t)out[i] != (uint16_t)in[i]) {
            printf("FAIL cast mode=%d i=%d in=%d got=%u\n", LITHIUM_MODE,
                   i, in[i], (uint16_t)out[i]);
            return 1;
        }
    }

    printf("PASS test_freeze_avx mode=%d n=%d q=%d\n", LITHIUM_MODE, N, Q);
    return 0;
}
