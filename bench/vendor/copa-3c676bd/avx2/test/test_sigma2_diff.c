/*
 * test_sigma2_diff.c - differential test: AVX2 8-way sampler_sigma2 vs an
 * independent in-test SCALAR reference (a verbatim copy of the ref borrow
 * chain), over random AND crafted corner inputs.
 *
 * The chain is the exact 31-bit-limb borrow: borrow0 = (v0 <u Z0), then
 * borrow_k = bit 31 of (v_k - Z_k - borrow_{k-1}).  The corner that
 * matters is v_mid=0 with an incoming borrow: there the borrow DOES
 * propagate (the exact 93-bit compare), so the AVX2 kernel and the scalar
 * must agree.  We force v0=v1=0 (guaranteed low-limb borrow into a zero
 * mid) and also pin v2 to real Z-table values to exercise the v2==Z2
 * equality boundary.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "params.h"
#include "sampler.h"

/* verbatim scalar reference (== ref/sampler.c sampler_sigma2 body). */
static void ref_sigma2(int32_t *z_out, const uint8_t *rand)
{
    for (int s = 0; s < GAUSS_BATCH; s++) {
        int group = s >> 3, lane = s & 7;
        int bb = group * 96 + lane * 4;
        uint32_t v0 = load_le32(rand + bb) & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + bb + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + bb + 64) & 0x7FFFFFFFU;
        int32_t z = 0;
        for (int i = 0; i < RCDT_ENTRIES; i++) {
            uint32_t cc;
            cc = ct_lt_u32(v0, RCDT_Z[i][0]);
            cc = ((uint32_t)(v1 - RCDT_Z[i][1] - cc)) >> 31;
            cc = ((uint32_t)(v2 - RCDT_Z[i][2] - cc)) >> 31;
            z += (int32_t)cc;
        }
        z_out[s] = z;
    }
}

static int cmp_batch(const uint8_t *buf, int *fails, int it,
                     const char *tag)
{
    int32_t za[GAUSS_BATCH], zb[GAUSS_BATCH];
    ref_sigma2(za, buf);
    sampler_sigma2(zb, buf);
    if (memcmp(za, zb, sizeof za) != 0) {
        if (*fails < 8)
            for (int s = 0; s < GAUSS_BATCH; s++)
                if (za[s] != zb[s])
                    printf("  [%s] it=%d s=%d ref=%d avx=%d\n", tag, it, s,
                           za[s], zb[s]);
        (*fails)++;
        return 1;
    }
    return 0;
}

int main(void)
{
    srand(20260602);
    int fails = 0;
    uint8_t buf[192];

    /* (a) fully random */
    for (int it = 0; it < 200000; it++) {
        for (int i = 0; i < 192; i++)
            buf[i] = (uint8_t)rand();
        cmp_batch(buf, &fails, it, "rand");
    }

    /* (b) forced corner: zero ALL mid limbs (v1=0) and ALL low limbs
     * (v0=0) so the low-limb compare borrows into a zero mid limb for
     * every entry with Z0>0 -- the borrow-into-zero-mid propagation case.
     * v2 random. */
    for (int it = 0; it < 100000; it++) {
        for (int i = 0; i < 192; i++)
            buf[i] = (uint8_t)rand();
        memset(buf + 0, 0, 32);   /* group0 low  (v0 = 0) */
        memset(buf + 32, 0, 32);  /* group0 mid  (v1 = 0) */
        memset(buf + 96, 0, 32);  /* group1 low  (v0 = 0) */
        memset(buf + 128, 0, 32); /* group1 mid  (v1 = 0) */
        cmp_batch(buf, &fails, it, "v1=0");
    }

    /* (d) v0=v1=0 (borrow into zero mid) AND v2 pinned to a real Z2 table
     * value, so the high limb hits the v2==Z2 equality boundary under an
     * incoming borrow -- the exact spot where a wrap chain (no propagate)
     * and the exact chain (propagate) disagree.  Pins down that the AVX2
     * kernel implements the same exact borrow as the scalar reference. */
    for (int it = 0; it < RCDT_ENTRIES; it++) {
        for (int i = 0; i < 192; i++)
            buf[i] = (uint8_t)rand();
        memset(buf + 0, 0, 32);
        memset(buf + 32, 0, 32);
        memset(buf + 96, 0, 32);
        memset(buf + 128, 0, 32);
        for (int g = 0; g < 2; g++)
            for (int l = 0; l < 8; l++) {
                uint32_t z2 = RCDT_Z[(it + l) % RCDT_ENTRIES][2];
                uint8_t *v2 = buf + g * 96 + 64 + l * 4;
                v2[0] = (uint8_t)z2;
                v2[1] = (uint8_t)(z2 >> 8);
                v2[2] = (uint8_t)(z2 >> 16);
                v2[3] = (uint8_t)(z2 >> 24);
            }
        cmp_batch(buf, &fails, it, "v2==Z2");
    }

    /* (c) mid limbs near boundaries (0,1) with random low/high to vary
     * borrow. */
    for (int it = 0; it < 100000; it++) {
        for (int i = 0; i < 192; i++)
            buf[i] = (uint8_t)rand();
        for (int g = 0; g < 2; g++)
            for (int l = 0; l < 8; l++) {
                uint8_t *mid = buf + g * 96 + 32 + l * 4;
                uint32_t mv = (uint32_t)(rand() & 3); /* 0,1,2,3 */
                mid[0] = (uint8_t)mv;
                mid[1] = mid[2] = mid[3] = 0;
            }
        cmp_batch(buf, &fails, it, "mid<4");
    }

    printf("%s test_sigma2_diff mode=%d  fails=%d (random+corner)\n",
           fails ? "FAIL" : "PASS", LITHIUM_MODE, fails);
    return fails ? 1 : 0;
}
