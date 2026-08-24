/*
 * test_sigma2_diff.c - differential test: AVX-512 16-lane CDT kernels
 * (sampler_sigma2 over RCDT_Z, noise_magnitude_batch over
 * RCDT_NOISE) vs an independent in-test SCALAR reference (a
 * verbatim copy of the ref borrow chain), over random AND crafted corner
 * inputs. GAUSS_BATCH == NOISE_BATCH == 32 (= two 16-lane streams), so the
 * mini-batch CDT block is 4 groups x 96 = 384 bytes.
 *
 * The corner that matters is v_mid=0 with an incoming low-limb borrow:
 * there the borrow DOES propagate (the exact multi-precision compare), so
 * the AVX-512 kernel and the scalar must agree.  We force v0=v1=0 and pin
 * v2 to real Z2 table values to exercise the v2==Z2 equality boundary.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "params.h"
#include "sampler.h"

#define BATCH 32
#define CDT_BYTES (BATCH * 12) /* 384: 4 groups x 96 */

/* verbatim scalar reference grouped CDT scan (== ref/sampler.c body). */
static void ref_scan(int32_t *out, const uint8_t *rand,
                     const uint32_t Z[][3], int entries)
{
    for (int s = 0; s < BATCH; s++) {
        int group = s >> 3, lane = s & 7;
        int bb = group * 96 + lane * 4;
        uint32_t v0 = load_le32(rand + bb) & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + bb + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + bb + 64) & 0x7FFFFFFFU;
        int32_t z = 0;
        for (int i = 0; i < entries; i++) {
            uint32_t cc;
            cc = ct_lt_u32(v0, Z[i][0]);
            cc = ((uint32_t)(v1 - Z[i][1] - cc)) >> 31;
            cc = ((uint32_t)(v2 - Z[i][2] - cc)) >> 31;
            z += (int32_t)cc;
        }
        out[s] = z;
    }
}

static int cmp(const uint8_t *buf, int which, int *fails, int it,
               const char *tag)
{
    int32_t za[BATCH], zb[BATCH];
    if (which == 0) { /* sigma2 */
        ref_scan(za, buf, RCDT_Z, RCDT_ENTRIES);
        sampler_sigma2(zb, buf);
    } else { /* noise */
        ref_scan(za, buf, RCDT_NOISE, NOISE_RCDT_ENTRIES);
        noise_magnitude_batch(zb, buf);
    }
    if (memcmp(za, zb, sizeof za) != 0) {
        if (*fails < 8)
            for (int s = 0; s < BATCH; s++)
                if (za[s] != zb[s])
                    printf("  [%s] it=%d s=%d ref=%d avx512=%d\n", tag, it,
                           s, za[s], zb[s]);
        (*fails)++;
        return 1;
    }
    return 0;
}

static void run(int which, const char *name, int *fails)
{
    uint8_t buf[CDT_BYTES];
    const uint32_t(*Z)[3] = (which == 0) ? RCDT_Z : RCDT_NOISE;
    int entries = (which == 0) ? RCDT_ENTRIES : NOISE_RCDT_ENTRIES;

    /* (a) fully random */
    for (int it = 0; it < 150000; it++) {
        for (int i = 0; i < CDT_BYTES; i++)
            buf[i] = (uint8_t)rand();
        cmp(buf, which, fails, it, name);
    }
    /* (b) forced v0=v1=0 across all 4 groups (borrow into zero mid) */
    for (int it = 0; it < 80000; it++) {
        for (int i = 0; i < CDT_BYTES; i++)
            buf[i] = (uint8_t)rand();
        for (int g = 0; g < 4; g++) {
            memset(buf + g * 96, 0, 32);      /* v0 = 0 */
            memset(buf + g * 96 + 32, 0, 32); /* v1 = 0 */
        }
        cmp(buf, which, fails, it, name);
    }
    /* (c) v0=v1=0 AND v2 pinned to real Z2 table values (v2==Z2 boundary
     * under incoming borrow -- where a wrap chain would disagree). */
    for (int it = 0; it < entries; it++) {
        for (int i = 0; i < CDT_BYTES; i++)
            buf[i] = (uint8_t)rand();
        for (int g = 0; g < 4; g++) {
            memset(buf + g * 96, 0, 32);
            memset(buf + g * 96 + 32, 0, 32);
            for (int l = 0; l < 8; l++) {
                uint32_t z2 = Z[(it + l) % entries][2];
                uint8_t *v2 = buf + g * 96 + 64 + l * 4;
                v2[0] = (uint8_t)z2;
                v2[1] = (uint8_t)(z2 >> 8);
                v2[2] = (uint8_t)(z2 >> 16);
                v2[3] = (uint8_t)(z2 >> 24);
            }
        }
        cmp(buf, which, fails, it, name);
    }
}

int main(void)
{
    srand(20260603);
    int fails = 0;
    run(0, "sigma2", &fails);
    run(1, "noise", &fails);
    printf(
        "%s test_sigma2_diff mode=%d  fails=%d (sigma2+noise, 32x==1x)\n",
        fails ? "FAIL" : "PASS", LITHIUM_MODE, fails);
    return fails ? 1 : 0;
}
