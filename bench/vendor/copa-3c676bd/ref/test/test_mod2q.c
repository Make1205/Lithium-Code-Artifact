/* test_mod2q.c - gate for rounding.c (mod-2q HighBits / hint).
 *
 *  (1) HighBits/LowBits: |LowBits| <= tau/2 + 2 (Lemma 4).
 *  (2) MakeHint/UseHint: usehint(wtilde, makehint(w, wtilde)) ==
 * highbits(w) for wtilde = (w - 2 z2) mod 2q with bounded z2. (3) tau |
 * 2(q-1); W1_MAX = 2(q-1)/tau.
 */
#include <stdint.h>
#include <stdio.h>

#include "params.h"
#include "rounding.h"

static uint64_t st = 0x9E3779B97F4A7C15ULL;
static uint32_t xr(void)
{
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    return (uint32_t)(st >> 27);
}

int main(void)
{
    int fails = 0, lowbits_max = 0;

    if ((2 * (Q - 1)) % TAU != 0) {
        printf("FAIL test_mod2q: tau!|2(q-1)\n");
        return 1;
    }
    if (W1_MAX != (2 * (Q - 1)) / TAU) {
        printf("FAIL test_mod2q: W1_MAX\n");
        return 1;
    }

    for (int t = 0; t < 200000; t++) {
        int32_t w = (int32_t)(xr() % DQ);

        int32_t b = highbits(w);
        int32_t hv = hbvalue(b);
        int32_t low = reduce2q(w) - hv;
        int32_t mod2 = 2 * (Q - 1);
        low %= mod2;
        if (low < 0)
            low += mod2;
        if (low > mod2 / 2)
            low -= mod2;
        int a = low < 0 ? -low : low;
        if (a > lowbits_max)
            lowbits_max = a;
        if (a > TAU / 2 + 2) {
            fails++;
            if (fails <= 3)
                printf("  lowbits %d > tau/2+2 (w=%d)\n", a, w);
        }

        int32_t z2 = (int32_t)(xr() % 6001) - 3000;
        int32_t wt = reduce2q(w - 2 * z2);
        int32_t h = makehint(w, wt);
        int32_t rec = usehint(wt, h);
        if (rec != b) {
            fails++;
            if (fails <= 3)
                printf("  hint mismatch w=%d z2=%d h=%d rec=%d b=%d\n", w,
                       z2, h, rec, b);
        }
    }

    int ok = (fails == 0);
    printf(
        "%s test_mod2q mode=%d  W1_MAX=%d tau=%d  lowbits_max=%d "
        "(<=tau/2+2=%d)  fails=%d\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, W1_MAX, TAU, lowbits_max,
        TAU / 2 + 2, fails);
    return ok ? 0 : 1;
}
