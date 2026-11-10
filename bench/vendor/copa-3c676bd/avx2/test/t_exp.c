/* t_exp.c - C approx_exp vs math.exp; matches tools/gen_approx_exp.py
 * gates: approx_exp      (Q60, sampler scheme): relative error over [0,8]
 * < 2^-44 (achieves ~2^-49 with the B=8 degree-6 minimax + lookup),
 *   approx_exp_wide (Q52, IRS table-reduced scheme): absolute error over
 * [0,40] < 2^-50. (Beyond a~30 exp(-a) nears the Q63 floor 2^-63, so
 * relative bits there are intentionally not gated; those terms are summed
 * negligibly in irs.c.)
 *   exp_small_signed (Q52/Q62): relative error over the actual IRS domain
 * [0, max BK_SQ/r^2] < 2^-44. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "approx_exp.h"
#include "params.h"

int main(void)
{
    const int niter = 200000;
    double rel8 = 0.0, abs40 = 0.0, rels = 0.0;
    /* sampler routine: Q60 input, relative error over the accept range
     * [0,8] */
    for (int i = 0; i <= niter; i++) {
        double a = 8.0 * i / niter;
        uint64_t q = (uint64_t)(a * (double)(1ULL << 60) + 0.5);
        double got = (double)approx_exp(q) / (double)(1ULL << 63);
        double want = exp(-a);
        double dr = fabs(got - want) / want;
        if (dr > rel8)
            rel8 = dr;
    }
    /* wide routine: Q52 input, absolute error over the irs.c range [0,40]
     */
    for (int i = 0; i <= niter; i++) {
        double a = 40.0 * i / niter;
        uint64_t q = (uint64_t)(a * (double)(1ULL << 52) + 0.5);
        double got = (double)approx_exp_wide(q) / (double)(1ULL << 63);
        double da = fabs(got - exp(-a));
        if (da > abs40)
            abs40 = da;
    }
    /* IRS clamps a_c to 2*gamma, whose maximum over accepted keys is
     * BK_SQ/r^2.  This used to exceed the old documented small-exp gate.
     */
    const double dmax = (double)BK_SQ / ((double)R * (double)R);
    for (int i = 0; i <= niter; i++) {
        double d = dmax * i / niter;
        int64_t q = (int64_t)(d * (double)(1ULL << 52) + 0.5);
        double got = (double)exp_small_signed(q) / (double)(1ULL << 62);
        double want = exp(d);
        double dr = fabs(got - want) / want;
        if (dr > rels)
            rels = dr;
    }
    int ok = (rel8 < 5.7e-14) && (abs40 < 8.9e-16) && (rels < 5.7e-14);
    printf(
        "%s t_exp mode=%d rel[0,8]=%.3e abs[0,40]=%.3e "
        "small[0,%.5f]=%.3e\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, rel8, abs40, dmax, rels);
    return ok ? 0 : 1;
}
