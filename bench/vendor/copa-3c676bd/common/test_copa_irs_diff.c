/* Ref and vector IRS objects are compiled separately; compare actual
 * outputs. */
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>

#include "copa_irs64.h"
#include "copa_params.h"
void copa_ref_irs64_prepare(copa_irs64_context *, uint64_t);
void copa_ref_irs64_probabilities(copa_fp63 *, copa_fp63 *, int64_t,
                                  const copa_irs64_context *);
int main(void)
{
    unsigned points = 0;
    for (unsigned norm = COPA_MIN_NORM; norm <= COPA_BK_SQ; ++norm) {
        copa_irs64_context ref, simd;
        copa_ref_irs64_prepare(&ref, norm);
        copa_irs64_prepare(&simd, norm);
        int64_t dots[] = {0,
                          1,
                          -1,
                          -(int64_t)norm - 1,
                          -(int64_t)norm,
                          1 - (int64_t)norm,
                          (int64_t)norm - 1,
                          norm,
                          (int64_t)norm + 1,
                          -1000000,
                          1000000};
        for (unsigned j = 0; j < sizeof(dots) / sizeof(dots[0]); ++j) {
            copa_fp63 a, b, c, d;
            copa_ref_irs64_probabilities(&a, &b, dots[j], &ref);
            copa_irs64_probabilities(&c, &d, dots[j], &simd);
            assert(a.m == c.m && a.e == c.e && b.m == d.m && b.e == d.e);
            ++points;
        }
    }
    printf(
        "{\"mode\":%d,\"all_norms\":%u,\"points\":%u,\"bitwise_equal\":"
        "true}\n",
        COPA_MODE, COPA_BK_SQ, points);
    return 0;
}
