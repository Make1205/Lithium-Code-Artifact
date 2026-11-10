#ifndef COPA_IRS64_H
#define COPA_IRS64_H
#include "copa_scalar64.h"
#include "copa_irs_limits.h"
typedef struct {
 uint64_t norm;
 /* Ephemeral coefficients: scalar canonical order, SIMD even/odd planes.
  * No context field is part of the serialized secret key. */
 uint64_t c[COPA64_IRS_TERMS] __attribute__((aligned(64)));
 uint64_t exp_gamma;
 copa64_ratio_ctx norm_ratio;
 copa_fp63 gamma,bump,gap;
} copa_irs64_context;
void copa_irs64_prepare(copa_irs64_context *ctx,uint64_t norm);
/* All arithmetic inside this step has a fixed schedule. Outcome probabilities
 * are Gartner f,g,1-f-g, evaluated with stable, equivalent factorizations. */
void copa_irs64_probabilities(copa_fp63 *abort,copa_fp63 *conditional_small,int64_t dot,const copa_irs64_context *ctx);
unsigned copa_irs64_step(unsigned *minus,int64_t dot,const copa_irs64_context *ctx,const uint64_t random[6]);
#endif
