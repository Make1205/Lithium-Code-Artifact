#ifndef COPA_Q120_SAMPLER_H
#define COPA_Q120_SAMPLER_H
#include "copa_params.h"
#include "copa_fixed.h"
void copa_expand_y(int32_t y[COPA_K][COPA_N],const uint8_t seed[32],uint64_t attempt);
void copa_irs_thresholds(copa_i128 *minus,copa_i128 *accept,int64_t dot,uint64_t norm_sq,const copa_i128 table[COPA_IRS_TERMS]);
void copa_irs_table(copa_i128 table[COPA_IRS_TERMS],uint64_t norm_sq);
/* Precondition: 1 <= norm_sq <= COPA_BK_SQ and norm_sq is the squared
 * norm of s. Fixed sequence for secret data. Challenge and its support are public.
 * Every challenge step runs even after rejection, then the aggregate
 * acceptance bit is released to the outer signing retry loop. */
int copa_irs(int32_t z[COPA_K][COPA_N],const uint8_t challenge[COPA_N],
             const int32_t s[COPA_K][COPA_N],uint64_t norm_sq,
             const copa_i128 table[COPA_IRS_TERMS],const uint8_t seed[32],uint64_t attempt);
#endif
