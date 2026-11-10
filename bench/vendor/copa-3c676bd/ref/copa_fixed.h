#ifndef COPA_FIXED_H
#define COPA_FIXED_H
#include <stdint.h>
#include "copa_params.h"
typedef __uint128_t copa_u128;
typedef __int128 copa_i128;
#include "copa_fixed_tables.h"
copa_i128 copa_mul_fixed(copa_i128 a,copa_i128 b);
/* floor((num/den)*2^120), capped at 88*2^120. den is public. */
copa_u128 copa_ratio_fixed(uint64_t num,uint32_t den);
copa_u128 copa_exp_negative(copa_u128 a);
/* Caller proves 0 <= a <= 1. */
copa_u128 copa_exp_positive(copa_u128 a);
#endif
