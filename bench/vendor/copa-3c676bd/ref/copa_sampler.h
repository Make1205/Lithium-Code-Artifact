#ifndef COPA_SAMPLER_H
#define COPA_SAMPLER_H
#include "copa_params.h"
#include "copa_irs64.h"
#include "copa_keygen.h"
void copa_sample_batch(uint32_t candidate[32],uint32_t accepted[32],const uint64_t u[][32],const uint8_t low[32],const uint8_t tail[256]);
void copa_expand_y(int32_t y[COPA_K][COPA_N],const uint8_t seed[32],uint64_t attempt);
#if defined(COPA_PARALLEL_XOF) && (COPA_MODE==2 || COPA_MODE==3 || defined(COPA_ENABLE_IRS_POOL))
#define COPA_USE_IRS_POOL 1
#endif
#ifdef COPA_USE_IRS_POOL
#include "copa_xof8.h"
#if COPA_SIMD_BITS==256
#define COPA_IRS_POOL_LANES 4
#else
#define COPA_IRS_POOL_LANES 8
#endif
#define COPA_IRS_POOL_BYTES (((48*COPA_TAU+135)/136)*136)
typedef struct{copa_xof8 engine;uint64_t group;unsigned ready,filled;uint8_t data[COPA_IRS_POOL_LANES][COPA_IRS_POOL_BYTES];} copa_irs_pool;
int copa_irs_with_pool(int32_t z[COPA_K][COPA_N],const uint8_t challenge[COPA_N],const copa_secret_key *sk,const uint32_t z_ntt[COPA_K][2][COPA_N],const uint8_t seed[32],uint64_t attempt,unsigned *steps,copa_irs_pool *pool);
#endif
int copa_irs(int32_t z[COPA_K][COPA_N],const uint8_t challenge[COPA_N],
             const copa_secret_key *sk,const uint32_t z_ntt[COPA_K][2][COPA_N],
             const uint8_t seed[32],uint64_t attempt,unsigned *steps);
#endif
