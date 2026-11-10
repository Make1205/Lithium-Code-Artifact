#ifndef COPA_SAMPLER_TABLES_H
#define COPA_SAMPLER_TABLES_H
#include <stdint.h>
#include "copa_params.h"
/* Output-changing research path. No parameter or wire-codec constants change. */
#ifdef COPA_SAMPLER_EXP
#if COPA_SAMPLER_EXP!=8 && COPA_SAMPLER_EXP!=16 && COPA_SAMPLER_EXP!=32 && COPA_SAMPLER_EXP!=64 && COPA_SAMPLER_EXP!=128
#error COPA_SAMPLER_EXP must be 8, 16, 32, 64 or 128
#endif
#include "copa_sampler_exp_tables.h"
#else
#include "copa_base_v3_tables.h"
#define COPA_BASE_LIMBS COPA_BASE_LIMBS32
#endif
#endif
