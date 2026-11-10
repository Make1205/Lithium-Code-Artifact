#ifndef COPA_SAMPLER_PROFILE_H
#define COPA_SAMPLER_PROFILE_H
#include <stdint.h>
enum {SP_Y_SHAKE,SP_BASE,SP_EXP,SP_SELECT,SP_IRS_DOT,SP_IRS_PROB,SP_IRS_UPDATE,SP_IRS_MUL,SP_IRS_SHAKE,SP_DECOMPOSE,SP_CHASH,SP_SETUP,SP_COUNT};
enum {SC_Y_CALLS,SC_CANDIDATES,SC_ACCEPTED,SC_EXP_ACCEPTED,SC_Y_USED,SC_Y_SQUEEZED,SC_WHOLE_RETRIES,SC_COUNT};
#ifdef COPA_SAMPLER_PROFILE
#include <x86intrin.h>
extern uint64_t copa_sp_cycles[SP_COUNT],copa_sp_counts[SC_COUNT];
static inline uint64_t copa_sp_tick(void){unsigned a;_mm_lfence();uint64_t t=__rdtscp(&a);_mm_lfence();return t;}
#define SP_BEGIN(v) uint64_t v=copa_sp_tick()
#define SP_END(id,v) (copa_sp_cycles[id]+=copa_sp_tick()-(v))
#define SC_ADD(id,n) (copa_sp_counts[id]+=(n))
#else
#define SP_BEGIN(v) ((void)0)
#define SP_END(id,v) ((void)0)
#define SC_ADD(id,n) ((void)0)
#endif
#endif
