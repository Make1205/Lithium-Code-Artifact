#ifndef COPA_PROF_H
#define COPA_PROF_H
#include <stdint.h>
enum{CP_PARSE_SK,CP_Y,CP_AY,CP_HIGH_CH,CP_IRS,CP_BS,CP_HINT_BV,CP_RANS,CP_PK_UNLAB,CP_A,CP_D,CP_MATRIX_NTT,CP_DECODE,CP_VF_NTT,CP_CHALLENGE,CP_RECONSTRUCT,CP_S,CP_INVERSE,CP_QUANT,CP_PACK,CP_NTT_FORWARD,CP_NTT_POINTWISE,CP_NTT_INVERSE,CP_SK_DECODE,CP_SK_HASH,CP_SK_NTT,CP_IRS_PREP,CP_SIGN_SETUP,CP_B_REBUILD,CP_Y_SHAKE,CP_IRS_PROJECT,CP_IRS_TRANSITION,CP_IRS_ROTATE,CP_UNLAB,CP_B_NTT,CP_A_NTT,CP_VF_NORM,CP_KG_HASH,CP_KG_NORM,CP_KG_B,CP_KG_MATRIX,CP_COUNT};
#ifdef COPA_PROFILE
#include <x86intrin.h>
extern uint64_t copa_profile_cycles[CP_COUNT],copa_profile_calls[CP_COUNT];
static inline uint64_t copa_profile_tick(void){unsigned aux;_mm_lfence();return __rdtscp(&aux);}
#define CP_BEGIN(v) uint64_t v=copa_profile_tick()
#define CP_END(id,v) do{copa_profile_cycles[id]+=copa_profile_tick()-(v);++copa_profile_calls[id];}while(0)
#else
#define CP_BEGIN(v) ((void)0)
#define CP_END(id,v) ((void)0)
#endif
#endif
