#ifndef COPA_CDT_VECTORS_H
#define COPA_CDT_VECTORS_H
#include "copa_simd.h"
#define CG_LANES (COPAV_W/2)
static inline copav_vec cg_ult(copav_vec a,copav_vec b){copav_vec sign=cv_set64(UINT64_C(1)<<63);return cv_gt64(cv_xor(b,sign),cv_xor(a,sign));}
static inline copav_vec cg_eq(copav_vec a,copav_vec b){
#if COPA_SIMD_BITS==256
 return _mm256_cmpeq_epi64(a,b);
#else
 return _mm512_maskz_set1_epi64(_mm512_cmpeq_epi64_mask(a,b),-1);
#endif
}
#if defined(COPA_CDT_SINGLE_GROUP) || (COPA_SIMD_BITS==256 && COPA_MODE!=2 && !defined(COPA_CDT_FORCE_PAIR))
static inline void copa_cdt_vectors(uint32_t out[32],const uint64_t u[3][32],const uint64_t (*cdf)[3],unsigned count){
 for(unsigned g=0;g<32;g+=CG_LANES){copav_vec u0=cv_load(u[0]+g),u1=cv_load(u[1]+g),u2=cv_load(u[2]+g),sum=cv_set(0);
  for(unsigned k=0;k<count;++k){copav_vec t0=cv_set64(cdf[k][0]),t1=cv_set64(cdf[k][1]),t2=cv_set64(cdf[k][2]);copav_vec lt=cv_or(cg_ult(u2,t2),cv_and(cg_eq(u2,t2),cv_or(cg_ult(u1,t1),cv_and(cg_eq(u1,t1),cg_ult(u0,t0)))));sum=cv_add64(sum,cv_and(lt,cv_set64(1)));}
  uint64_t words[CG_LANES];cv_store(words,sum);for(unsigned j=0;j<CG_LANES;++j)out[g+j]=(uint32_t)words[j];
 }
}
#else
static inline void copa_cdt_vectors(uint32_t out[32],const uint64_t u[3][32],const uint64_t (*cdf)[3],unsigned count){
 const uint64_t sign=UINT64_C(1)<<63;
 for(unsigned g=0;g<32;g+=2*CG_LANES){copav_vec u0[2],u1[2],u2[2],sum[2];
  for(unsigned j=0;j<2;++j){u0[j]=cv_xor(cv_load(u[0]+g+j*CG_LANES),cv_set64(sign));u1[j]=cv_xor(cv_load(u[1]+g+j*CG_LANES),cv_set64(sign));u2[j]=cv_xor(cv_load(u[2]+g+j*CG_LANES),cv_set64(sign));sum[j]=cv_set(0);}
  for(unsigned k=0;k<count;++k){copav_vec t0=cv_set64(cdf[k][0]^sign),t1=cv_set64(cdf[k][1]^sign),t2=cv_set64(cdf[k][2]^sign);
   for(unsigned j=0;j<2;++j){copav_vec lt=cv_or(cv_gt64(t2,u2[j]),cv_and(cg_eq(u2[j],t2),cv_or(cv_gt64(t1,u1[j]),cv_and(cg_eq(u1[j],t1),cv_gt64(t0,u0[j])))));sum[j]=cv_add64(sum[j],cv_and(lt,cv_set64(1)));}
  }
  for(unsigned j=0;j<2;++j){uint64_t words[CG_LANES];cv_store(words,sum[j]);for(unsigned k=0;k<CG_LANES;++k)out[g+j*CG_LANES+k]=(uint32_t)words[k];}
 }
}
#endif
#endif
