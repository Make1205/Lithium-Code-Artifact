#ifndef COPA_IRS_SETUP_VECTORS_H
#define COPA_IRS_SETUP_VECTORS_H
#include "copa_simd.h"
#include "copa_irs64_soa.h"
static inline void copa_backend_irs_constants(uint64_t selected[5],uint64_t norm){
 const unsigned lanes=COPAV_W/2;uint64_t indices[COPAV_W/2];for(unsigned j=0;j<lanes;++j)indices[j]=j+1;
 copav_vec sums[5];for(unsigned c=0;c<5;++c)sums[c]=cv_set(0);copav_vec index=cv_load(indices),target=cv_set64(norm);
 for(unsigned j=0;j<COPA_BK_SQ;j+=lanes){
#if COPA_SIMD_BITS==256
  copav_vec mask=_mm256_cmpeq_epi64(index,target);
#else
  copav_vec mask=_mm512_maskz_set1_epi64(_mm512_cmpeq_epi64_mask(index,target),-1);
#endif
  for(unsigned c=0;c<5;++c)sums[c]=cv_or(sums[c],cv_and(mask,cv_load(copa_irs_soa[c]+j)));
  index=cv_add64(index,cv_set64(lanes));
 }
 for(unsigned c=0;c<5;++c){uint64_t words[COPAV_W/2];cv_store(words,sums[c]);selected[c]=0;for(unsigned j=0;j<lanes;++j)selected[c]|=words[j];}
}
#endif
