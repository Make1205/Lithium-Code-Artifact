#include "copa_simd.h"
static void copa_expand_d_vectors(int32_t *out,const uint8_t rho[32],unsigned mode){
 unsigned n=copa_label_width(mode),total=copa_rows(mode)*256,per=total/8,s=sigma(mode),inverse[16],range[16],threshold[16],reciprocal[16];
 if(mode==3||mode==5){
  /* Every range is 1/2/4, so every 16-bit candidate is accepted. */
  uint8_t raw[2*1024+4]={0};uint8_t *ptr[8];copa_xof8 engine;
  for(unsigned lane=0;lane<8;++lane)ptr[lane]=raw+2*lane*per;
  copa_xof8_seed_init(&engine,rho,0x14);copa_xof8_prefix(&engine,ptr,2*per);
  for(unsigned j=0;j<n;++j){unsigned row=pi(mode,j),mask=diagonal(mode,row)-1;
   for(unsigned block=0;block<64;block+=COPAV_W){uint32_t index[COPAV_W];for(unsigned k=0;k<COPAV_W;++k)index[k]=2*((block+k)*n+row);
#if COPA_SIMD_BITS==256
    copav_vec v=_mm256_i32gather_epi32((const int *)raw,cv_load(index),1);
#else
    copav_vec v=_mm512_i32gather_epi32(cv_load(index),(const void *)raw,1);
#endif
    cv_store(out+(j/4)*256+(j%4)*64+block,cv_and(v,cv_set(mask)));
   }
  }
  return;
 }

 for(unsigned r=0;r<n;++r)inverse[pi(mode,r)]=r;
 for(unsigned r=0;r<n;++r){range[r]=s*diagonal(mode,r);threshold[r]=65536/range[r]*range[r];reciprocal[r]=(uint32_t)(UINT64_C(0x100000000)/range[r]);}
 copa_xof8 engine;uint8_t blocks[8][136];unsigned counts[8]={0},done=0;
 copa_xof8_seed_init(&engine,rho,0x14);
 while(done<total){copa_xof8_block(&engine,blocks);for(unsigned lane=0;lane<8;++lane){unsigned j=0;
  while(j+2*COPAV_W<=136&&counts[lane]+COPAV_W<=per){uint32_t rr[COPAV_W],tt[COPAV_W],ri[COPAV_W],values[COPAV_W];unsigned first=lane*per+counts[lane];
   for(unsigned k=0;k<COPAV_W;++k){unsigned idx=first+k,rm=idx%n;rr[k]=idx<64*n?range[rm]:3;tt[k]=idx<64*n?threshold[rm]:65535;ri[k]=idx<64*n?reciprocal[rm]:UINT32_C(1431655765);}
#if COPA_SIMD_BITS==256
   copav_vec v=_mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(blocks[lane]+j)));
   unsigned valid=(unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(cv_gt(cv_load(tt),v)))==255;
#else
   copav_vec v=_mm512_cvtepu16_epi32(_mm256_loadu_si256((const __m256i *)(blocks[lane]+j)));
   unsigned valid=_mm512_cmpgt_epi32_mask(cv_load(tt),v)==65535;
#endif
   if(!valid){unsigned idx=first,r=rr[0],value=blocks[lane][j]+256u*blocks[lane][j+1];j+=2;if(value>=tt[0])continue;value%=r;++counts[lane];++done;if(idx<64*n){unsigned row=inverse[idx%n];out[(row/4)*256+idx/n+(row%4)*64]=(int32_t)value;}else out[512+idx-64*n]=(int32_t)value;continue;}
   copav_vec rv=cv_load(rr),q=cv_high32(v,cv_load(ri)),rem=cv_sub(v,cv_mul32(q,rv));rem=cv_sub(rem,cv_and(rv,cv_gt(rem,cv_sub(rv,cv_set(1)))));rem=cv_and(rem,cv_xor(cv_eq(rv,cv_set(1)),cv_set(-1)));cv_store(values,rem);
   for(unsigned k=0;k<COPAV_W;++k){unsigned idx=first+k;if(idx<64*n){unsigned row=inverse[idx%n];out[(row/4)*256+idx/n+(row%4)*64]=(int32_t)values[k];}else out[512+idx-64*n]=(int32_t)values[k];}counts[lane]+=COPAV_W;done+=COPAV_W;j+=2*COPAV_W;
  }
  for(;j<136&&counts[lane]<per;j+=2){unsigned idx=lane*per+counts[lane],r=idx<64*n?range[idx%n]:3,v=blocks[lane][j]+256u*blocks[lane][j+1];if(v>=65536/r*r)continue;v%=r;++counts[lane];++done;if(idx<64*n){unsigned row=inverse[idx%n];out[(row/4)*256+idx/n+(row%4)*64]=(int32_t)v;}else out[512+idx-64*n]=(int32_t)v;}
 }}
}
