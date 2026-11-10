#include <string.h>
/* Optional Verify-only auxiliary prime. No parameter/wire changes.
 * For each valid z1, ||z1||<=Bv and |A_ij|<q imply
 * |(A1*z1)_i| <= q*sqrt((1+ell)*n)*Bv < P/2.
 * Entry points REQUIRE the norm gate before multiplication. */
#ifdef COPA_VEC_LATTICE
#include "copa_simd.h"
#endif
#include "copa_verify31_tables.h"
_Static_assert(UINT64_C(1)*(1+COPA_ELL)*COPA_N*COPA_Q*COPA_Q*COPA_BV*COPA_BV < UINT64_C(1006632960)*1006632960,"Verify-only single-prime Cauchy bound");
#ifdef COPA_VEC_LATTICE
#include "copa_simd.h"
#include "copa_ntt_zip.h"
static inline copav_vec q31_vreduce(copav_vec x){
#if COPAV_W==8
 return _mm256_min_epu32(x,cv_sub(x,cv_set(Q31_P)));
#else
 return _mm512_min_epu32(x,cv_sub(x,cv_set(Q31_P)));
#endif
}
static inline copav_vec q31_vmul(copav_vec x,copav_vec w,copav_vec ws){return cv_sub(cv_mul32(x,w),cv_mul32(cv_high32(x,ws),cv_set(Q31_P)));}
static inline void q31_butterfly(copav_vec *a,copav_vec *b,copav_vec w,copav_vec ws,unsigned inverse){
 copav_vec x=q31_vreduce(*a),y;
 if(!inverse){y=q31_vreduce(q31_vmul(*b,w,ws));*a=cv_add(x,y);*b=cv_sub(cv_add(x,cv_set(Q31_P)),y);}
 else{y=q31_vreduce(*b);*a=cv_add(x,y);*b=q31_vmul(cv_sub(cv_add(x,cv_set(Q31_P)),y),w,ws);}
}
#else
static inline uint32_t q31_reduce(uint32_t x){return x-Q31_P*(x>=Q31_P);}
static inline uint32_t q31_mul(uint32_t x,uint32_t w,uint32_t ws){return (uint32_t)(x*w-(uint32_t)(((uint64_t)x*ws)>>32)*Q31_P);}
#endif
static void q31_stages(uint32_t x[COPA_N],unsigned inverse){
 const uint32_t *roots=inverse?q31_iroot:q31_root,*shoup=inverse?q31_iroot_shoup:q31_root_shoup;
 for(unsigned stage=0;stage<8;++stage){unsigned h=inverse?(1u<<stage):(128u>>stage);
#ifdef COPA_VEC_LATTICE
  if(h<COPAV_W){unsigned index=(unsigned)__builtin_ctz(h)*128;const uint32_t *w=inverse?q31_inv_packed:q31_fwd_packed,*ws=inverse?q31_inv_packed_shoup:q31_fwd_packed_shoup;
   for(unsigned b=0;b<256;b+=2*COPAV_W){copav_vec lo=cv_load(x+b),hi=cv_load(x+b+COPAV_W),a=ntt_zip_low(lo,hi,h),c=ntt_zip_high(lo,hi,h);q31_butterfly(&a,&c,cv_load(w+index),cv_load(ws+index),inverse);cv_store(x+b,ntt_unzip_low(a,c,h));cv_store(x+b+COPAV_W,ntt_unzip_high(a,c,h));index+=COPAV_W;}
  }else for(unsigned b=0;b<256;b+=2*h){unsigned root=256/(2*h)+b/(2*h);copav_vec w=cv_set(roots[root]),ws=cv_set(shoup[root]);for(unsigned j=0;j<h;j+=COPAV_W){copav_vec a=cv_load(x+b+j),c=cv_load(x+b+j+h);q31_butterfly(&a,&c,w,ws,inverse);cv_store(x+b+j,a);cv_store(x+b+j+h,c);}}
#else
  for(unsigned b=0;b<256;b+=2*h){unsigned root=256/(2*h)+b/(2*h);for(unsigned j=0;j<h;++j){uint32_t a=q31_reduce(x[b+j]),c;
   if(!inverse){c=q31_reduce(q31_mul(x[b+j+h],roots[root],shoup[root]));x[b+j]=a+c;x[b+j+h]=a+Q31_P-c;}
   else{c=q31_reduce(x[b+j+h]);x[b+j]=a+c;x[b+j+h]=q31_mul(a+Q31_P-c,roots[root],shoup[root]);}
  }}
#endif
 }
}
static void q31_forward_ntt(uint32_t out[COPA_N],const int32_t in[COPA_N]){
#ifdef COPA_VEC_LATTICE
 for(unsigned j=0;j<256;j+=COPAV_W){copav_vec x=cv_load(in+j);x=cv_add(x,cv_and(cv_sra32(x,31),cv_set(Q31_P)));cv_store(out+j,x);}
#else
 for(unsigned j=0;j<256;++j){uint32_t x=(uint32_t)in[j]+(Q31_P & (0u-(uint32_t)(in[j]<0)));out[j]=x;}
#endif
 q31_stages(out,0);
}
static void q31_cache_impl(uint32_t out[2][COPA_N],const int32_t in[COPA_N],unsigned positive){
 if(positive){memcpy(out[0],in,COPA_N*sizeof(uint32_t));q31_stages(out[0],0);}else q31_forward_ntt(out[0],in);
#ifdef COPA_VEC_LATTICE
 const uint64_t reciprocal=UINT64_C(9162596893);
 for(unsigned j=0;j<256;j+=COPAV_W){copav_vec x=q31_vreduce(cv_load(out[0]+j)),q=cv_add(cv_high32(x,cv_set((uint32_t)reciprocal)),cv_sll32(x,1));
  copav_vec rem=cv_sub(cv_set(0),cv_mul32(q,cv_set(Q31_P))),r=q31_vreduce(rem),correct=cv_and(cv_xor(cv_eq(rem,r),cv_set(-1)),cv_set(1));cv_store(out[0]+j,x);cv_store(out[1]+j,cv_add(q,correct));}
#else
 for(unsigned j=0;j<256;++j){out[0][j]=q31_reduce(out[0][j]);out[1][j]=(uint32_t)(((uint64_t)out[0][j]<<32)/Q31_P);}
#endif
}
static void q31_cache(uint32_t out[2][COPA_N],const int32_t in[COPA_N]){q31_cache_impl(out,in,0);}
#ifndef COPA_VERIFY31_CORE_ONLY
int copa_verify31_expand(copa_public_ntt matrix,const uint8_t pk[COPA_PK_BYTES]){
 uint8_t rho[32];int32_t bp[COPA_M][COPA_N],d[COPA_M][COPA_N],a[COPA_M][COPA_ELL][COPA_N];
 CP_BEGIN(pu);int ok=copa_unpack_pk(rho,&bp[0][0],pk,COPA_PK_BYTES,COPA_MODE);CP_END(CP_PK_UNLAB,pu);if(!ok)return 0;
 CP_BEGIN(pa);copa_expand_a(a,rho);CP_END(CP_A,pa);CP_BEGIN(pd);copa_expand_d(&d[0][0],rho,COPA_MODE);CP_END(CP_D,pd);CP_BEGIN(pntt);
 for(unsigned r=0;r<COPA_M;++r){int32_t first[COPA_N];for(unsigned i=0;i<256;++i){int32_t t=bp[r][i]+d[r][i];t-=COPA_Q*(t>=COPA_Q);t=-t;first[i]=t;}
  CP_BEGIN(pb);q31_cache(matrix[r][0],first);CP_END(CP_B_NTT,pb);CP_BEGIN(pan);
  for(unsigned j=0;j<COPA_ELL;++j){q31_cache_impl(matrix[r][1+j],a[r][j],1);}CP_END(CP_A_NTT,pan);
 }CP_END(CP_MATRIX_NTT,pntt);return 1;
}
void copa_verify31_apply(int32_t out[COPA_M][COPA_N],const copa_public_ntt matrix,const int32_t z[1+COPA_ELL][COPA_N]){
 uint32_t transformed[1+COPA_ELL][COPA_N];for(unsigned j=0;j<1+COPA_ELL;++j)q31_forward_ntt(transformed[j],z[j]);
 for(unsigned r=0;r<COPA_M;++r){uint32_t value[COPA_N];
#ifdef COPA_VEC_LATTICE
  for(unsigned j=0;j<256;j+=COPAV_W){copav_vec sum=cv_set(0);for(unsigned k=0;k<1+COPA_ELL;++k)sum=q31_vreduce(cv_add(sum,q31_vreduce(q31_vmul(cv_load(transformed[k]+j),cv_load(matrix[r][k][0]+j),cv_load(matrix[r][k][1]+j)))));cv_store(value+j,sum);}
#else
  for(unsigned j=0;j<256;++j){uint32_t sum=0;for(unsigned k=0;k<1+COPA_ELL;++k)sum=q31_reduce(sum+q31_reduce(q31_mul(transformed[k][j],matrix[r][k][0][j],matrix[r][k][1][j])));value[j]=sum;}
#endif
  q31_stages(value,1);
#ifdef COPA_VEC_LATTICE
  for(unsigned j=0;j<256;j+=COPAV_W){copav_vec x=q31_vreduce(q31_vmul(cv_load(value+j),cv_load(q31_post+j),cv_load(q31_post_shoup+j))),negative=cv_gt(x,cv_set(Q31_P/2));copav_vec mag=cv_select(cv_sub(cv_set(Q31_P),x),x,negative),rem=cv_sub(mag,cv_mul32(qv_div_positive(mag,COPA_Q),cv_set(COPA_Q)));copav_vec neg=cv_and(cv_sub(cv_set(COPA_Q),rem),cv_xor(cv_eq(rem,cv_set(0)),cv_set(-1)));cv_store(out[r]+j,cv_select(neg,rem,negative));}
#else
  for(unsigned j=0;j<256;++j){uint32_t v=q31_reduce(q31_mul(value[j],q31_post[j],q31_post_shoup[j]));int64_t x=(int64_t)v-(int64_t)Q31_P*(v>Q31_P/2);out[r][j]=(int32_t)copa_mod_i64(x,COPA_Q);}
#endif
 }
}

#endif /* production entry points */
