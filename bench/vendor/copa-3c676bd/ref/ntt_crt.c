#include "ntt_crt.h"
#include "copa_scalar64.h"
#include <assert.h>
#include <string.h>
#define P1 UINT32_C(8380417)
#define P2 UINT32_C(7340033)
#ifndef COPA_EXTERNAL_NTT
typedef struct {uint32_t p,qinv;const uint32_t *twist,*post,*forward,*inverse,*twist_shoup,*post_shoup,*forward_shoup,*inverse_shoup,*zeta,*izeta,*zeta_shoup,*izeta_shoup,*scale,*scale_shoup;} ntt_prime;
#include "copa_ntt_cache_tables.h"
static inline uint32_t mont(uint32_t a,uint32_t b,const ntt_prime *p){
 uint64_t z=(uint64_t)a*b;uint32_t t=(uint32_t)z*p->qinv;
 uint32_t u=(uint32_t)((z+(uint64_t)t*p->p)>>32);return u-(p->p & (uint32_t)copa64_mask(u>=p->p));
}
static inline uint32_t shoup(uint32_t a,uint32_t w,uint32_t ws,uint32_t p){
 uint32_t q=(uint32_t)(((uint64_t)a*ws)>>32),r=a*w-q*p;return r-p*(uint32_t)(r>=p);
}
#endif
static uint64_t reciprocal(uint32_t p){
 /* These public moduli cover the scheme and norm-tower arithmetic. */
 switch(p){
 case 7680:return UINT64_C(2401919801264264);
 case 12584:return UINT64_C(1465888753473422);
 case 6048:return UINT64_C(3050056890494304);
 case 8400:return UINT64_C(2196040961155899);
 default:{uint64_t q=UINT64_MAX/p;return q+(UINT64_MAX%p==p-1);}
 }
}
uint32_t copa_mod_i64(int64_t x,uint32_t p){
 uint64_t mask=copa64_mask(x<0),a=((uint64_t)x^mask)-mask,lo;
 uint64_t q=copa64_product(a,reciprocal(p),&lo),r=a-q*p;r-=(uint64_t)p & copa64_mask(r>=p);
 uint32_t neg=(p-(uint32_t)r)&(uint32_t)copa64_mask(r!=0);return (neg&(uint32_t)mask)|((uint32_t)r&~(uint32_t)mask);
}
#ifndef COPA_EXTERNAL_NTT
void copa_ntt_forward(uint32_t out[2][COPA_N],const int32_t *in,size_t n){
 assert(n&&n<=COPA_N&&!(n&(n-1)));
 for(unsigned k=0;k<2;++k){const ntt_prime *p=&ntt_primes[k];
  for(size_t j=0;j<n;++j)out[k][j]=(uint32_t)in[j]+(p->p&(uint32_t)copa64_mask(in[j]<0));
  for(size_t len=n/2;len>=8;len>>=1){size_t root=n/(2*len);
   for(size_t b=0;b<n;b+=2*len,++root){uint32_t w=p->zeta[n-1+root],ws=p->zeta_shoup[n-1+root];
    for(size_t j=b;j<b+len;++j){uint32_t a=out[k][j],v=shoup(out[k][j+len],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][j]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][j+len]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
   }
  }
  if(n>=8)for(size_t b=0,root=n/8;b<n;b+=8,++root){uint32_t w=p->zeta[n-1+root],ws=p->zeta_shoup[n-1+root];
   {uint32_t a=out[k][b+0],v=shoup(out[k][b+4],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+0]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+4]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
   {uint32_t a=out[k][b+1],v=shoup(out[k][b+5],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+1]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+5]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
   {uint32_t a=out[k][b+2],v=shoup(out[k][b+6],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+2]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+6]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
   {uint32_t a=out[k][b+3],v=shoup(out[k][b+7],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+3]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+7]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
  }
  if(n>=4)for(size_t b=0,root=n/4;b<n;b+=4,++root){uint32_t w=p->zeta[n-1+root],ws=p->zeta_shoup[n-1+root];
   {uint32_t a=out[k][b+0],v=shoup(out[k][b+2],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+0]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+2]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
   {uint32_t a=out[k][b+1],v=shoup(out[k][b+3],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+1]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+3]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
  }
  if(n>=2)for(size_t b=0,root=n/2;b<n;b+=2,++root){uint32_t w=p->zeta[n-1+root],ws=p->zeta_shoup[n-1+root];
   {uint32_t a=out[k][b+0],v=shoup(out[k][b+1],w,ws,p->p),sum=a+v,diff=a+p->p-v;out[k][b+0]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));out[k][b+1]=diff-(p->p & (uint32_t)copa64_mask(diff>=p->p));}
  }
 }
}
void copa_ntt_dot(uint32_t out[2][COPA_N],const uint32_t a[][2][COPA_N],const uint32_t b[][2][COPA_N],unsigned count,size_t n){
 for(unsigned k=0;k<2;++k){const ntt_prime *p=&ntt_primes[k];memset(out[k],0,n*sizeof(uint32_t));
  for(unsigned r=0;r<count;++r)for(size_t j=0;j<n;++j){uint32_t sum=out[k][j]+mont(a[r][k][j],b[r][k][j],p);out[k][j]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));}
 }
}
void copa_ntt_dot_reversed(uint32_t out[2][COPA_N],const uint32_t a[][2][COPA_N],const uint32_t b[][2][COPA_N],unsigned count,size_t n){
 for(unsigned k=0;k<2;++k){const ntt_prime *p=&ntt_primes[k];memset(out[k],0,n*sizeof(uint32_t));
  for(unsigned r=0;r<count;++r)for(size_t j=0;j<n;++j){uint32_t sum=out[k][j]+mont(a[r][k][j],b[r][k][n-1-j],p);out[k][j]=sum-(p->p & (uint32_t)copa64_mask(sum>=p->p));}
 }
}
static void inverse_residues(uint32_t v[2][COPA_N],const uint32_t in[2][COPA_N],size_t n){
 assert(n&&n<=COPA_N&&!(n&(n-1)));unsigned logn=(unsigned)__builtin_ctz((unsigned)n);
 for(unsigned k=0;k<2;++k){const ntt_prime *p=&ntt_primes[k];memcpy(v[k],in[k],n*sizeof(uint32_t));
  if(n>=2)for(size_t b=0,root=n/2;b<n;b+=2,++root){uint32_t w=p->izeta[n-1+root],ws=p->izeta_shoup[n-1+root];
   {uint32_t a=v[k][b+0],t=v[k][b+1],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+0]=sum;v[k][b+1]=shoup(diff,w,ws,p->p);}
  }
  if(n>=4)for(size_t b=0,root=n/4;b<n;b+=4,++root){uint32_t w=p->izeta[n-1+root],ws=p->izeta_shoup[n-1+root];
   {uint32_t a=v[k][b+0],t=v[k][b+2],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+0]=sum;v[k][b+2]=shoup(diff,w,ws,p->p);}
   {uint32_t a=v[k][b+1],t=v[k][b+3],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+1]=sum;v[k][b+3]=shoup(diff,w,ws,p->p);}
  }
  if(n>=8)for(size_t b=0,root=n/8;b<n;b+=8,++root){uint32_t w=p->izeta[n-1+root],ws=p->izeta_shoup[n-1+root];
   {uint32_t a=v[k][b+0],t=v[k][b+4],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+0]=sum;v[k][b+4]=shoup(diff,w,ws,p->p);}
   {uint32_t a=v[k][b+1],t=v[k][b+5],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+1]=sum;v[k][b+5]=shoup(diff,w,ws,p->p);}
   {uint32_t a=v[k][b+2],t=v[k][b+6],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+2]=sum;v[k][b+6]=shoup(diff,w,ws,p->p);}
   {uint32_t a=v[k][b+3],t=v[k][b+7],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][b+3]=sum;v[k][b+7]=shoup(diff,w,ws,p->p);}
  }
  for(size_t len=8;len<n;len<<=1){size_t root=n/(2*len);
   for(size_t b=0;b<n;b+=2*len,++root){uint32_t w=p->izeta[n-1+root],ws=p->izeta_shoup[n-1+root];
    for(size_t j=b;j<b+len;++j){uint32_t a=v[k][j],t=v[k][j+len],sum=a+t,diff=a+p->p-t;sum-=(p->p & (uint32_t)copa64_mask(sum>=p->p));diff-=(p->p & (uint32_t)copa64_mask(diff>=p->p));v[k][j]=sum;v[k][j+len]=shoup(diff,w,ws,p->p);}
   }
  }
  for(size_t j=0;j<n;++j)v[k][j]=shoup(v[k][j],p->scale[logn],p->scale_shoup[logn],p->p);
 }
}
void copa_ntt_inverse_crt(int64_t *out,const uint32_t in[2][COPA_N],size_t n){
 uint32_t v[2][COPA_N];inverse_residues(v,in,n);
 const uint64_t product=(uint64_t)P1*P2;
 for(size_t j=0;j<n;++j){int64_t d=(int64_t)v[1][j]-v[0][j];
  uint32_t u=(uint32_t)(d+(int64_t)((2*P2)&copa64_mask(d<0)));u-=(P2 & (uint32_t)copa64_mask(u>=P2));
  uint32_t t=mont(u,NTT_P1_INV_MONT_P2,&ntt_primes[1]);uint64_t x=v[0][j]+(uint64_t)P1*t;
  out[j]=(int64_t)x-(int64_t)(product & copa64_ct_mask(x>product/2));
 }
}
void copa_ntt_crt_mul(int64_t *out,const int32_t *a,const int32_t *b,size_t n){
 uint32_t aa[1][2][COPA_N],bb[1][2][COPA_N],cc[2][COPA_N];
 copa_ntt_forward(aa[0],a,n);copa_ntt_forward(bb[0],b,n);copa_ntt_dot(cc,aa,bb,1,n);copa_ntt_inverse_crt(out,cc,n);
}

/* Fused centered CRT and public-q reduction; all products fit in uint64_t. */
void copa_ntt_inverse_modq(int32_t *out,const uint32_t in[2][COPA_N],size_t n,uint32_t q){
 uint32_t v[2][COPA_N];inverse_residues(v,in,n);
 const uint64_t product=(uint64_t)P1*P2;
 const uint32_t recip=(uint32_t)(UINT64_C(0x100000000)/q),p1q=P1%q,pq=(uint32_t)(product%q);
 for(size_t j=0;j<n;++j){
  int64_t d=(int64_t)v[1][j]-v[0][j];
  uint32_t u=(uint32_t)(d+(int64_t)((2*P2)&copa64_mask(d<0)));u-=(P2 & (uint32_t)copa64_mask(u>=P2));
  uint32_t t=mont(u,NTT_P1_INV_MONT_P2,&ntt_primes[1]);
  uint32_t negative=(uint32_t)(copa64_mask(t>P2/2)|(copa64_mask(t==P2/2)&copa64_mask(v[0][j]>P1/2)));
  uint32_t a=v[0][j]-(uint32_t)(((uint64_t)v[0][j]*recip)>>32)*q;a-=(q & (uint32_t)copa64_mask(a>=q));
  uint32_t b=t-(uint32_t)(((uint64_t)t*recip)>>32)*q;b-=(q & (uint32_t)copa64_mask(b>=q));
  uint32_t c=a+p1q*b;c-=(uint32_t)(((uint64_t)c*recip)>>32)*q;c-=(q & (uint32_t)copa64_mask(c>=q));
  int32_t z=(int32_t)c-(int32_t)(pq&negative);out[j]=z+(int32_t)(q&(uint32_t)copa64_mask(z<0));
 }
}

#endif /* COPA_EXTERNAL_NTT: shared reduction, replaceable transform engine */
