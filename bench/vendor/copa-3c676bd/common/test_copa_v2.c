/* Exhaust the installed rounding domain; also exercise hint reconstruction. */
#include "copa_round.h"
#include <assert.h>
#include <stdio.h>
#ifdef COPA_VEC_POLY_HEADER
#include COPA_VEC_POLY_HEADER
#endif
static int high_oracle(int x){return ((x+COPA_ALPHA/2)/COPA_ALPHA)%(COPA_Q/COPA_ALPHA);}
static int wrap_oracle(int x,int q){return (x%q+q)%q;}
static int center_oracle(int x,int q){return x-q*(x>q/2);}
int main(void){
 for(int x=0;x<COPA_Q;++x)assert(copa_highbits(x)==high_oracle(x));
#ifdef COPA_VEC_POLY_HEADER
 for(int first=0;first<COPA_Q;first+=COPAV_W){int32_t x[COPAV_W],h[COPAV_W];for(int j=0;j<COPAV_W;++j)x[j]=(first+j)%COPA_Q;cv_store(h,copa_poly_high(cv_load(x)));for(int j=0;j<COPAV_W;++j)assert(h[j]==high_oracle(x[j]));}
 int32_t w[COPA_M*COPA_N],w1[COPA_M*COPA_N],z2[COPA_M*COPA_N],h[COPA_M*COPA_N],u[COPA_M*COPA_N],reconstructed[COPA_M*COPA_N];uint8_t c[COPA_N]={0};
 for(int first=0;first<COPA_Q;first+=COPA_M*COPA_N){uint64_t want_norm=0;
  for(unsigned j=0;j<COPA_M*COPA_N;++j){w[j]=(first+j)%COPA_Q;w1[j]=high_oracle(w[j]);u[j]=(17*w[j]+COPA_ALPHA/2)%COPA_Q;z2[j]=center_oracle(wrap_oracle(w[j]-u[j],COPA_Q),COPA_Q);int zt=center_oracle(wrap_oracle(COPA_ALPHA*w1[j]-u[j],COPA_Q),COPA_Q);want_norm+=(int64_t)zt*zt;}
  assert(copa_backend_make_hint(h,w,w1,z2)==want_norm);
  for(unsigned j=0;j<COPA_M*COPA_N;++j)assert(h[j]==center_oracle(wrap_oracle(w1[j]-high_oracle(u[j]),COPA_HINT_LEVELS),COPA_HINT_LEVELS));
  assert(copa_backend_reconstruct(reconstructed,u,h,c)==want_norm);
  for(unsigned j=0;j<COPA_M*COPA_N;++j)assert(reconstructed[j]==w1[j]);
 }
#endif
 printf("PASS mode=%d q=%d alpha=%d: all HighBits inputs and SIMD hint roundtrips\n",COPA_MODE,COPA_Q,COPA_ALPHA);
 return 0;
}
