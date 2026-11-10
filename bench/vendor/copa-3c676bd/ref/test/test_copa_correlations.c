#include "copa_sampler.h"
#include "copa_keygen.h"
#include "ntt_crt.h"
#include "fips202.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
static uint64_t random_state=1;
static uint64_t rnd(void){random_state^=random_state<<13;random_state^=random_state>>7;random_state^=random_state<<17;return random_state;}
static uint64_t load(const uint8_t *p){uint64_t v=0;for(unsigned j=0;j<8;++j)v|=(uint64_t)p[j]<<(8*j);return v;}
static int full_vector(int32_t z[COPA_K][COPA_N],const uint8_t ch[COPA_N],const copa_secret_key *sk,const uint8_t seed[32],uint64_t attempt,unsigned *steps){
 uint8_t ext[41];memcpy(ext,seed,32);ext[32]=0x15;for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));keccak_state rng;shake256_absorb_once(&rng,ext,sizeof(ext));*steps=0;
 for(unsigned j=0;j<COPA_N;++j)if(ch[j]){
  int64_t dot=0;int32_t v[COPA_K][COPA_N];for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){v[r][i]=i>=j?sk->s[r][i-j]:-sk->s[r][i+COPA_N-j];dot+=(int64_t)z[r][i]*v[r][i];}
  uint8_t rb[48];uint64_t words[6];shake256_squeeze(rb,sizeof(rb),&rng);for(unsigned b=0;b<6;++b)words[b]=load(rb+8*b);
  unsigned minus,ok=copa_irs64_step(&minus,dot,&sk->irs_context,words);++*steps;if(!ok)return 0;
  int32_t sign=1-2*(int32_t)minus;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)z[r][i]+=sign*v[r][i];
 }
 return 1;
}
int main(void){
 unsigned accepted=0,rejected=0;
 for(unsigned t=0;t<100;++t){copa_secret_key sk;memset(&sk,0,sizeof(sk));
  uint64_t norm=0;for(unsigned j=0;j<1+t%64;++j){unsigned i=(unsigned)(rnd()%(COPA_K*COPA_N));if((&sk.s[0][0])[i]==0){(&sk.s[0][0])[i]=(rnd()&1)?1:-1;++norm;}}
  assert(norm);for(unsigned r=0;r<COPA_K;++r)copa_ntt_forward(sk.secret_ntt[r],sk.s[r],COPA_N);copa_irs64_prepare(&sk.irs_context,norm);
  int32_t z[COPA_K][COPA_N],copy[COPA_K][COPA_N];uint32_t zntt[COPA_K][2][COPA_N],ct[2][COPA_N];int64_t dots[COPA_N];
  for(unsigned r=0;r<COPA_K;++r){for(unsigned i=0;i<COPA_N;++i)z[r][i]=(int32_t)(rnd()%201)-100;copa_ntt_forward(zntt[r],z[r],COPA_N);}
  copa_ntt_dot_reversed(ct,zntt,sk.secret_ntt,COPA_K,COPA_N);copa_ntt_inverse_crt(dots,ct,COPA_N);
  copa_ntt_dot_reversed(ct,sk.secret_ntt,sk.secret_ntt,COPA_K,COPA_N);copa_ntt_inverse_crt(sk.autocorrelation,ct,COPA_N);
  for(unsigned j=0;j<COPA_N;++j){int64_t d=0,a=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){int32_t v=i>=j?sk.s[r][i-j]:-sk.s[r][i+COPA_N-j];d+=(int64_t)z[r][i]*v;a+=(int64_t)sk.s[r][i]*v;}assert(d==dots[j]&&a==sk.autocorrelation[j]);}
  uint8_t ch[COPA_N]={0},seed[32];for(unsigned i=0;i<32;++i)seed[i]=(uint8_t)rnd();for(unsigned i=COPA_N-COPA_TAU;i<COPA_N;++i){unsigned j=(unsigned)(rnd()%(i+1));ch[i]=ch[j];ch[j]=1;}
  memcpy(copy,z,sizeof(z));unsigned a,b;int full=full_vector(copy,ch,&sk,seed,t,&a),fast=copa_irs(z,ch,&sk,zntt,seed,t,&b);assert(full==fast&&a==b);if(full){assert(!memcmp(z,copy,sizeof(z)));++accepted;}else ++rejected;
 }assert(accepted&&rejected);printf("PASS correlations mode=%d: all initial/autocorrelation dots and 100 full-vector/early-IRS equivalence cases; accept=%u reject=%u\n",COPA_MODE,accepted,rejected);return 0;
}
