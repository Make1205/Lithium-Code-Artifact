/* Paired one-zeta attempts. A second randombytes means the candidate restarted.
 * setjmp observes the retry without changing either tested implementation. */
#include "copa_api.h"
#include "poly_inv.h"
#include "fips202.h"
#include <assert.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static jmp_buf retry;
static uint8_t zeta[32],pk[2][COPA_PK_BYTES],wire[2][COPA_SK_BYTES];
static copa_secret_key key[2];static copa_keygen_stats stats[2];
static unsigned draws,accepted;
int copa_keygen_paper_raw(uint8_t *,copa_secret_key *,copa_keygen_stats *);
void __wrap_randombytes(uint8_t *out,size_t n){assert(n==32);if(draws++)longjmp(retry,1);memcpy(out,zeta,32);}
static unsigned capture(unsigned which){
 draws=accepted=0;
 if(!setjmp(retry)){int r=which?copa_keygen_raw(pk[1],&key[1],&stats[1]):copa_keygen_paper_raw(pk[0],&key[0],&stats[0]);assert(r==0);accepted=1;}
 return accepted;
}
int main(int argc,char **argv){
 uint64_t count=argc>1?strtoull(argv[1],0,10):1000000,start=argc>2?strtoull(argv[2],0,10):0;
 int32_t f[COPA_N]={0},inv[COPA_N];f[0]=2;assert(!copa_poly_inv(inv,f,COPA_Q));
 f[0]=-1;f[64]=COPA_MODE==3?3:1;f[128]=1;assert(!copa_poly_inv(inv,f,COPA_Q));
 uint32_t parity=0;for(unsigned j=0;j<COPA_N;++j)parity^=(uint32_t)f[j];assert(parity&1u);
 uint64_t yes=0,even=0,odd=0,norm=0,other=0,early=0;keccak_state transcript;shake256_init(&transcript);
 for(uint64_t i=start;i<start+count;++i){uint8_t seed[24]={0};memcpy(seed,"KeyGen-order-v31",16);for(unsigned j=0;j<8;++j)seed[16+j]=(uint8_t)(i>>(8*j));shake256(zeta,32,seed,sizeof(seed));
  unsigned a=capture(0),b=capture(1);assert(a==b);early+=stats[1].partial_norm_reject;
  uint32_t sum=0;for(unsigned j=0;j<COPA_N;++j)sum^=(uint32_t)key[0].s[0][j];
  if(stats[0].noninvertible){if(sum&1u)++odd;else ++even;}
  norm+=stats[0].norm_reject;other+=stats[0].lower_norm_reject+stats[0].packing_failure;
  shake256_absorb(&transcript,zeta,32);uint8_t decision=(uint8_t)a;shake256_absorb(&transcript,&decision,1);
  if(a){++yes;assert(!copa_pack_sk(wire[0],&key[0]));assert(!copa_pack_sk(wire[1],&key[1]));assert(!memcmp(pk[0],pk[1],COPA_PK_BYTES));assert(!memcmp(wire[0],wire[1],COPA_SK_BYTES));shake256_absorb(&transcript,wire[0],COPA_SK_BYTES);}
  if((i-start+1)%10000==0){fprintf(stderr,"mode=%d checked=%llu accepted=%llu odd_nonunit=%llu\n",COPA_MODE,(unsigned long long)(i-start+1),(unsigned long long)yes,(unsigned long long)odd);fflush(stderr);}
 }
 uint8_t hash[32];shake256_finalize(&transcript);shake256_squeeze(hash,32,&transcript);
 printf("{\"new_partial_rejects\":%llu,",(unsigned long long)early);
 printf("\"mode\":%d,\"attempts\":%llu,\"start\":%llu,\"accepted\":%llu,\"nonunit_even\":%llu,\"nonunit_odd\":%llu,\"norm_reject\":%llu,\"other\":%llu,\"known_even_and_odd_nonunits\":true,\"decisions_and_accepted_keys_identical\":true,\"transcript_shake256_256\":\"",COPA_MODE,(unsigned long long)count,(unsigned long long)start,(unsigned long long)yes,(unsigned long long)even,(unsigned long long)odd,(unsigned long long)norm,(unsigned long long)other);
 for(unsigned j=0;j<32;++j){printf("%02x",hash[j]);}puts("\"}");return 0;
}
