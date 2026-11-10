/* Reproducible acceptance audit: 1000 independent keys, 100 messages each. */
#include "copa_api.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "copa_kat/nist_rng.h"
int main(int argc,char **argv){
 unsigned trials=argc>1?(unsigned)strtoul(argv[1],0,10):100000;
 uint8_t entropy[48];for(unsigned j=0;j<48;++j)entropy[j]=(uint8_t)(j+COPA_MODE);copa_nist_seed(entropy);
 uint64_t attempts=0,accepted=0,bs=0,bv=0,over=0,length=0,kgtries=0,keys=0;size_t minlen=COPA_SIG_BYTES,maxlen=0;
 copa_secret_key sk;copa_public_key pub;uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES];
 for(unsigned t=0;t<trials;++t){
  if(t%100==0){copa_keygen_stats kg;assert(!copa_keygen(pk,&sk,&kg));assert(!copa_prepare_public(&pub,pk,sizeof(pk)));kgtries+=kg.attempts;++keys;}
  uint8_t msg[64]={0};for(unsigned j=0;j<4;++j)msg[j]=(uint8_t)(t>>(8*j));copa_sign_stats st;
  assert(!copa_sign(sig,msg,sizeof(msg),&sk,&st));assert(!copa_verify_prepared(sig,sizeof(sig),msg,sizeof(msg),&pub));
  attempts+=st.attempts;accepted+=st.irs_accepted;bs+=st.bs_restarts;bv+=st.bv_restarts;over+=st.encode_overflows;size_t used=32+st.payload_used;length+=used;if(used<minlen)minlen=used;if(used>maxlen)maxlen=used;
 }
 assert(accepted==trials+bs+bv+over);
 printf("{\"mode\":%d,\"signatures\":%u,\"keys\":%"PRIu64",\"keygen_attempts\":%"PRIu64",\"pk_bytes\":%d,\"sig_bytes\":%d,\"attempts\":%"PRIu64",\"irs_accepted\":%"PRIu64",\"irs_aborts\":%"PRIu64",\"bs_restarts\":%"PRIu64",\"bv_restarts\":%"PRIu64",\"encoder_overflows\":%"PRIu64",\"irs_acceptance\":%.12f,\"signature_acceptance\":%.12f,\"mean_iterations\":%.12f,\"mean_variable_sig_bytes\":%.9f,\"min_variable_sig_bytes\":%zu,\"max_variable_sig_bytes\":%zu}\n",COPA_MODE,trials,keys,kgtries,COPA_PK_BYTES,COPA_SIG_BYTES,attempts,accepted,attempts-accepted,bs,bv,over,(double)accepted/attempts,(double)trials/attempts,(double)attempts/trials,(double)length/trials,minlen,maxlen);
 return 0;
}
