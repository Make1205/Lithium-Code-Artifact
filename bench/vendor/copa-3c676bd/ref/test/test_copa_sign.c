#include "copa_sign.h"
#include "copa_rans.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc,char **argv)
{
    unsigned trials=argc>1?(unsigned)strtoul(argv[1],NULL,10):100;
    uint64_t minnorm=UINT64_MAX,maxnorm=0;
    uint64_t kg_attempts=0,kg_noninv=0,kg_norm=0,attempts=0,accepted=0,bs=0,bv=0,over=0,used=0;
    for(unsigned t=0;t<trials;++t){
        copa_secret_key sk;uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES],copy[COPA_SIG_BYTES],badpk[COPA_PK_BYTES],msg[64]={0};
        msg[0]=(uint8_t)t;msg[1]=(uint8_t)(t>>8);size_t mlen=t%65;
        copa_keygen_stats kg;copa_sign_stats st;assert(!copa_keygen(pk,&sk,&kg));uint64_t norm=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)norm+=(uint64_t)((int64_t)sk.s[r][i]*sk.s[r][i]);if(norm<minnorm)minnorm=norm;if(norm>maxnorm)maxnorm=norm;assert(!copa_sign(sig,msg,mlen,&sk,&st));
        assert(!copa_verify(sig,sizeof(sig),msg,mlen,pk,sizeof(pk)));
        memcpy(copy,sig,sizeof(sig));copy[0]^=1;assert(copa_verify(copy,sizeof(copy),msg,mlen,pk,sizeof(pk)));
        memcpy(badpk,pk,sizeof(pk));badpk[0]^=1;assert(copa_verify(sig,sizeof(sig),msg,mlen,badpk,sizeof(badpk)));
        uint8_t other[65];memcpy(other,msg,64);other[64]=42;assert(copa_verify(sig,sizeof(sig),other,mlen+1,pk,sizeof(pk)));
        assert(copa_verify(sig,sizeof(sig)-1,msg,mlen,pk,sizeof(pk)));
        int32_t z1[1+COPA_ELL][COPA_N],h[COPA_M][COPA_N];assert(copa_decode(z1,h,sig+32));
        size_t u;memcpy(copy,sig,32);++z1[0][0];
        if(copa_encode(copy+32,&u,z1,h))assert(copa_verify(copy,sizeof(copy),msg,mlen,pk,sizeof(pk)));
        --z1[0][0];h[0][0]=(h[0][0]==COPA_Q/COPA_ALPHA/2)?h[0][0]-1:h[0][0]+1;
        if(copa_encode(copy+32,&u,z1,h))assert(copa_verify(copy,sizeof(copy),msg,mlen,pk,sizeof(pk)));
        memcpy(copy,sig,sizeof(sig));memset(copy+32,0,4);assert(!copa_decode(z1,h,copy+32));
        memcpy(copy,sig,sizeof(sig));copy[35]|=0x80;assert(!copa_decode(z1,h,copy+32));
        if(st.payload_used<COPA_PAYLOAD){memcpy(copy,sig,sizeof(sig));copy[sizeof(copy)-1]=1;assert(!copa_decode(z1,h,copy+32));}
        kg_attempts+=kg.attempts;kg_noninv+=kg.noninvertible;kg_norm+=kg.norm_reject;
        attempts+=st.attempts;accepted+=st.irs_accepted;bs+=st.bs_restarts;bv+=st.bv_restarts;over+=st.encode_overflows;used+=st.payload_used+32;
    }
    assert(trials);
    printf("{\"mode\":%d,\"roundtrips\":%u,\"failures\":0,\"pk_bytes\":%d,\"sig_bytes\":%d,\"min_key_norm_sq\":%"PRIu64",\"max_key_norm_sq\":%"PRIu64",\"keygen_attempts\":%"PRIu64",\"keygen_noninvertible\":%"PRIu64",\"keygen_norm_restarts\":%"PRIu64",\"sign_attempts\":%"PRIu64",\"irs_accepted\":%"PRIu64",\"irs_acceptance\":%.9f,\"bs_restarts\":%"PRIu64",\"bv_restarts\":%"PRIu64",\"encoder_overflows\":%"PRIu64",\"mean_actual_sig_bytes\":%.6f}\n",COPA_MODE,trials,COPA_PK_BYTES,COPA_SIG_BYTES,minnorm,maxnorm,kg_attempts,kg_noninv,kg_norm,attempts,accepted,(double)accepted/attempts,bs,bv,over,(double)used/trials);
    return 0;
}
