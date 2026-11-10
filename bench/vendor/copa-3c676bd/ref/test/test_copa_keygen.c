#include "copa_keygen.h"
#include "pack_pk.h"
#include "quant.h"
#include "ntt_crt.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
static void add_product(int64_t acc[COPA_N],const int32_t a[COPA_N],const int32_t b[COPA_N],int sign)
{for(unsigned i=0;i<COPA_N;++i)for(unsigned j=0;j<COPA_N;++j){int64_t v=sign*(int64_t)a[i]*b[j];if(i+j<COPA_N)acc[i+j]+=v;else acc[i+j-COPA_N]-=v;}}
int main(void)
{
    uint64_t attempts=0,noninv=0,normreject=0;
    for(unsigned trial=0;trial<100;++trial){
        copa_secret_key sk;uint8_t pk[COPA_PK_BYTES],rho[32];copa_keygen_stats stats;
        int32_t b[COPA_M][COPA_N],d[COPA_M][COPA_N],a[COPA_M][COPA_ELL][COPA_N];
        assert(!copa_keygen(pk,&sk,&stats));assert(copa_unpack_pk(rho,&b[0][0],pk,sizeof(pk),COPA_MODE));
        copa_expand_a(a,rho);copa_expand_d(&d[0][0],rho,COPA_MODE);
        for(unsigned r=0;r<COPA_M;++r){
            int32_t bp[COPA_N];int64_t acc[COPA_N]={0};
            for(unsigned i=0;i<COPA_N;++i)bp[i]=(int32_t)copa_mod_i64((int64_t)b[r][i]+d[r][i],COPA_Q);
            add_product(acc,bp,sk.s[0],-1);
            for(unsigned j=0;j<COPA_ELL;++j)add_product(acc,a[r][j],sk.s[1+j],1);
            for(unsigned i=0;i<COPA_N;++i)assert(copa_mod_i64(acc[i]+sk.s[1+COPA_ELL+r][i],COPA_Q)==(r==0&&i==0?COPA_Q/2:0));
        }
        int64_t norm=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)norm+=(int64_t)sk.s[r][i]*sk.s[r][i];
        assert(norm<=COPA_BK_SQ);assert(stats.attempts==stats.noninvertible+stats.norm_reject+1);
        attempts+=stats.attempts;noninv+=stats.noninvertible;normreject+=stats.norm_reject;
    }
    printf("PASS Copa-%d KeyGen: keys=100 pk_bytes=%d attempts=%"PRIu64" noninvertible=%"PRIu64" norm_reject=%"PRIu64" key_relation=all norm_gate=all\n",COPA_MODE,COPA_PK_BYTES,attempts,noninv,normreject);
    return 0;
}
