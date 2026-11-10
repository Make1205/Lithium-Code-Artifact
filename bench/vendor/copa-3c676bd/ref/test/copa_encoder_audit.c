/* Test-only linker wrapper. Production encoder, parameters and bytes are unchanged. */
#include "copa_sign.h"
#include "copa_rans.h"
#include "copa_kat/nist_rng.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "models.h"
#ifdef NDEBUG
#error "The encoder audit requires assertions to execute signing and verification"
#endif
#if COPA_MODE != 5
#error "This experiment measures Copa-5a"
#endif
#define ZCOUNT ((1+COPA_ELL)*COPA_N)
#define HCOUNT (COPA_M*COPA_N)
#define HIST_SIZE 16384
#define POPS 2 /* all encode-eligible candidates; candidates emitted by baseline */
typedef struct {
    uint64_t count,over;
    long double required,ideal,low,quant,flush,state_residual,finite_state;
    uint64_t histogram[HIST_SIZE];
} totals;
static totals stats[POPS][AUD_VARIANTS];
static uint64_t candidates,baseline_over;
static long double max_identity_error;
int __real_copa_encode(uint8_t *,size_t *,const int32_t [1+COPA_ELL][COPA_N],const int32_t [COPA_M][COPA_N]);

static size_t dry_run(const aud_model *m,const unsigned *high,const unsigned *hint,uint32_t *final)
{
    uint64_t lower=UINT64_C(1)<<m->lower_bits,state=lower;size_t bytes=0;
    for(int i=ZCOUNT+HCOUNT-1;i>=0;--i){
        const uint32_t *cdf=i>=ZCOUNT?m->hint:m->high;
        unsigned s=i>=ZCOUNT?hint[i-ZCOUNT]:high[i];
        uint64_t freq=cdf[s+1]-cdf[s];
        uint64_t maximum=((lower>>m->precision)<<8)*freq;
        while(state>=maximum){++bytes;state>>=8;}
        state=((state/freq)<<m->precision)+(state%freq)+cdf[s];
        assert(state>=lower&&state<(lower<<8)&&state<=UINT32_MAX);
    }
    *final=(uint32_t)state;
    return 6+ZCOUNT*COPA_BLO/8+bytes;
}

int __wrap_copa_encode(uint8_t out[COPA_PAYLOAD],size_t *used,const int32_t z1[1+COPA_ELL][COPA_N],const int32_t h[COPA_M][COPA_N])
{
    unsigned zh[ZCOUNT],hh[HCOUNT];int64_t square=0;
    long double high_ideal=0,hint_ideal=0,quant[AUD_VARIANTS]={0};
    for(unsigned i=0;i<ZCOUNT;++i){
        int32_t v=(&z1[0][0])[i];assert(v>=-COPA_BS&&v<=COPA_BS);
        int high=v>=0?v/32:-((-v+31)/32),index=high-AUD_HIGH_MIN;
        assert(index>=0&&index<AUD_HIGH_COUNT);zh[i]=(unsigned)index;
        square+=(int64_t)v*v;high_ideal+=aud_high_ideal_bits[index];
        for(unsigned m=0;m<AUD_VARIANTS;++m)quant[m]+=aud_models[m].high_bits[index];
    }
    for(unsigned i=0;i<HCOUNT;++i){
        int index=(&h[0][0])[i]-AUD_HINT_MIN;assert(index>=0&&index<AUD_HINT_COUNT);hh[i]=(unsigned)index;
        hint_ideal+=aud_hint_ideal_bits[index];
        for(unsigned m=0;m<AUD_VARIANTS;++m)quant[m]+=aud_models[m].hint_bits[index];
    }
    long double z_ideal=ZCOUNT*AUD_Z_LOGZ+square*AUD_Z_SQUARE_SCALE;
    long double ideal=(z_ideal+hint_ideal)/8+32;
    long double low=(ZCOUNT*COPA_BLO+high_ideal-z_ideal)/8;
    size_t lengths[AUD_VARIANTS];uint32_t states[AUD_VARIANTS];
    for(unsigned m=0;m<AUD_VARIANTS;++m)lengths[m]=dry_run(&aud_models[m],zh,hh,&states[m]);
    size_t real_used=SIZE_MAX;int ok=__real_copa_encode(out,&real_used,z1,h);
    assert(ok==(lengths[0]<=COPA_PAYLOAD));
    if(ok){
        assert(real_used==lengths[0]);
        uint32_t actual=0;for(unsigned j=0;j<4;++j)actual|=(uint32_t)out[2+j]<<(8*j);
        assert(actual==states[0]);if(used)*used=real_used;
    }
    ++candidates;baseline_over+=!ok;
    for(unsigned m=0;m<AUD_VARIANTS;++m){
        long double qloss=(quant[m]-high_ideal-hint_ideal)/8;
        long double flush=(long double)lengths[m]-(ZCOUNT*COPA_BLO+quant[m])/8;
        long double residual=((long double)aud_models[m].lower_bits-log2l(states[m]))/8;
        long double finite=flush-6-residual;
        long double err=fabsl(low+qloss+flush-((long double)lengths[m]+32-ideal));
        assert(err<1e-8L);if(err>max_identity_error)max_identity_error=err;
        assert(lengths[m]<HIST_SIZE);
        for(unsigned pop=0;pop<POPS;++pop){
            if(pop==1&&!ok)continue;
            totals *t=&stats[pop][m];++t->count;t->over+=lengths[m]>COPA_PAYLOAD;
            t->required+=lengths[m]+32;t->ideal+=ideal;t->low+=low;t->quant+=qloss;t->flush+=flush;
            t->state_residual+=residual;t->finite_state+=finite;++t->histogram[lengths[m]];
        }
    }
    return ok;
}

static void selftest(void)
{
    int32_t z[1+COPA_ELL][COPA_N]={{0}},h[COPA_M][COPA_N]={{0}};uint8_t out[COPA_PAYLOAD];size_t used;
    assert(__wrap_copa_encode(out,&used,z,h));
    for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i)h[r][i]=20;
    assert(!__wrap_copa_encode(out,&used,z,h));
    memset(stats,0,sizeof(stats));candidates=baseline_over=0;max_identity_error=0;
}

int main(int argc,char **argv)
{
    uint64_t runs=argc>1?strtoull(argv[1],NULL,10):20;
    unsigned worker=argc>2?(unsigned)strtoul(argv[2],NULL,10):0;
    assert(runs&&worker<256);selftest();copa_nist_selftest();
    uint8_t seed[48];for(unsigned j=0;j<48;++j)seed[j]=(uint8_t)(j^(31*worker));copa_nist_seed(seed);
    copa_secret_key sk;uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES];
    uint64_t keys=0,attempts=0,irs=0,bs=0,bv=0,over=0;
    for(uint64_t i=0;i<runs;++i){
        if(i%100==0){assert(!copa_keygen(pk,&sk,NULL));++keys;}
        uint8_t msg[32]={0};for(unsigned j=0;j<8;++j)msg[j]=(uint8_t)(i>>(8*j));msg[8]=(uint8_t)worker;msg[31]=0xa5;
        copa_sign_stats s;assert(!copa_sign(sig,msg,sizeof(msg),&sk,&s));assert(!copa_verify(sig,sizeof(sig),msg,sizeof(msg),pk,sizeof(pk)));
        attempts+=s.attempts;irs+=s.irs_accepted;bs+=s.bs_restarts;bv+=s.bv_restarts;over+=s.encode_overflows;
        if((i+1)%1000==0)fprintf(stderr,"worker=%u verified=%"PRIu64" overflows=%"PRIu64"\n",worker,i+1,over);
    }
    assert(candidates==runs+over&&over==baseline_over&&stats[1][0].count==runs);
    printf("{\"worker\":%u,\"signatures\":%"PRIu64",\"verified\":%"PRIu64",\"keys\":%"PRIu64",\"sign_attempts\":%"PRIu64",\"irs_accepted\":%"PRIu64",\"bs_restarts\":%"PRIu64",\"bv_restarts\":%"PRIu64",\"encode_candidates\":%"PRIu64",\"baseline_overflows\":%"PRIu64",\"max_decomposition_identity_error_bytes\":%.15Lg,\"expected_ideal_signature_bytes\":%.12Lf,\"models\":[",worker,runs,runs,keys,attempts,irs,bs,bv,candidates,over,max_identity_error,AUD_EXPECTED_IDEAL_SIG_BYTES);
    for(unsigned m=0;m<AUD_VARIANTS;++m){
        if(m)putchar(',');
        printf("{\"precision_bits\":%u,\"lower_state_bits\":%u,\"populations\":[",aud_models[m].precision,aud_models[m].lower_bits);
        for(unsigned pop=0;pop<POPS;++pop){
            totals *t=&stats[pop][m];if(pop)putchar(',');
            printf("{\"name\":\"%s\",\"count\":%"PRIu64",\"overflows\":%"PRIu64",\"sum_required_signature_bytes\":%.12Lf,\"sum_ideal_signature_bytes\":%.12Lf,\"sum_low_overhead_bytes\":%.12Lf,\"sum_frequency_overhead_bytes\":%.12Lf,\"sum_flush_overhead_bytes\":%.12Lf,\"sum_state_residual_bytes\":%.12Lf,\"sum_finite_state_overhead_bytes\":%.12Lf,\"payload_length_histogram\":{",pop?"baseline_emitted":"all_candidates",t->count,t->over,t->required,t->ideal,t->low,t->quant,t->flush,t->state_residual,t->finite_state);
            unsigned first=1;for(unsigned len=0;len<HIST_SIZE;++len)if(t->histogram[len]){if(!first)putchar(',');first=0;printf("\"%u\":%"PRIu64,len,t->histogram[len]);}
            printf("}}");
        }
        printf("]}");
    }
    printf("]}\n");return 0;
}
