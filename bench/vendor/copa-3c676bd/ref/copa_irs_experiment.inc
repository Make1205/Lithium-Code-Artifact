/* Fixed-schedule incremental/full-scan controls, opt-in only. */
static int irs_run(int32_t z[COPA_K][COPA_N],const uint8_t challenge[COPA_N],
             const copa_secret_key *sk,const uint32_t z_ntt[COPA_K][2][COPA_N],
             const uint8_t seed[32],uint64_t attempt,unsigned *steps
#ifdef COPA_USE_IRS_POOL
,copa_irs_pool *pool
#endif
)
{
#if !defined(COPA_IRS_EXP) || COPA_IRS_EXP!=2
    SP_BEGIN(pdot);uint32_t correlation[2][COPA_N];int64_t projection[COPA_N];
    copa_ntt_dot_reversed(correlation,z_ntt,sk->secret_ntt,COPA_K,COPA_N);copa_ntt_inverse_crt(projection,correlation,COPA_N);SP_END(SP_IRS_DOT,pdot);
#endif
    unsigned position[COPA_TAU],count=0;int32_t masks[COPA_TAU];
    for(unsigned j=0;j<COPA_N;++j)if(challenge[j])position[count++]=j;
    uint8_t ext[41];memcpy(ext,seed,32);ext[32]=0x15;for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));
    copa_xofbuf rng;
#ifdef COPA_USE_IRS_POOL
    if(pool){uint64_t group=attempt/COPA_IRS_POOL_LANES;if(!pool->ready||pool->group!=group){uint8_t inputs[COPA_IRS_POOL_LANES][41];const uint8_t *ptr[COPA_IRS_POOL_LANES];for(unsigned lane=0;lane<COPA_IRS_POOL_LANES;++lane){memcpy(inputs[lane],ext,33);uint64_t nonce=group*COPA_IRS_POOL_LANES+lane;for(unsigned j=0;j<8;++j)inputs[lane][33+j]=(uint8_t)(nonce>>(8*j));ptr[lane]=inputs[lane];}
#if COPA_IRS_POOL_LANES==4
     copa_xof4_init(&pool->engine,ptr,41);
#else
     copa_xof8_init(&pool->engine,ptr,41);
#endif
     pool->ready=1;pool->group=group;pool->filled=0;}}
    else
#endif
    copa_xof_init(&rng,ext,sizeof(ext));
#if defined(COPA_IRS_EXP) && (COPA_IRS_EXP<1 || COPA_IRS_EXP>2)
#error COPA_IRS_EXP must be 1 (incremental fixed schedule) or 2 (full scan fixed schedule)
#endif
#ifdef COPA_IRS_EXP
    unsigned alive=1;*steps=0;
#endif
    for(unsigned k=0;k<count;++k){
        unsigned j=position[k];uint8_t rb[48];uint64_t words[6];
        SP_BEGIN(prng);
        /* Step k is the fixed slice [48*k,48*(k+1)) of this attempt's
         * domain-15 stream; another attempt never shares its XOF state. */
#ifdef COPA_USE_IRS_POOL
        if(pool){while(pool->filled<48*(k+1)){uint8_t blocks[COPA_IRS_POOL_LANES][136];
#if COPA_IRS_POOL_LANES==4
         copa_xof4_block(&pool->engine,blocks);
#else
         copa_xof8_block(&pool->engine,blocks);
#endif
         for(unsigned lane=0;lane<COPA_IRS_POOL_LANES;++lane)memcpy(pool->data[lane]+pool->filled,blocks[lane],136);pool->filled+=136;}
         memcpy(rb,pool->data[attempt%COPA_IRS_POOL_LANES]+48*k,48);
        }else
#endif
        copa_xof_read(&rng,rb,sizeof(rb));for(unsigned b=0;b<6;++b)words[b]=le64(rb+8*b);
        SP_END(SP_IRS_SHAKE,prng);
#if defined(COPA_IRS_EXP) && COPA_IRS_EXP==2
        SP_BEGIN(pfull);int64_t dot=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){int32_t v=i>=j?sk->s[r][i-j]:-sk->s[r][i+COPA_N-j];dot+=(int64_t)z[r][i]*v;}SP_END(SP_IRS_DOT,pfull);
#else
        int64_t dot=projection[j];
#endif
        SP_BEGIN(pprob);unsigned minus,accepted=copa_irs64_step(&minus,dot,&sk->irs_context,words);SP_END(SP_IRS_PROB,pprob);
#ifdef COPA_IRS_EXP
        /* Opaque results keep every probability evaluation live even after
         * the first rejection. The loop trip count stays public. */
        __asm__ volatile("":"+r"(accepted),"+r"(minus)::"memory");
        *steps+=alive;alive&=accepted;alive=(unsigned)(copa64_ct_mask(alive)&1);
#else
        *steps=k+1;if(!accepted)return 0;
#endif
         /* permitted step-boundary early rejection */
        masks[k]=-(int32_t)minus;int64_t mask=-(int64_t)minus;
        SP_BEGIN(pupdate);
#if defined(COPA_IRS_EXP) && COPA_IRS_EXP==2
        for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){int32_t v=i>=j?sk->s[r][i-j]:-sk->s[r][i+COPA_N-j];z[r][i]+=(v^masks[k])-masks[k];}
#else
        for(unsigned next=k+1;next<count;++next){unsigned target=position[next];int64_t c=sk->autocorrelation[target-j];projection[target]+=(c^mask)-mask;}
#endif
        SP_END(SP_IRS_UPDATE,pupdate);
    }
#ifdef COPA_IRS_EXP
#ifdef COPA_IRS_CT_CHECK
    /* Only the final public success/failure status is declassified. All
     * loop operands, epsilon masks and projections remain tainted. */
    VALGRIND_MAKE_MEM_DEFINED(&alive,sizeof(alive));
#endif
    if(!alive)return 0;
#endif
    SP_BEGIN(pmul);
#if !defined(COPA_IRS_EXP) || COPA_IRS_EXP!=2
    /* All updates commute. Rejected candidates need no vector writes. */
    for(unsigned r=0;r<COPA_K;++r)for(unsigned k=0;k<count;++k){unsigned j=position[k];int32_t mask=masks[k];
#ifdef COPA_VEC_POLY_HEADER
        copa_backend_rotate_add(z[r],sk->s[r],j,mask);
#else
        for(unsigned i=j;i<COPA_N;++i)z[r][i]+=(sk->s[r][i-j]^mask)-mask;
        for(unsigned i=0;i<j;++i)z[r][i]-=(sk->s[r][i+COPA_N-j]^mask)-mask;
#endif
    }
    #endif
    SP_END(SP_IRS_MUL,pmul);return 1;
}
