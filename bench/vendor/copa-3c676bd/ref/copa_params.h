#ifndef COPA_PARAMS_H
#define COPA_PARAMS_H
#ifndef COPA_MODE
#define COPA_MODE 2
#endif
#define COPA_N 256
#if COPA_MODE == 2
#define COPA_Q 7680
#define COPA_ELL 1
#define COPA_M 2
#define COPA_TAU 58
#define COPA_ALPHA 192
#define COPA_ALPHA_RECIP 87382u
#define COPA_R100 9130
#define COPA_SIGMA100 150
#define COPA_BK100 6960
#define COPA_BK_SQ 4844
#define COPA_BS 3717
#define COPA_BV 3839
#define COPA_PK_BYTES 648
#define COPA_SIG_BYTES 700
#define COPA_T 20
#elif COPA_MODE == 3
#define COPA_Q 12584
#define COPA_ELL 3
#define COPA_M 2
#define COPA_TAU 80
#define COPA_ALPHA 88
#define COPA_ALPHA_RECIP 190651u
#define COPA_R100 3504
#define COPA_SIGMA100 70
#define COPA_BK100 2620
#define COPA_BK_SQ 686
#define COPA_BS 1677
#define COPA_BV 1742
#define COPA_PK_BYTES 840
#define COPA_SIG_BYTES 1063
#define COPA_T 4
#elif COPA_MODE == 5
#define COPA_Q 6048
#define COPA_ELL 4
#define COPA_M 4
#define COPA_TAU 128
#define COPA_ALPHA 168
#define COPA_ALPHA_RECIP 99865u
#define COPA_R100 4466
#define COPA_SIGMA100 70
#define COPA_BK100 3260
#define COPA_BK_SQ 1062
#define COPA_BS 2528
#define COPA_BV 3007
#define COPA_PK_BYTES 1480
#define COPA_SIG_BYTES 1391
#define COPA_T 4
#elif COPA_MODE == 6
#define COPA_Q 8400
#define COPA_ELL 4
#define COPA_M 3
#define COPA_TAU 128
#define COPA_ALPHA 120
#define COPA_ALPHA_RECIP 139811u
#define COPA_R100 4392
#define COPA_SIGMA100 70
#define COPA_BK100 3340
#define COPA_BK_SQ 1115
#define COPA_BS 2367
#define COPA_BV 2530
#define COPA_PK_BYTES 1160
#define COPA_SIG_BYTES 1378
#define COPA_T 12
#else
#error "COPA_MODE must be 2, 3, 5, or 6"
#endif
/* In every mode B_k is approximately sqrt(E||s-hat||^2). A key norm
 * rejection rate near 50% is intentional. Only the upper norm gate is used. */
#define COPA_BLO (30-__builtin_clz((unsigned)(COPA_R100/100)))
#define COPA_KAPPA_NUM COPA_R100
#define COPA_KAPPA_DEN COPA_BK100
_Static_assert(COPA_BK_SQ==COPA_BK100*COPA_BK100/10000,"integer upper norm gate");
#define COPA_HINT_LEVELS (COPA_Q/COPA_ALPHA)
#define COPA_K (1+COPA_ELL+COPA_M)
_Static_assert(COPA_Q < (1<<15),"canonical coefficients fit int16_t");
_Static_assert(1024ULL*COPA_N*(COPA_Q/2)*(COPA_Q/2) < 8380417ULL*7340033ULL,"CRT product has over 1024-fold centered-product margin");
_Static_assert(128ULL*COPA_N*(1+COPA_ELL)*COPA_Q*COPA_Q < 8380417ULL*7340033ULL,"CRT matrix dot products have over 128-fold margin");
_Static_assert(100*(1<<(COPA_BLO+1))<=COPA_R100 && COPA_R100<100*(1<<(COPA_BLO+2)),"b_lo=floor(log2(r))-1");
_Static_assert(COPA_Q%4==0,"q divisible by four");
_Static_assert(COPA_Q%COPA_ALPHA==0,"alpha divides q");
_Static_assert(2*COPA_BV<COPA_Q,"verification norm below q/2");
_Static_assert(COPA_Q%(2*COPA_T)==0,"2t divides q");
#endif
