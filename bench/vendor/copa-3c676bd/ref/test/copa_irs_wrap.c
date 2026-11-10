#include "copa_q120_sampler.h"
void thresholds(uint64_t out[4],uint64_t norm,int64_t dot)
{
    static copa_i128 table[COPA_IRS_TERMS];static uint64_t previous=UINT64_MAX;
    if(norm!=previous){copa_irs_table(table,norm);previous=norm;}
    copa_i128 a,ab;copa_irs_thresholds(&a,&ab,dot,norm,table);
    out[0]=(uint64_t)a;out[1]=(uint64_t)((copa_u128)a>>64);out[2]=(uint64_t)ab;out[3]=(uint64_t)((copa_u128)ab>>64);
}

void get_m(uint64_t out[2]){out[0]=(uint64_t)COPA_IRS_M;out[1]=(uint64_t)(COPA_IRS_M>>64);}
