#include "copa_irs64.h"
void probabilities64(uint64_t out[4],uint64_t norm,int64_t dot){
 static copa_irs64_context ctx;static uint64_t last=0;
 if(norm!=last){copa_irs64_prepare(&ctx,norm);last=norm;}
 copa_fp63 a,s;copa_irs64_probabilities(&a,&s,dot,&ctx);
 out[0]=a.m;out[1]=(uint32_t)a.e;out[2]=s.m;out[3]=(uint32_t)s.e;
}
