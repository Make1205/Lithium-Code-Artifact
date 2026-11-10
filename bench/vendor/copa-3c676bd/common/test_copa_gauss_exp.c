#include "copa_scalar64.h"
#include "copa_scalar64_tables.h"
#include "copa_sampler_tables.h"
#include "copa_gauss_simd.h"
#include <assert.h>
#include <stdio.h>
int main(void){uint64_t nums[CG_LANES],values[CG_LANES];unsigned cases=0;
 for(unsigned x=0;x<=COPA_BASE_TAIL;++x)for(unsigned y=0;y<(1u<<COPA_BASE_BITS);++y){for(unsigned j=0;j<CG_LANES;++j){unsigned v=(x*(1u<<COPA_BASE_BITS)+y+j)%((COPA_BASE_TAIL+1)*(1u<<COPA_BASE_BITS)),xx=v>>COPA_BASE_BITS,yy=v&((1u<<COPA_BASE_BITS)-1);nums[j]=(uint64_t)yy*(yy+(2u<<COPA_BASE_BITS)*xx)*10000;}cv_store(values,copa_gauss_exp_vector(cv_load(nums)));for(unsigned j=0;j<CG_LANES;++j)if(values[j]!=copa_gaussian_threshold(nums[j])){printf("FAIL x=%u low=%u num=%llu got=%llu expected=%llu reciprocal=%llu log=%llu shift=%u\n",x,y,(unsigned long long)nums[j],(unsigned long long)values[j],(unsigned long long)copa_gaussian_threshold(nums[j]),(unsigned long long)copa64_r2.reciprocal,(unsigned long long)copa64_r2.log_reciprocal,copa64_r2.shift);return 1;};++cases;}
 printf("PASS mode=%d Gaussian vector thresholds: %u complete proposal pairs\n",COPA_MODE,cases);return 0;}
