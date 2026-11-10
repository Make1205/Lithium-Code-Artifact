#include "copa_prof.h"
#include "copa_matrix.h"
#include "copa_keygen.h"
#include "ntt_crt.h"
#include "pack_pk.h"
#include "quant.h"
int copa_matrix_expand(copa_public_ntt matrix,const uint8_t pk[COPA_PK_BYTES]){
 uint8_t rho[32];int32_t bp[COPA_M][COPA_N],d[COPA_M][COPA_N],a[COPA_M][COPA_ELL][COPA_N];
 CP_BEGIN(pu);int ok=copa_unpack_pk(rho,&bp[0][0],pk,COPA_PK_BYTES,COPA_MODE);CP_END(CP_PK_UNLAB,pu);if(!ok)return 0;
 CP_BEGIN(pa);copa_expand_a(a,rho);CP_END(CP_A,pa);CP_BEGIN(pd);copa_expand_d(&d[0][0],rho,COPA_MODE);CP_END(CP_D,pd);CP_BEGIN(pntt);
 for(unsigned r=0;r<COPA_M;++r){int32_t first[COPA_N];
  for(unsigned i=0;i<COPA_N;++i){int32_t t=bp[r][i]+d[r][i];first[i]=-(t-COPA_Q*(int32_t)(t>=COPA_Q));}
  CP_BEGIN(pbntt);copa_ntt_forward(matrix[r][0],first,COPA_N);CP_END(CP_B_NTT,pbntt);
  CP_BEGIN(pantt);
  for(unsigned j=0;j<COPA_ELL;++j){copa_ntt_forward(matrix[r][1+j],a[r][j],COPA_N);}CP_END(CP_A_NTT,pantt);
 }CP_END(CP_MATRIX_NTT,pntt);return 1;
}
void copa_matrix_apply_ntt(int32_t out[COPA_M][COPA_N],const copa_public_ntt matrix,const uint32_t input[][2][COPA_N]){
 for(unsigned r=0;r<COPA_M;++r){uint32_t product[2][COPA_N];
  copa_ntt_dot(product,matrix[r],input,1+COPA_ELL,COPA_N);copa_ntt_inverse_modq(out[r],product,COPA_N,COPA_Q);
 }
}
void copa_matrix_apply(int32_t out[COPA_M][COPA_N],const copa_public_ntt matrix,const int32_t input[1+COPA_ELL][COPA_N]){
 uint32_t transformed[1+COPA_ELL][2][COPA_N];for(unsigned j=0;j<1+COPA_ELL;++j)copa_ntt_forward(transformed[j],input[j],COPA_N);
 copa_matrix_apply_ntt(out,matrix,transformed);
}

#ifdef COPA_VERIFY31
#include "../common/copa_verify31.h"
#endif
