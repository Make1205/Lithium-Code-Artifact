#include "copa_sampler_tables.h"
#include "copa_gaussian_fp.h"
uint64_t threshold56(uint64_t num){copa_gaussian_env old=copa_gaussian_enter();uint64_t p=copa_gaussian_threshold56(num);copa_gaussian_leave(old);return p;}
unsigned table_limb(unsigned row,unsigned limb){return copa_base_rcdt[row][limb];}
unsigned table_width(void){return COPA_BASE_WIDTH;}
unsigned table_tail(void){return COPA_BASE_TAIL;}
unsigned sampler_k(void){return 1u<<COPA_BASE_BITS;}
unsigned sampler_r100(void){return COPA_R100;}
