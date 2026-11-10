#ifndef COPA_PACK_PK_H
#define COPA_PACK_PK_H
#include <stddef.h>
#include <stdint.h>
size_t copa_pk_bytes(unsigned mode);
/* b is row-major coefficient-domain b in [0,q). */
int copa_pack_pk(uint8_t *pk,size_t len,const uint8_t rho[32],const int32_t *b,unsigned mode);
int copa_unpack_pk(uint8_t rho[32],int32_t *b,const uint8_t *pk,size_t len,unsigned mode);
#endif
