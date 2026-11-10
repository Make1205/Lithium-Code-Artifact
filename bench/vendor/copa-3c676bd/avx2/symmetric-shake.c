/*
 * symmetric-shake.c - SHAKE-256 stream initialization for IRS.
 */

#include <string.h>

#include "symmetric.h"

void shake256_stream_init(stream256_state *state,
                          const uint8_t seed[SEEDBYTES], uint64_t nonce)
{
    uint8_t buf[SEEDBYTES + 8];
    memcpy(buf, seed, SEEDBYTES);
    buf[SEEDBYTES + 0] = (uint8_t)(nonce);
    buf[SEEDBYTES + 1] = (uint8_t)(nonce >> 8);
    buf[SEEDBYTES + 2] = (uint8_t)(nonce >> 16);
    buf[SEEDBYTES + 3] = (uint8_t)(nonce >> 24);
    buf[SEEDBYTES + 4] = (uint8_t)(nonce >> 32);
    buf[SEEDBYTES + 5] = (uint8_t)(nonce >> 40);
    buf[SEEDBYTES + 6] = (uint8_t)(nonce >> 48);
    buf[SEEDBYTES + 7] = (uint8_t)(nonce >> 56);
    shake256_absorb_once(state, buf, sizeof(buf));
}
