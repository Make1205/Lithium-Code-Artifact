/*
 * symmetric.h - Symmetric primitives for IRS.
 *
 * Wraps FIPS 202 SHAKE-256 as the stream cipher / XOF.
 * Uses keccak_state (stack-allocated, no heap).
 */

#ifndef LITHIUM_SYMMETRIC_H
#define LITHIUM_SYMMETRIC_H

#include <stddef.h>
#include <stdint.h>

#include "fips202.h"
#include "params.h"

/* Stream256: SHAKE-256 based XOF */
typedef keccak_state stream256_state;

#define STREAM256_BLOCKBYTES SHAKE256_RATE /* 136 */

/*
 * Initialize SHAKE-256 stream from seed (32 bytes) + nonce (8 bytes LE).
 */
void shake256_stream_init(stream256_state *state,
                          const uint8_t seed[SEEDBYTES], uint64_t nonce);

/* Convenience macros */
#define stream256_init(STATE, SEED, NONCE) \
    shake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, NBLOCKS, STATE) \
    shake256_squeezeblocks(OUT, NBLOCKS, STATE)

#endif /* LITHIUM_SYMMETRIC_H */
