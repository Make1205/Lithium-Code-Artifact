/*
 * fips202x8.h -- 8-way SHAKE128/SHAKE256 (AVX-512), declarations.
 *
 * Derived from the XKCP 8-way Keccak-f[1600] code path
 *   https://github.com/XKCP/XKCP .
 * Re-namespaced lithium_fips202x8_avx512_* for this tree (matching the
 * lithium_fips202x4_avx2_* convention of the AVX2 build).
 */
#ifndef LITHIUM_FIPS202X8_H
#define LITHIUM_FIPS202X8_H

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

#define FIPS202X8_NAMESPACE(s) lithium_fips202x8_avx512_##s

typedef struct {
    __m512i s[25];
} keccakx8_state;

#define shake128x8_absorb_once FIPS202X8_NAMESPACE(shake128x8_absorb_once)
void shake128x8_absorb_once(keccakx8_state *state, const uint8_t *in0,
                            const uint8_t *in1, const uint8_t *in2,
                            const uint8_t *in3, const uint8_t *in4,
                            const uint8_t *in5, const uint8_t *in6,
                            const uint8_t *in7, size_t inlen);

#define shake128x8_squeezeblocks \
    FIPS202X8_NAMESPACE(shake128x8_squeezeblocks)
void shake128x8_squeezeblocks(uint8_t *out0, uint8_t *out1, uint8_t *out2,
                              uint8_t *out3, uint8_t *out4, uint8_t *out5,
                              uint8_t *out6, uint8_t *out7, size_t nblocks,
                              keccakx8_state *state);

#define shake256x8_absorb_once FIPS202X8_NAMESPACE(shake256x8_absorb_once)
void shake256x8_absorb_once(keccakx8_state *state, const uint8_t *in0,
                            const uint8_t *in1, const uint8_t *in2,
                            const uint8_t *in3, const uint8_t *in4,
                            const uint8_t *in5, const uint8_t *in6,
                            const uint8_t *in7, size_t inlen);

#define shake256x8_squeezeblocks \
    FIPS202X8_NAMESPACE(shake256x8_squeezeblocks)
void shake256x8_squeezeblocks(uint8_t *out0, uint8_t *out1, uint8_t *out2,
                              uint8_t *out3, uint8_t *out4, uint8_t *out5,
                              uint8_t *out6, uint8_t *out7, size_t nblocks,
                              keccakx8_state *state);

#define shake128x8 FIPS202X8_NAMESPACE(shake128x8)
void shake128x8(uint8_t *out0, uint8_t *out1, uint8_t *out2, uint8_t *out3,
                uint8_t *out4, uint8_t *out5, uint8_t *out6, uint8_t *out7,
                size_t outlen, const uint8_t *in0, const uint8_t *in1,
                const uint8_t *in2, const uint8_t *in3, const uint8_t *in4,
                const uint8_t *in5, const uint8_t *in6, const uint8_t *in7,
                size_t inlen);

#define shake256x8 FIPS202X8_NAMESPACE(shake256x8)
void shake256x8(uint8_t *out0, uint8_t *out1, uint8_t *out2, uint8_t *out3,
                uint8_t *out4, uint8_t *out5, uint8_t *out6, uint8_t *out7,
                size_t outlen, const uint8_t *in0, const uint8_t *in1,
                const uint8_t *in2, const uint8_t *in3, const uint8_t *in4,
                const uint8_t *in5, const uint8_t *in6, const uint8_t *in7,
                size_t inlen);

#endif /* LITHIUM_FIPS202X8_H */
