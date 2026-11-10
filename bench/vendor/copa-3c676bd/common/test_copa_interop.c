/* Portable corpus with explicit headers, deterministic NIST DRBG and fixed
 * messages. Emit/verify binaries are built separately with each backend.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "copa_kat/nist_rng.h"
#include "copa_sign.h"
static void message(uint8_t m[64], unsigned key, unsigned index)
{
    for (unsigned i = 0; i < 64; ++i)
        m[i] = (uint8_t)(key * 37 + index * 11 + i * 19);
}
int main(int argc, char **argv)
{
    assert(argc == 3 || argc == 4);
    unsigned chunk = argc == 4 ? (unsigned)strtoul(argv[3], NULL, 10) : 0;
    unsigned keys = argc == 4 ? 10 : 100;
    assert(chunk < 256);
    int emit = !strcmp(argv[1], "emit");
    assert(emit || !strcmp(argv[1], "verify"));
    FILE *f = fopen(argv[2], emit ? "wb" : "rb");
    assert(f);
    uint8_t seed[48];
    for (unsigned i = 0; i < 48; ++i)
        seed[i] = (uint8_t)i;
    if (argc == 4) {
        seed[0] ^= (uint8_t)chunk;
        seed[1] ^= 0xa5;
    }
    copa_nist_seed(seed);
    uint8_t hdr[9] = {'C',  'O', 'P',  'A', COPA_MODE, argc == 4 ? 2 : 1,
                      keys, 100, chunk},
            readhdr[9];
    size_t header_bytes = argc == 4 ? 9 : 8;
    if (emit)
        assert(fwrite(hdr, 1, header_bytes, f) == header_bytes);
    else {
        assert(fread(readhdr, 1, header_bytes, f) == header_bytes);
        assert(!memcmp(hdr, readhdr, header_bytes));
    }
    for (unsigned key = 0; key < keys; ++key) {
        copa_secret_key sk;
        uint8_t pk[COPA_PK_BYTES], sig[COPA_SIG_BYTES], msg[64];
        if (emit) {
            assert(!copa_keygen(pk, &sk, NULL));
            assert(fwrite(pk, 1, sizeof(pk), f) == sizeof(pk));
        } else
            assert(fread(pk, 1, sizeof(pk), f) == sizeof(pk));
        for (unsigned i = 0; i < 100; ++i) {
            message(msg, key, i);
            size_t len = i % 65;
            if (emit) {
                assert(!copa_sign(sig, msg, len, &sk, NULL));
                assert(fwrite(sig, 1, sizeof(sig), f) == sizeof(sig));
            } else
                assert(fread(sig, 1, sizeof(sig), f) == sizeof(sig));
            assert(
                !copa_verify(sig, sizeof(sig), msg, len, pk, sizeof(pk)));
        }
    }
    if (!emit)
        assert(fgetc(f) == EOF);
    assert(!fclose(f));
    printf("PASS mode=%d %s corpus: %u keys, %u signatures, chunk=%u\n",
           COPA_MODE, argv[1], keys, 100 * keys, chunk);
    return 0;
}
