/* test_keccakx8.c -- verify the 8-way SHAKE matches scalar SHAKE per lane.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "fips202x8.h"

#define OUTLEN 400 /* spans multiple rate blocks for both 128 and 256 */

static int run128(void)
{
    uint8_t in[8][37], o8[8][OUTLEN], o1[OUTLEN];
    int fails = 0;
    for (int j = 0; j < 8; j++)
        for (int i = 0; i < 37; i++)
            in[j][i] = (uint8_t)(j * 41 + i * 7 + 3);
    shake128x8(o8[0], o8[1], o8[2], o8[3], o8[4], o8[5], o8[6], o8[7],
               OUTLEN, in[0], in[1], in[2], in[3], in[4], in[5], in[6],
               in[7], 37);
    for (int j = 0; j < 8; j++) {
        shake128(o1, OUTLEN, in[j], 37);
        if (memcmp(o1, o8[j], OUTLEN)) {
            fails++;
            printf("  shake128 lane %d MISMATCH\n", j);
        }
    }
    return fails;
}

static int run256(void)
{
    uint8_t in[8][49], o8[8][OUTLEN], o1[OUTLEN];
    int fails = 0;
    for (int j = 0; j < 8; j++)
        for (int i = 0; i < 49; i++)
            in[j][i] = (uint8_t)(j * 13 + i * 5 + 1);
    shake256x8(o8[0], o8[1], o8[2], o8[3], o8[4], o8[5], o8[6], o8[7],
               OUTLEN, in[0], in[1], in[2], in[3], in[4], in[5], in[6],
               in[7], 49);
    for (int j = 0; j < 8; j++) {
        shake256(o1, OUTLEN, in[j], 49);
        if (memcmp(o1, o8[j], OUTLEN)) {
            fails++;
            printf("  shake256 lane %d MISMATCH\n", j);
        }
    }
    return fails;
}

int main(void)
{
    int f = run128() + run256();
    printf("test_keccakx8: %s (8x==1x, shake128+shake256, %d fails)\n",
           f ? "FAIL" : "PASS", f);
    return f ? 1 : 0;
}
