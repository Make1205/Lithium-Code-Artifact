/* test_rans.c - gate for irs_rans.c.
 *
 *  (1) roundtrip: decode(encode(zhead,hint)) == (zhead,hint).
 *  (2) canonical stream: trailing or missing bytes are rejected.
 *  (3) deterministic dump: fixed input -> FNV-1a hash + length, reproduced
 *      byte-for-byte by tools/check_rans.py (C<->Python lock).
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "irs_rans.h"
#include "params.h"

static uint64_t st = 0x243F6A8885A308D3ULL;
static uint32_t xr(void)
{
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    return (uint32_t)(st >> 27);
}

int main(void)
{
    enum { NZ = 64, NH = 48 };
    int32_t zh[NZ], hi[NH], zh2[NZ], hi2[NH];
    uint8_t buf[8192];
    size_t len;

    int rt_fail = 0;
    int canon_fail = 0;
    for (int t = 0; t < 3000; t++) {
        for (int i = 0; i < NZ; i++)
            zh[i] = RANS_ZH_LO + (int)(xr() % RANS_ZH_N);
        for (int i = 0; i < NH; i++)
            hi[i] = RANS_HINT_LO + (int)(xr() % RANS_HINT_N);
        if (rans_encode(buf, &len, sizeof buf, zh, NZ, hi, NH) != 0) {
            rt_fail++;
            continue;
        }
        if (rans_decode(zh2, NZ, hi2, NH, buf, len) != 0) {
            rt_fail++;
            continue;
        }
        if (memcmp(zh, zh2, sizeof zh) || memcmp(hi, hi2, sizeof hi))
            rt_fail++;
        if (len < sizeof buf) {
            buf[len] = 0x42;
            if (rans_decode(zh2, NZ, hi2, NH, buf, len + 1) == 0)
                canon_fail++;
        }
        if (len > 0 && rans_decode(zh2, NZ, hi2, NH, buf, len - 1) == 0)
            canon_fail++;
    }

    st = 0xCAFEBABEDEADBEEFULL;
    for (int i = 0; i < NZ; i++)
        zh[i] = RANS_ZH_LO + (int)(xr() % RANS_ZH_N);
    for (int i = 0; i < NH; i++)
        hi[i] = RANS_HINT_LO + (int)(xr() % RANS_HINT_N);
    rans_encode(buf, &len, sizeof buf, zh, NZ, hi, NH);
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < len; i++) {
        h ^= buf[i];
        h *= 1099511628211ULL;
    }

    int ok = (rt_fail == 0) && (canon_fail == 0);
    printf(
        "%s test_rans mode=%d ZH_N=%d HINT_N=%d roundtrip_fail=%d "
        "canonical_fail=%d DUMP len=%zu hash=%llu\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, RANS_ZH_N, RANS_HINT_N,
        rt_fail, canon_fail, len, (unsigned long long)h);
    return ok ? 0 : 1;
}
