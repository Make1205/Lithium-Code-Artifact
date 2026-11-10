/* test_irs.c - gate for the IRS exp-recurrence rejection sampler.
 *
 *  (A) bit-exact cross-check of decide() vs tools/irs_ref.py over a
 * fixed grid of (nrm, ip, r63): prints an FNV-1a hash of all decision bits
 * and the count of "minus" decisions.  The Python reference prints the
 * same. (B) sign_irs self-consistency: z == y + sum_j (+/-)(x^j s)
 * reconstructed from the recorded sign bits.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "irs.h"
#include "params.h"
#include "polyvec.h"
#include "symmetric.h"

static int64_t gamma_from_nrm(int64_t nrm)
{
    return (int64_t)(((__int128)nrm << 52) / (2 * (int64_t)R * R));
}

static uint64_t st = 0xDEADBEEF12345678ULL;
static uint32_t xr(void)
{
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    return (uint32_t)(st >> 27);
}

int main(void)
{
    /* ---- (A) decision grid cross-check ---- */
    uint64_t h = 1469598103934665603ULL; /* FNV-1a 64 */
    long minus = 0, total = 0;
    const int64_t nrms[3] = {500, 2000, 4000};
    const uint64_t r63s[5] = {1ULL << 10, 1ULL << 40, 1ULL << 62,
                              ((1ULL << 63) - 1), 0x123456789ULL};
    for (int ni = 0; ni < 3; ni++) {
        int64_t g = gamma_from_nrm(nrms[ni]);
        int64_t ctab[REJ_T];
        precompute_c(ctab, g); /* once per gamma */
        for (int64_t ip = -3000; ip <= 3000; ip += 50) {
            for (int ri = 0; ri < 5; ri++) {
                int b = decide(ip, g, ctab, r63s[ri]);
                h ^= (uint64_t)(b & 1);
                h *= 1099511628211ULL;
                minus += b;
                total++;
            }
        }
    }
    printf("GRID mode=%d hash=%llu minus=%ld total=%ld\n", LITHIUM_MODE,
           (unsigned long long)h, minus, total);

    /* ---- (B) sign_irs reconstruction ---- */
    poly sk[K], y[K], z[K], c;
    for (int comp = 0; comp < K; comp++) {
        for (int i = 0; i < N; i++)
            sk[comp].coeffs[i] = (int32_t)(xr() % 19) - 9;
        for (int i = 0; i < N; i++)
            y[comp].coeffs[i] = (int32_t)(xr() % (4 * R)) - 2 * R;
    }
    /* challenge: kappa distinct positions = 1 */
    for (int i = 0; i < N; i++)
        c.coeffs[i] = 0;
    {
        int placed = 0;
        while (placed < KAPPA) {
            int j = xr() % N;
            if (!c.coeffs[j]) {
                c.coeffs[j] = 1;
                placed++;
            }
        }
    }

    memcpy(z, y, sizeof y);
    uint8_t seed[32];
    for (int i = 0; i < 32; i++)
        seed[i] = (uint8_t)(i + 1);
    keccak_state rng;
    shake256_init(&rng);
    shake256_absorb(&rng, seed, 32);
    shake256_finalize(&rng);
    uint8_t signs[(N + 7) / 8];
    memset(signs, 0, sizeof signs);
    sign_irs(z, &c, sk, signs, &rng);

    /* reconstruct y + sum_j sign_j (x^j s) and compare to z */
    poly rec[K];
    memcpy(rec, y, sizeof y);
    unsigned cnt = 0;
    int ok_recon = 1;
    for (unsigned j = 0; j < N; j++) {
        if (!c.coeffs[j])
            continue;
        int sb = (signs[cnt / 8] >> (cnt % 8)) & 1;
        for (int comp = 0; comp < K; comp++) {
            poly v;
            poly_cyclic_shift(&v, &sk[comp], j);
            for (int i = 0; i < N; i++)
                rec[comp].coeffs[i] += sb ? -v.coeffs[i] : v.coeffs[i];
        }
        cnt++;
    }
    for (int comp = 0; comp < K && ok_recon; comp++)
        for (int i = 0; i < N; i++)
            if (rec[comp].coeffs[i] != z[comp].coeffs[i]) {
                ok_recon = 0;
                break;
            }

    printf("%s test_irs mode=%d recon=%d (kappa=%d)\n",
           ok_recon ? "PASS" : "FAIL", LITHIUM_MODE, ok_recon, KAPPA);
    return ok_recon ? 0 : 1;
}
