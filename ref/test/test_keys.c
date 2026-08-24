/* test_keys.c - NTWE key identity gate: for generated keys, A s == q j
 * (mod 2q).
 *
 * Inline keygen (no packing): sample noise s0,e,f0 ~ D_{Z,1}; f = 2 f0 +
 * 1; invert f (resample if not invertible); b = (A0 s0 + e) f^{-1}; s =
 * [f,s0,e]. Then w = A s mod 2q via mat_mul_2q must equal q j: w[0] =
 * [q,0,...,0],  w[i]=0 for i>0  (all in [0,2q)). Also enforces the keygen
 * norm window B_k' <= ||s|| < B_k and reports accept rate.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "polyvec.h"
#include "randombytes.h"
#include "reduce.h"
#include "sampler.h"

static void noise_poly(poly *p, const uint8_t seed[SEEDBYTES],
                       uint64_t nonce)
{
    int32_t tmp[N];
    sample_noise_N(tmp, seed, nonce, N);
    for (int i = 0; i < N; i++)
        p->coeffs[i] = tmp[i];
}

int main(void)
{
    int good = 0, tries = 0, identity_ok = 1, window_ok = 1, noninv = 0;
    long acc_in_window = 0, total_keys = 0;

    for (int trial = 0; trial < 40 && good < 20; trial++) {
        uint8_t rho[SEEDBYTES], nseed[SEEDBYTES];
        randombytes(rho, SEEDBYTES);
        randombytes(nseed, SEEDBYTES);

        poly s0[ELL], e[M], f0, f, b[M];
        poly16 A0h[M][ELL];
        expand_A0_ntt(A0h, rho);
        for (int j = 0; j < ELL; j++)
            noise_poly(&s0[j], nseed, 100 + j);
        for (int i = 0; i < M; i++)
            noise_poly(&e[i], nseed, 200 + i);
        noise_poly(&f0, nseed, 300);
        for (int k = 0; k < N; k++)
            f.coeffs[k] = 2 * f0.coeffs[k] + (k == 0 ? 1 : 0);

        tries++;
        if (compute_b(b, A0h, s0, e, &f) != 0) {
            noninv++;
            continue;
        }

        /* s = [f, s0, e]; ||s||^2 */
        poly s[K];
        s[0] = f;
        for (int j = 0; j < ELL; j++)
            s[1 + j] = s0[j];
        for (int i = 0; i < M; i++)
            s[1 + ELL + i] = e[i];
        int64_t nrm = 0;
        for (int c = 0; c < K; c++)
            nrm += poly_sqnorm(&s[c]);
        total_keys++;
        int in_win = (nrm >= BK_LOW_SQ && nrm <= BK_SQ);
        if (in_win)
            acc_in_window++;

        /* identity check (independent of the norm window) */
        poly w[M];
        poly16 bh[M];
        for (int i = 0; i < M; i++)
            poly_ntt(&bh[i], &b[i]);
        mat_mul_2q(w, A0h, bh, s);
        for (int i = 0; i < M && identity_ok; i++)
            for (int k = 0; k < N; k++) {
                int32_t want = (i == 0 && k == 0) ? Q : 0;
                if (w[i].coeffs[k] != want) {
                    identity_ok = 0;
                    printf("  mismatch i=%d k=%d got=%d want=%d\n", i, k,
                           w[i].coeffs[k], want);
                    break;
                }
            }
        good++;
    }

    double accrate = total_keys ? (double)acc_in_window / total_keys : 0;
    int ok = identity_ok && window_ok && good > 0;
    printf(
        "%s test_keys mode=%d  identity(As=qj mod2q)=%d  keys=%d "
        "noninv=%d "
        "window_accept=%.0f%% (B_k'^2=%lld..B_k^2=%lld)\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, identity_ok, good, noninv,
        accrate * 100, (long long)BK_LOW_SQ, (long long)BK_SQ);
    return ok ? 0 : 1;
}
