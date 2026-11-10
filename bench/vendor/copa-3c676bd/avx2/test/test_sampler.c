/* test_sampler.c - statistical + consistency gate for the Lithium
 * samplers.
 *
 *  (1) D_{Z,r} sampler:  mean ~ 0, stddev ~ r
 *  (2) D_{Z,1} noise:    mean ~ 0, stddev ~ 1
 *  (3) 4x == 1x:         sample_gauss_N_4x lane k == sample_gauss_N (KAT
 * lock)
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "sampler.h"

int main(void)
{
    uint8_t seed[SEEDBYTES];
    for (int i = 0; i < SEEDBYTES; i++)
        seed[i] = (uint8_t)(i * 7 + 1);

    const size_t NB = 60;
    static int32_t buf[N];

    double sum = 0, sumsq = 0;
    long cnt = 0;
    int maxabs = 0;
    for (size_t b = 0; b < NB; b++) {
        sample_gauss_N(buf, seed, (uint64_t)b, N);
        for (int i = 0; i < N; i++) {
            sum += buf[i];
            sumsq += (double)buf[i] * buf[i];
            cnt++;
            int a = buf[i] < 0 ? -buf[i] : buf[i];
            if (a > maxabs)
                maxabs = a;
        }
    }
    double mean = sum / cnt, sd = sqrt(sumsq / cnt - mean * mean);
    int gok = fabs(mean) < 5.0 * R / sqrt((double)cnt) &&
              fabs(sd / (double)R - 1.0) < 0.05;

    double ns = 0, nsq = 0;
    long nc = 0;
    int nmax = 0;
    for (size_t b = 0; b < NB; b++) {
        sample_noise_N(buf, seed, 1000 + (uint64_t)b, N);
        for (int i = 0; i < N; i++) {
            ns += buf[i];
            nsq += (double)buf[i] * buf[i];
            nc++;
            int a = buf[i] < 0 ? -buf[i] : buf[i];
            if (a > nmax)
                nmax = a;
        }
    }
    double nmean = ns / nc, nsd = sqrt(nsq / nc - nmean * nmean);
    int nok = fabs(nmean) < 0.08 && fabs(nsd - 1.0) < 0.06;

    static int32_t a0[N], a1[N], a2[N], a3[N], ref[N];
    sample_gauss_N_4x(a0, a1, a2, a3, seed, 5, 6, 7, 8, N, N, N, N);
    sample_gauss_N(ref, seed, 5, N);
    int kok = (memcmp(a0, ref, sizeof ref) == 0);
    sample_gauss_N(ref, seed, 8, N);
    kok &= (memcmp(a3, ref, sizeof ref) == 0);

    int ok = gok && nok && kok;
    printf(
        "%s test_sampler mode=%d gauss(mean=%.3f sd=%.2f r=%d max=%d) "
        "noise(mean=%.3f sd=%.3f max=%d) 4x==1x=%d\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, mean, sd, R, maxabs, nmean,
        nsd, nmax, kok);
    return ok ? 0 : 1;
}
