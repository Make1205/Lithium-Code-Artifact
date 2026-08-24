/* ifma_ref.c -- runnable scalar reference for the Q51 integer ApproxExp
 * (AVX-512 IFMA path).  Fills the 2^{-j/32} table, computes exp(-a) via
 * fixed-point Q51 Horner + table + 2^{-n}, and validates against libm exp.
 * Build & run:  gcc -O2 -o ifma_ref ifma_ref.c -lm && ./ifma_ref
 * Measured max relative error over [0, 7.369]: 2^-48.78.
 * See docs/design-notes/AVXExp.md, Appendix B.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>

// Q51 degree-5 coefficients from tools/run_polyapprox.sh (true value =
// C[i]*2^-51)
static const uint64_t C[6] = {
    2251231781277982ULL, 2256720101499690ULL, 1108248280973626ULL,
    408721085595457ULL,  58973132430121ULL,   37125483605256ULL,
};
static uint64_t T[32];  // T[j] = round(2^{-j/32} * 2^51)

static inline uint64_t mulsh(uint64_t a, uint64_t b, int s)
{
    return (uint64_t)(((unsigned __int128)a * b) >> s);
}

// exp(-a) in the Q51 integer scheme (a in [0, a_max])
static double approx_exp_q51(double a)
{
    const int B = 32, Q = 51;
    const double ln2 = 0.69314718055994530942;
    long m = (long)ceil(B * a / ln2) + (B - 1);
    long n = m >> 5, j = m & 31;
    double rl = ((double)m / B) * ln2 - a;  // r_l in [31/32 ln2, ln2)
    uint64_t X = (uint64_t)llround(ldexp(rl, Q));  // r_l in Q51
    uint64_t p = C[5];
    for (int i = 4; i >= 0; i--)
        p = mulsh(p, X, Q) + C[i];  // Q51 Horner
    uint64_t s = T[j];              // 2^{-j/B} in Q51
    uint64_t ps =
        mulsh(p, s, Q);  // 2^{-j/B} exp(r_l) in Q51, value in [1,2)
    return ldexp((double)ps, -Q - (int)n);  // * 2^{-n}
}

int main(void)
{
    for (int j = 0; j < 32; j++)
        T[j] = (uint64_t)llround(ldexp(pow(2.0, -(double)j / 32.0), 51));
    // validate over the IRS range
    double amax = 7.369, worst = 0;
    double wa = 0;
    for (double a = 0; a <= amax; a += 0.0007) {
        double ref = exp(-a), got = approx_exp_q51(a);
        double rel = fabs((got - ref) / ref);
        if (rel > worst) {
            worst = rel;
            wa = a;
        }
    }
    printf(
        "max relative error over [0,%.3f] = %.3e = 2^%.2f  (at a=%.4f)\n",
        amax, worst, log2(worst), wa);
    printf("\n// T[j] = round(2^{-j/32} * 2^51), Q51, j=0..31\n");
    printf("static const uint64_t TAB_Q51[32] = {\n");
    for (int j = 0; j < 32; j++)
        printf("  %16lluULL,%s", (unsigned long long)T[j],
               (j % 2) ? "\n" : " ");
    printf("};\n\n// double table T_d[j] = 2^{-j/32}\n");
    printf("static const double TAB_D[32] = {\n");
    for (int j = 0; j < 32; j++)
        printf("  %.17g,%s", pow(2.0, -(double)j / 32.0),
               (j % 2) ? "\n" : " ");
    printf("};\n");
    return 0;
}
