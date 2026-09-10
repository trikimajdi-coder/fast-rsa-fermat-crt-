#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <gmp.h>

typedef struct {
    unsigned long long iterations;
    double execution_time;
    int success;
} CRTResult;

CRTResult fast_crt_factorize_31_omp(
    const char *N_str,
    int *r0, int len0,
    int *r1, int len1,
    int *r2, int len2,
    int *r3, int len3,
    int *r4, int len4,
    int *r5, int len5,
    int *r6, int len6,
    int *r7, int len7,
    unsigned long long *W,
    unsigned long long M,
    int k_max,
    char *out_x_str
) {
    CRTResult res = {0, 0.0, 0};
    double start_time = omp_get_wtime();

    mpz_t N, sqrtN;
    mpz_init(N);
    mpz_init(sqrtN);

    mpz_set_str(N, N_str, 10);
    mpz_sqrt(sqrtN, N);

    volatile int found = 0;
    unsigned long long total_iter = 0;

    #pragma omp parallel reduction(+:total_iter)
    {
        mpz_t t_X, t_Y_sq;
        mpz_init(t_X);
        mpz_init(t_Y_sq);

        unsigned long long local_iter = 0;

        #pragma omp for collapse(2) schedule(dynamic, 1)
        for (int i0 = 0; i0 < len0; i0++) {
            for (int i1 = 0; i1 < len1; i1++) {
                if (found) continue;
                for (int i2 = 0; i2 < len2; i2++) {
                    if (found) continue;
                    for (int i3 = 0; i3 < len3; i3++) {
                        if (found) continue;
                        for (int i4 = 0; i4 < len4; i4++) {
                            if (found) continue;
                            for (int i5 = 0; i5 < len5; i5++) {
                                if (found) continue;
                                for (int i6 = 0; i6 < len6; i6++) {
                                    if (found) continue;
                                    for (int i7 = 0; i7 < len7; i7++) {
                                        if (found) continue;

                                        unsigned long long rem = (
                                            (unsigned long long)r0[i0] * W[0] +
                                            (unsigned long long)r1[i1] * W[1] +
                                            (unsigned long long)r2[i2] * W[2] +
                                            (unsigned long long)r3[i3] * W[3] +
                                            (unsigned long long)r4[i4] * W[4] +
                                            (unsigned long long)r5[i5] * W[5] +
                                            (unsigned long long)r6[i6] * W[6] +
                                            (unsigned long long)r7[i7] * W[7]
                                        ) % M;

                                        for (int k = 0; k < k_max; k++) {
                                            local_iter++;
                                            unsigned long long x_val = rem + (unsigned long long)k * M;

                                            mpz_set_ui(t_X, x_val);
                                            mpz_add(t_X, sqrtN, t_X);

                                            mpz_mul(t_Y_sq, t_X, t_X);
                                            mpz_sub(t_Y_sq, t_Y_sq, N);

                                            if (mpz_sgn(t_Y_sq) >= 0 && mpz_perfect_square_p(t_Y_sq)) {
                                                #pragma omp critical
                                                {
                                                    if (!found) {
                                                        found = 1;
                                                        mpz_get_str(out_x_str, 10, t_X);
                                                    }
                                                }
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        total_iter += local_iter;

        mpz_clear(t_X);
        mpz_clear(t_Y_sq);
    }

    res.iterations = total_iter;
    res.execution_time = omp_get_wtime() - start_time;
    res.success = found;

    mpz_clear(N);
    mpz_clear(sqrtN);

    return res;
}
