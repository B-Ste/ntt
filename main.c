#define N 4096
#define BITLENGTH 12    // 2 ^ 12 = 4096
#define NUM_Q 13
#define L 30            // bit-length of moduli Q

#ifdef BAR
    #define MOD_ADD(a, q) mod_addition(a, q)
    #define MOD_SUB(a, q) mod_subtraction(a, q)
    #define MOD_MUL(a, q) mod_barrett(a, q)
#else
    #define MOD_ADD(a, q) mod(a, q)
    #define MOD_SUB(a, q) mod(a, q)
    #define MOD_MUL(a, q) mod(a, q)
#endif

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

// extra primes: 1063321601, 1063452673, 1064697857, 1065484289, 1065811969, 1068236801, 1068433409

const int32_t modulus[NUM_Q] = {1063321601, 1063452673, 1064697857, 1065484289, 1065811969, 1068236801, 
    1068433409, 1068564481, 1069219841, 1070727169, 1071513601, 1072496641, 1073479681};
const int32_t psi[NUM_Q] = {210627128, 528517222, 152301028, 433219053, 644958739, 932161456, 
    839237106, 485176247, 429380462, 1026271992, 250036504, 378716967, 371836615};
const int32_t psi_neg[NUM_Q] = {806968250, 496328859, 1011460660, 468376471, 167589381, 773385751, 
    280720051, 320928324, 491001446, 874272385, 284619963, 371608897, 730652244};
const int32_t n_neg[NUM_Q] = {1063062001, 1063193041, 1064437921, 1065224161, 1065551761, 1067976001, 
    1068172561, 1068303601, 1068958801, 1070465761, 1071252001, 1072234801, 1073217601};

void nwc_naive(int32_t, int32_t*, int32_t*, int32_t*);
void nwc_ntt(int32_t, int32_t*, int32_t*, int32_t, int32_t*, int32_t*, int32_t*);
void ntt(int32_t, int32_t*, int32_t*);
void intt(int32_t, int32_t, int32_t*, int32_t*);
int32_t* brv_powers(int32_t, int32_t, int32_t);
int32_t brv(int32_t, int32_t);
int32_t mod(int64_t, int32_t);
int32_t mod_barrett(int64_t, int32_t); 
int32_t mod_addition(int32_t, int32_t);
int32_t mod_subtraction(int32_t, int32_t);

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        printf("Synopsis: ntt_bench <correctness runs> <timing runs>\n");
        return 1;
    }

    unsigned int seed1 = clock();
    int k = 0;
    for (int j = 0; j < strtol(argv[1], NULL, 10); j++) {
        int32_t* a = calloc(N, sizeof(int32_t));
        int32_t* b = calloc(N, sizeof(int32_t));
        int32_t* c1 = calloc(N, sizeof(int32_t));
        int32_t* c2 = calloc(N, sizeof(int32_t));
        int32_t q = modulus[k];
        int32_t psi_p = psi[k];
        int32_t psi_n = psi_neg[k];
        int32_t n_n = n_neg[k];
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            a[i] = mod_barrett(rand_r(&seed1), q);
            b[i] = mod_barrett(rand_r(&seed1), q);
        }
        nwc_naive(q, a, b, c1);
        int32_t* psis = brv_powers(psi_p, q, N);
        int32_t* psis_ns = brv_powers(psi_n, q, N);
        nwc_ntt(q, psis, psis_ns, n_n, a, b, c2);
        for (int i = 0; i < N; i++) {
            if (c1[i] != c2[i]) {
                printf("Failed correctness run %i with k = %i at i = %i\n", j, k, i);
                break;
            }
        }
        free(a); free(b); free(c1); free(c2); free(psis); free(psis_ns);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    printf("Passed all correctness runs\n");

    unsigned int seed2 = clock();
    int timing_runs = strtol(argv[2], NULL, 10);
    double total_time = 0;
    k = 0;
    clock_t t;
    for (int j = 0; j < timing_runs; j++) {
        int32_t* a = calloc(N, sizeof(int32_t));
        int32_t* b = calloc(N, sizeof(int32_t));
        int32_t* c = calloc(N, sizeof(int32_t));
        int32_t q = modulus[k];
        int32_t psi_p = psi[k];
        int32_t psi_n = psi_neg[k];
        int32_t n_n = n_neg[k];
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            a[i] = mod_barrett(rand_r(&seed2), q);
            b[i] = mod_barrett(rand_r(&seed2), q);
        }
        int32_t* psis = brv_powers(psi_p, q, N);
        int32_t* psis_ns = brv_powers(psi_n, q, N);

        t = clock();
        nwc_ntt(q, psis, psis_ns, n_n, a, b, c);
        t = clock() - t;
        total_time += (double) t / CLOCKS_PER_SEC;

        free(a); free(b); free(c); free(psis); free(psis_ns);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    double average_time = total_time / timing_runs;
    printf("Average time in %i timing-runs: %fs\n", timing_runs, average_time);
    return 0;
}

/**
 * Calculates the NWC of a and b with modulus q naively in O(n^2) and saves it in c.
 */
void nwc_naive(int32_t q, int32_t* a, int32_t* b, int32_t* c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int32_t c_k = MOD_MUL((int64_t) a[i] * b[j], q);
            if (i + j < N) c[i + j] = MOD_ADD(c_k + c[i + j], q);
            else c[i + j - N] = MOD_SUB(c[i + j - N] - c_k, q);
        }
    }
}

/**
 * Calculates the NWC of a and b with modulus q, tables psi and psi^-1 and multiplicative inverse of n
 * and saves it in c.
 */
void nwc_ntt(int32_t q, int32_t* psis, int32_t* psi_ns, int32_t n_neg, int32_t* a, int32_t* b, int32_t* c) {
    ntt(q, psis, a);
    ntt(q, psis, b);
    for (int i = 0; i < N; i++) c[i] = MOD_MUL((int64_t) a[i] * b[i], q);
    intt(n_neg, q, psi_ns, c);
}

/**
 * Calculates the forward NTT of a with modulus q and bit-reverse ordered list psis of powers of
 * a 2n-th root of unity. The result is saved in a int bit-reversed order.
 */
void ntt(int32_t q, int32_t* psis, int32_t* a) {
    int t = N;
    for (int m = 1; m < N; m *= 2) {
        t /= 2;
        for (int i = 0; i < m; i++) {
            int j_1 = 2 * i * t;
            int j_2 = j_1 + t - 1;
            int32_t w = psis[m + i];
            for (int j = j_1; j <= j_2; j++) {
                int32_t U = a[j];
                int32_t V = MOD_MUL((int64_t) w * a[j + t], q);
                a[j] = MOD_ADD(U + V, q);
                a[j + t] = MOD_SUB(U - V, q);
            }
        }
    }
}

/**
 * Calculates the inverse NTT (INTT) of bit-reversed a with modulus q, multiplicative inverse of N, bit-reverse 
 * ordered list psis of powers of the multiplicative inverse of a 2n-th root of unity. The result is saved in a 
 * in normal order.
 */
void intt(int32_t inv_n, int32_t q, int32_t* psis, int32_t* a) {
    int t = 1;
    for (int m = N; m > 1; m /= 2) {
        int j_1 = 0;
        int h = m / 2;
        for (int i = 0; i < h; i++) {
            int j_2 = j_1 + t - 1;
            int32_t w = psis[h + i];
            for (int j = j_1; j <= j_2; j++) {
                int32_t U = a[j];
                int32_t V = a[j + t];
                int32_t inner = MOD_SUB(U - V, q);
                a[j] = MOD_ADD(U + V, q);
                a[j + t] = MOD_MUL((int64_t) w * inner, q); 
            }
            j_1 += 2 * t;
        }
        t *= 2;
    }
    for (int j = 0; j < N; j++) a[j] = MOD_MUL((int64_t) a[j] * inv_n, q);
}

/**
 * Bit-reverses i, which is of length n
 */
int32_t brv(int32_t i, int32_t n) {
    int32_t brv = 0;
    while (n--) {
        brv = (brv << 1) | (i & 1);
        i >>= 1;
    }
    return brv;
}

/**
 * Creates a list with length n of 1,p,p^2, ..., p^(n-1) mod q in bit-reversed order.
 */
int32_t* brv_powers(int32_t p, int32_t q, int32_t n) {
    int32_t* brv_pwr = malloc(n * (sizeof(int32_t)));
    int32_t tmp = 1;
    brv_pwr[0] = 1;
    for (int i = 1; i < n; i++) {
        int32_t r = mod_barrett((int64_t) tmp * p, q);
        brv_pwr[brv(i, BITLENGTH)] = r;
        tmp = r;
    }
    return brv_pwr;
}

/**
 * Calculates a mod q naively.
*/
int32_t mod(int64_t a, int32_t q) {
    int32_t r = a % q;
    return r < 0 ? r + q : r;
}

/**
 * Calculates a mod q, where a is in [0, 2q - 1].
*/
int32_t mod_addition(int32_t a, int32_t q) {
    return a > q ? a - q : a;
}

/**
 * Calculates a mod q, where a is in [-q + 1, q - 1].
*/
int32_t mod_subtraction(int32_t a, int32_t q) {
    return a < 0 ? a + q : a;
}

/**
 * Calculates a mod q by Barrett-Reduction, where a is in [0, q^2 - 2q + 1].
 * Algorithm by: https://doi.org/10.1109/ICSICT.2012.6467863
*/
int32_t mod_barrett(int64_t a, int32_t q) {
    const uint64_t lpower_2Lp3 = 1LL << (2 * L + 3);
    const uint64_t lpower_Lm2 = 1LL << (L - 2);
    const uint64_t lpower_Lp5 = 1LL << (L +  5);
    const uint64_t lpower_Lp1 = 1LL << (L + 1);
    const uint64_t lambda = lpower_2Lp3 / q;
    uint64_t q_hat = floor(floor(a / lpower_Lm2) * lambda) / lpower_Lp5;
    int64_t r0 = (a - q_hat * q) & (lpower_Lp1 - 1);
    int64_t r1 = r0 - q;
    return r1 >= 0 ? r1 : r0;
}
