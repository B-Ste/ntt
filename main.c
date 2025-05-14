#define N 4096
#define BITLENGTH 12 // 2 ^ 12 = 4096
#define NUM_Q 6

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

const int32_t modulus[NUM_Q] = {1068564481, 1069219841, 1070727169, 1071513601, 1072496641, 1073479681};
const int32_t psi[NUM_Q] = {485176247, 429380462, 1026271992, 250036504, 378716967, 371836615};
const int32_t psi_neg[NUM_Q] = {320928324, 491001446, 874272385, 284619963, 371608897, 730652244};
const int32_t n_neg[NUM_Q] = {1068303601, 1068958801, 1070465761, 1071252001, 1072234801, 1073217601};

void nwc_naive(int32_t, int32_t*, int32_t*, int32_t*);
void nwc_ntt(int32_t, int32_t*, int32_t*, int32_t, int32_t*, int32_t*, int32_t*);
void ntt(int32_t, int32_t*, int32_t*);
void intt(int32_t, int32_t, int32_t*, int32_t*);
int32_t* brv_powers(int32_t, int32_t, int32_t);
int32_t brv(int32_t, int32_t);
int32_t mod(int64_t, int32_t);

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
            a[i] = mod(rand_r(&seed1), q);
            b[i] = mod(rand_r(&seed1), q);
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
            a[i] = mod(rand_r(&seed2), q);
            b[i] = mod(rand_r(&seed2), q);
        }
        int32_t* psis = brv_powers(psi_p, q, N);
        int32_t* psis_ns = brv_powers(psi_n, q, N);

        t = clock();
        nwc_ntt(q, psis, psis_ns, n_n, a, b, c);
        t = clock() - t;
        double time = (double) t / CLOCKS_PER_SEC;
        total_time += time;

        free(a); free(b); free(c); free(psis); free(psis_ns);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    double average = total_time / timing_runs;
    printf("Average time taken for %i timing-runs: %fs\n", timing_runs, average);
    return 0;
}

/**
 * Calculates the NWC of a and b with modulus q naively in O(n^2) and saves it in c.
 */
void nwc_naive(int32_t q, int32_t* a, int32_t* b, int32_t* c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int64_t c_k = mod((u_int64_t) a[i] * b[j], q);
            if (i + j < N) c[i + j] = mod(c_k + c[i + j], q);
            else c[i + j - N] = mod(c[i + j - N] - c_k, q);
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
    for (int i = 0; i < N; i++) c[i] = mod((int64_t) a[i] * b[i], q);
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
                int64_t U = a[j];
                int64_t V = (int64_t) w * a[j + t];
                a[j] = mod(U + V, q);
                a[j + t] = mod(U - V, q);
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
                int64_t U = a[j];
                int64_t V = a[j + t];
                a[j] = mod(U + V, q);
                a[j + t] = mod(w * (U - V), q); 
            }
            j_1 += 2 * t;
        }
        t *= 2;
    }
    for (int j = 0; j < N; j++) a[j] = mod((int64_t) a[j] * inv_n, q);
}

/**
 * bit-reverses i, which is of length n
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
        int32_t r = mod((int64_t) tmp * p, q);
        brv_pwr[brv(i, BITLENGTH)] = r;
        tmp = r;
    }
    return brv_pwr;
}

/**
 * Computes a mod q, with the result beeing in [0,q-1].
 */
int32_t mod(int64_t a, int32_t q) {
    int64_t r = a % q;
    return (int32_t) r < 0 ? r + q : r;
}
