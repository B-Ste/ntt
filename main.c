#define N 4096
#define BITLENGTH 12 // 2 ^ 12 = 4096
#define NUM_Q 6

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

const int64_t modulus[NUM_Q] = {1068564481, 1069219841, 1070727169, 1071513601, 1072496641, 1073479681};
const int64_t psi[NUM_Q] = {485176247, 429380462, 1026271992, 250036504, 378716967, 371836615};
const int64_t psi_neg[NUM_Q] = {320928324, 491001446, 874272385, 284619963, 371608897, 730652244};
const int64_t n_neg[NUM_Q] = {1068303601, 1068958801, 1070465761, 1071252001, 1072234801, 1073217601};

void nwc_naive(int64_t, int32_t*, int32_t*, int32_t*);
void nwc_ntt(int64_t, int64_t*, int64_t*, int64_t, int32_t*, int32_t*, int32_t*);
void ntt(int64_t, int64_t*, int32_t*);
void intt(int64_t, int64_t, int64_t*, int32_t*);
int64_t* brv_powers(int64_t, int64_t, int64_t);
int64_t brv(int64_t, int64_t);

int32_t mod(int64_t a, int64_t q) {
    int64_t r = a % q;
    return (int32_t) r < 0 ? r + q : r;
}

int main(int argc, char const *argv[]) {
    int k = 0;
    int32_t* a = malloc(N * sizeof(int32_t));
    int32_t* b = malloc(N * sizeof(int32_t));
    int32_t* c1 = malloc(N * sizeof(int32_t));
    int32_t* c2 = malloc(N * sizeof(int32_t));
    int64_t q = modulus[k];
    int64_t psi_p = psi[k];
    int64_t psi_n = psi_neg[k];
    int64_t n_n = n_neg[k];
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        a[i] = mod(rand(), q);
        b[i] = mod(rand(), q);
    }
    nwc_naive(q, a, b, c1);
    int64_t* psis = brv_powers(psi_p, q, N);
    int64_t* psis_ns = brv_powers(psi_n, q, N);
    nwc_ntt(q, psis, psis_ns, n_n, a, b, c2);
    for (int i = 0; i < N; i++) {
        if (c1[i] != c2[i]) {
            printf("not equal, i = %i\n", i);
            break;
        }
    }
    free(psis);
    free(psis_ns);
    return 0;
}

/**
 * Calculates the NWC of a and b with modulus q naively in O(n^2) and saves it in c.
 */
void nwc_naive(int64_t q, int32_t* a, int32_t* b, int32_t* c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int64_t c_k = mod((u_int64_t) a[i] * b[j], q);
            if (i + j < N) c[i + j] = mod(c_k + c[i + j], q);
            else c[i + j - N] = mod(c[i + j - N] - c_k, q);
        }
    }
}

/**
 * Calculates the NWV of a and b with modulus q, tables psi and psi^-1 and multiplicative inverse of n
 * and saves it in c.
 */
void nwc_ntt(int64_t q, int64_t* psis, int64_t* psi_ns, int64_t n_neg, int32_t* a, int32_t* b, int32_t* c) {
    ntt(q, psis, a);
    ntt(q, psis, b);
    for (int i = 0; i < N; i++) c[i] = mod((int64_t) a[i] * b[i], q);
    intt(n_neg, q, psi_ns, c);
}

/**
 * Calculates the forward NTT of a with modulus q and bit-reverse ordered list psis of powers of
 * a 2n-th root of unity. The result is saved in a int bit-reversed order.
 */
void ntt(int64_t q, int64_t* psis, int32_t* a) {
    int64_t t = N;
    for (int m = 1; m < N; m *= 2) {
        t /= 2;
        for (int i = 0; i < m; i++) {
            int64_t j_1 = 2 * i * t;
            int64_t j_2 = j_1 + t - 1;
            int64_t w = psis[m + i];
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
void intt(int64_t inv_n, int64_t q, int64_t* psis, int32_t* a) {
    int64_t t = 1;
    for (int m = N; m > 1; m /= 2) {
        int64_t j_1 = 0;
        int64_t h = m / 2;
        for (int i = 0; i < h; i++) {
            int64_t j_2 = j_1 + t - 1;
            int64_t w = psis[h + i];
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
    for (int j = 0; j < N; j++) {
        int64_t p = a[j] * inv_n;
        a[j] = mod(p, q);
    }
}

/**
 * bit-reverses i, which is of length n
 */
int64_t brv(int64_t i, int64_t n) {
    int64_t brv = 0;
    while (n--) {
        brv = (brv << 1) | (i & 1);
        i >>= 1;
    }
    return brv;
}

/**
 * Creates a list with length n of 1,p,p^2, ..., p^(n-1) mod q in bit-reversed order.
 */
int64_t* brv_powers(int64_t p, int64_t q, int64_t n) {
    int64_t* brv_pwr = malloc(n * (sizeof(int64_t)));
    int64_t tmp = 1;
    brv_pwr[0] = 1;
    for (int i = 1; i < n; i++) {
        int64_t r = mod(tmp * p, q);
        brv_pwr[brv(i, BITLENGTH)] = r;
        tmp = r;
    }
    return brv_pwr;
}
