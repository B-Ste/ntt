#define N 4096
#define BITLENGTH 12    // 2 ^ 12 = 4096
#define NUM_Q 13
#define L 30            // bit-length of moduli Q
#define C 32            // number of cores

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
void ntt_packed(int32_t, int32_t*, int64_t*);
void ntt_multi_packed(int32_t, int32_t*);
void intt(int32_t, int32_t, int32_t*, int32_t*);
void pack(int32_t*, int64_t*);
void unpack(int64_t*, int32_t*);
void unpack_after_ntt(int64_t*, int32_t*);
void load_memory(int64_t*);
void unload_memory(int64_t*);
void unload_after_ntt(int64_t*);
int64_t pack_single(int32_t, int32_t);
int32_t get_first(int64_t);
int32_t get_second(int64_t);
int32_t* brv_powers(int32_t, int32_t, int32_t);
int32_t brv(int32_t, int32_t);
int32_t mod(int64_t, int32_t);
int32_t mod_barrett(int64_t, int32_t); 
int32_t mod_addition(int32_t, int32_t);
int32_t mod_subtraction(int32_t, int32_t);
int performance_test();
int correctness_test_normal_ntt();
int packing_test();
int correctness_test_packed_ntt();
int correctness_test_multi_packed_ntt();
int memory_test();

int main(int argc, char const *argv[]) {
    if (argc != 1) {
        printf("This program is run without arguments.\n");
        return 1;
    }

    int program_selection;
    printf("Select program:\n");
    printf("1) performance-test\n");
    printf("2) correctness-test for normal NTT\n");
    printf("3) correctness-test for packing\n");
    printf("4) correctness-test for packed NTT\n");
    printf("5) correctness-test for multi-core packed NTT\n");
    printf("6) correctness-test for memory load & unload\n");
    printf("Selection: ");
    scanf("%d", &program_selection);

    switch (program_selection)
    {
    case 1: return performance_test();
    case 2: return correctness_test_normal_ntt();
    case 3: return packing_test();
    case 4: return correctness_test_packed_ntt();
    case 5: return correctness_test_multi_packed_ntt();
    case 6: return memory_test();
    default:
        printf("Unknown program. Exiting.\n");
        return 1;
    }
}

int performance_test() {
    int timing_runs;
    printf("Number of timing-runs: ");
    scanf("%d", &timing_runs);

    unsigned int seed2 = clock();
    double total_time = 0;
    int k = 0;
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

int correctness_test_normal_ntt() {
    int correctness_runs;
    printf("Number of correctness-runs: ");
    scanf("%d", &correctness_runs);

    unsigned int seed1 = clock();
    int k = 0;
    int fail = 0;
    for (int j = 0; j < correctness_runs; j++) {
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
                fail = 1;
                break;
            }
        }
        free(a); free(b); free(c1); free(c2); free(psis); free(psis_ns);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    if (fail) {
        printf("Failed tests on correctness test for normal NTT.\n");
        return 1;
    } else {
        printf("Passed all correctness runs on correctness test for normal NTT.\n");
        return 0;
    }
}

int packing_test() {
    int correctness_runs;
    printf("Number of correctness-runs: ");
    scanf("%d", &correctness_runs);

    unsigned int seed1 = clock();
    int k = 0;
    int fail = 0;
    for (int j = 0; j < correctness_runs; j++) {
        int32_t* a = calloc(N, sizeof(int32_t));
        int32_t* b = calloc(N, sizeof(int32_t));
        int64_t* g = calloc(N/2, sizeof(int64_t));
        int32_t q = modulus[k];
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            a[i] = mod_barrett(rand_r(&seed1), q);
        }
        pack(a, g);
        unpack(g, b);
        for (int i = 0; i < N; i++) {
            if (a[i] != b[i]) {
                printf("Failed correctness run %i with k = %i at i = %i with a[i] = %x b[i] = %x\n", j, k, i, a[i], b[i]);
                fail = 1;
                break;
            }
        }
        free(a); free(b); free(g);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    if (fail) {
        printf("Failed correctness runs on correctness test for packing.\n");
        return 1;
    } else {
        printf("Passed all correctness runs on correctness test for packing.\n");
        return 0;
    }
}

int correctness_test_packed_ntt() {
    int correctness_runs;
    printf("Number of correctness-runs: ");
    scanf("%d", &correctness_runs);

    unsigned int seed1 = clock();
    int k = 0;
    int fail = 0;
    for (int j = 0; j < correctness_runs; j++) {
        int32_t* a = calloc(N, sizeof(int32_t));
        int32_t* b = calloc(N, sizeof(int32_t));
        int32_t q = modulus[k];
        int32_t psi_p = psi[k];
        int32_t* psis = brv_powers(psi_p, q, N);
        int64_t* g = calloc(N/2, sizeof(int64_t));
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            a[i] = mod_barrett(rand_r(&seed1), q);
        }
        pack(a, g);
        ntt(q, psis, a);
        ntt_packed(q, psis, g);
        unpack_after_ntt(g, b);
        for (int i = 0; i < N; i++) {
            if (a[i] != b[i]) {
                printf("Failed correctness run %i with k = %i at i = %i\n", j, k, i);
                fail = 0;
                break;
            }
        }
        free(a); free(b); free(g); free(psis);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    if (fail) {
        printf("Failed correctness runs on correctness test for packed NTT.\n");
        return 1;
    } else {
        printf("Passed all correctness runs on correctness test for packed NTT.\n");
        return 0;
    }
}

int memory_test() {
    int correctness_runs;
    printf("Number of correctness-runs: ");
    scanf("%d", &correctness_runs);

    unsigned int seed1 = clock();
    int k = 0;
    int fail = 0;
    for (int j = 0; j < correctness_runs; j++) {
        int32_t* a = calloc(N, sizeof(int32_t));
        int32_t* b = calloc(N, sizeof(int32_t));
        int64_t* g = calloc(N/2, sizeof(int64_t));
        int32_t q = modulus[k];
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            a[i] = mod_barrett(rand_r(&seed1), q);
        }
        pack(a, g);
        load_memory(g);
        unload_memory(g);
        unpack(g, b);
        for (int i = 0; i < N; i++) {
            if (a[i] != b[i]) {
                printf("Failed correctness run %i with k = %i at i = %i with a[i] = %x b[i] = %x\n", j, k, i, a[i], b[i]);
                fail = 1;
                break;
            }
        }
        free(a); free(b); free(g);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    if (fail) {
        printf("Failed correctness runs on correctness test for memory.\n");
        return 1;
    } else {
        printf("Passed all correctness runs on correctness test for memory.\n");
        return 0;
    }
}

int correctness_test_multi_packed_ntt() {
    int correctness_runs;
    printf("Number of correctness-runs: ");
    scanf("%d", &correctness_runs);

    unsigned int seed1 = clock();
    int k = 0;
    int fail = 0;
    for (int j = 0; j < correctness_runs; j++) {
        int32_t* a = calloc(N, sizeof(int32_t));
        int32_t* b = calloc(N, sizeof(int32_t));
        int32_t q = modulus[k];
        int32_t psi_p = psi[k];
        int32_t* psis = brv_powers(psi_p, q, N);
        int64_t* g = calloc(N/2, sizeof(int64_t));
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            a[i] = mod_barrett(rand_r(&seed1), q);
        }
        pack(a, g);
        ntt(q, psis, a);
        load_memory(g);
        ntt_multi_packed(q, psis);
        unload_after_ntt(g);
        unpack_after_ntt(g, b);
        for (int i = 0; i < N; i++) {
            if (a[i] != b[i]) {
                printf("Failed correctness run %i with k = %i at i = %i\n", j, k, i);
                fail = 1;
                break;
            }
        }
        free(a); free(b); free(g); free(psis);
        if (k == NUM_Q - 1) k = 0;
        else k++;
    }
    if (fail) {
        printf("Failed correctness runs on correctness test for multi-core packed NTT.\n");
        return 1;
    } else {
        printf("Passed all correctness runs on correctness test for multi-core packed NTT.\n");
        return 0;
    }
}



///////////////////
// NAIVE VERSION //
///////////////////

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



//////////////////
// STANDARD NTT //
//////////////////

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
 * Calculates the forward NTT of a polynomial a with modulus q and bit-reverse ordered list psis of powers of
 * a 2n-th root of unity. The result is saved in a in bit-reversed order.
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



////////////////
// PACKED NTT //
////////////////

/**
 * Calculates the forward NTT of a packed poylnomial a with modulus q and bit-reverse ordered list psis of powers of
 * a 2n-th root of unity. The result is saved in a in packed bit-reversed order.
 */
void ntt_packed(int32_t q, int32_t* psis, int64_t* a) {
    int t = N/2;
    for (int m = 1; m < N/2; m *= 2) {
        t /= 2;
        for (int i = 0; i < m; i++) {
            int j1 = 2 * i * t;
            int j2 = j1 + t - 1;
            int32_t w = psis[m + i];
            for (int j = j1; j <= j2; j++) {
                int32_t u1 = get_first(a[j]);
                int32_t v1 = MOD_MUL((int64_t) w * get_second(a[j]), q);
                int32_t u2 = get_first(a[j + t]);
                int32_t v2 = MOD_MUL((int64_t) w * get_second(a[j + t]), q);
                int32_t r1 = MOD_ADD(u1 + v1, q);
                int32_t r2 = MOD_SUB(u1 - v1, q);
                int32_t r3 = MOD_ADD(u2 + v2, q);
                int32_t r4 = MOD_SUB(u2 - v2, q);
                a[j] = pack_single(r1, r3);
                a[j + t] = pack_single(r2, r4);
            }
        }
    }
    for (int j = 0; j < N/2; j++) {
        int32_t w = psis[N/2 + j];
        int32_t U = get_first(a[j]);
        int32_t V = MOD_MUL((int64_t) w * get_second(a[j]), q);
        int32_t r1 = MOD_ADD(U + V, q);
        int32_t r2 = MOD_SUB(U - V, q);
        a[j] = pack_single(r1, r2);
    }
}



/////////////////////////////////////
// MULTICORE PACKED NTT SIMULATION //
/////////////////////////////////////

int64_t memory[C][2][N / (4 * C)];
int32_t router_input[C][N / (4 * C)][4];

/**
 * Loads packed polynomial a into memory.
 */
void load_memory(int64_t* a) {
    for (int k = 0; k < C; k++) {
        for (int j = 0; j < N / (4 * C); j++) {
            memory[k][0][j] = a[k * (N / (4 * C)) + j];
            memory[k][1][j] = a[k * (N / (4 * C)) + j + (N / 4)];
        }
    }
}

/**
 * Unloads packed polynomial from memory.
 */
void unload_memory(int64_t* a) {
    for (int k = 0; k < C; k++) {
        for (int j = 0; j < N / (4 * C); j++) {
            a[k * (N / (4 * C)) + j] = memory[k][0][j];
            a[k * (N / (4 * C)) + j + (N / 4)] = memory[k][1][j];
        }
    }
}

/**
 * Unloads packed polynomial from memory after ntt.
 */
void unload_after_ntt(int64_t* a) {
    for (int k = 0; k < C; k++) {
        for (int j = 0; j < N / (4 * C); j++) {
            a[k * (N / (2 * C)) + 2 * j] = memory[k][0][j];
            a[k * (N / (2 * C)) + 2 * j + 1] = memory[k][1][j];
        }
    }
}

/**
 * Simulates multi-core packed NTT for contents in memory variable.
 */
void ntt_multi_packed(int32_t q, int32_t* psis) {
    int32_t t = N / 2;
    for (int m = 1; m < N/2; m *= 2) {
        t /= 2;

        /*  We only simulate the multi-core aspect. Therefor, the following for-loop iteratively 
            calculates the results for each core and saves them in the router_input. The router will write
            the results into their corresponding memory-location. 
        */
        for (int k = 0; k < C; k++) {

            if (t > N / (4 * C)) {
                int32_t w = psis[m + (k * C) / t];
                for (int j = 0; j < N / (4 * C); j++) {
                    int32_t u1 = get_first(memory[k][0][j]);
                    int32_t v1 = MOD_MUL((int64_t) w * get_second(memory[k][0][j]), q);
                    int32_t u2 = get_first(memory[k][1][j]);
                    int32_t v2 = MOD_MUL((int64_t) w * get_second(memory[k][1][j]), q);
                    int32_t r1 = MOD_ADD(u1 + v1, q);
                    int32_t r2 = MOD_SUB(u1 - v1, q);
                    int32_t r3 = MOD_ADD(u2 + v2, q);
                    int32_t r4 = MOD_SUB(u2 - v2, q);
                    router_input[k][j][0] = r1;
                    router_input[k][j][1] = r2;
                    router_input[k][j][2] = r3;
                    router_input[k][j][3] = r4;
                }
            } else {
                for (int i = 0; i < m / C; i++) {
                    int32_t w = psis[m + (k * m) / C + i];
                    int j1 = i * t;
                    int j2 = j1 + t - 1;
                    for (int j = j1; j <= j2; j++) {
                        int32_t u1 = get_first(memory[k][0][j]);
                        int32_t v1 = MOD_MUL((int64_t) w * get_second(memory[k][0][j]), q);
                        int32_t u2 = get_first(memory[k][1][j]);
                        int32_t v2 = MOD_MUL((int64_t) w * get_second(memory[k][1][j]), q);
                        int32_t r1 = MOD_ADD(u1 + v1, q);
                        int32_t r2 = MOD_SUB(u1 - v1, q);
                        int32_t r3 = MOD_ADD(u2 + v2, q);
                        int32_t r4 = MOD_SUB(u2 - v2, q);
                        router_input[k][j][0] = r1;
                        router_input[k][j][1] = r2;
                        router_input[k][j][2] = r3;
                        router_input[k][j][3] = r4;
                    }
                }
            }

        }

        // Routing algorithm. Will be its own module on the FPGA.
        for (int k = 0; k < C; k++) {
            for (int j = 0; j < N / (4 * C); j++) {
                if (t > N / (4 * C)) {
                    if (((k * 2 * m) / C) % 2 == 0) {
                        memory[k][0][j] = pack_single(router_input[k][j][0], router_input[k][j][2]);
                        memory[k][1][j] = pack_single(router_input[k + ((2 * t * C) / N)][j][0], router_input[k + ((2 * t * C) / N)][j][2]);
                    } else {
                        memory[k][0][j] = pack_single(router_input[k - ((2 * t * C) / N)][j][1], router_input[k - ((2 * t * C) / N)][j][3]);
                        memory[k][1][j] = pack_single(router_input[k][j][1], router_input[k][j][3]);    
                    }
                } else {
                    if (t != 1) {
                        if (j % t < t/2) {
                            memory[k][0][j] = pack_single(router_input[k][j][0], router_input[k][j][2]);
                            memory[k][0][j + t/2] = pack_single(router_input[k][j][1], router_input[k][j][3]);
                        } else {
                            memory[k][1][j - t/2] = pack_single(router_input[k][j][0], router_input[k][j][2]);
                            memory[k][1][j] = pack_single(router_input[k][j][1], router_input[k][j][3]);
                        }
                    } else {
                        memory[k][0][j] = pack_single(router_input[k][j][0], router_input[k][j][2]);
                        memory[k][1][j] = pack_single(router_input[k][j][1], router_input[k][j][3]);
                    }
                }
            }
        }
    }

    // Again for simulation purposes
    // We write diretly to memory, as no complex reordering is required
    for (int k = 0; k < C; k++) {

        for (int j = 0; j < N / (4 * C); j++) {
            int32_t w1 = psis[N/2 + k * N / (2 * C) + 2 * j];
            int32_t w2 = psis[N/2 + k * N / (2 * C) + 2 * j + 1];
            int32_t u1 = get_first(memory[k][0][j]);
            int32_t v1 = MOD_MUL((int64_t) w1 * get_second(memory[k][0][j]), q);
            int32_t u2 = get_first(memory[k][1][j]);
            int32_t v2 = MOD_MUL((int64_t) w2 * get_second(memory[k][1][j]), q);
            int32_t r1 = MOD_ADD(u1 + v1, q);
            int32_t r2 = MOD_SUB(u1 - v1, q);
            int32_t r3 = MOD_ADD(u2 + v2, q);
            int32_t r4 = MOD_SUB(u2 - v2, q);
            memory[k][0][j] = pack_single(r1, r2);
            memory[k][1][j] = pack_single(r3, r4);
        }

    }
}


//////////////////////////////
// PACKING HELPER FUNCTIONS //
//////////////////////////////

/**
 * Packs two 30-bit coefficients into one 60-bit word, with a1 in the lower significance and
 * a2 in the higher significance bits.
 */
int64_t pack_single(int32_t a1, int32_t a2) {
    return a1 | ((int64_t) a2 << 30);
}

/**
 * Returns first 30-bit coefficient in 60-bit word containing two coefficients.
 */
int32_t get_first(int64_t g) {
    return (int32_t) (g & 1073741823);
}

/**
 * Returns second 30-bit coefficient in 60-bit word containing two coefficients.
 */
int32_t get_second(int64_t g) {
    return (int32_t) ((g >> 30) & 1073741823);
}

/**
 * Packs polynomial a of N 30-bit coefficients into a N/2 long array of 60-bit words, where
 * word i contains coefficients a[i] and a[i + N/2].
 */
void pack(int32_t* a, int64_t* g) {
    for (int i = 0; i < N/2; i++) g[i] = pack_single(a[i], a[i + N/2]);
}

/**
 * Unpacks 30-bit coefficients from N/2 long 60-bit word array into length N polynomial, where
 * word i contains coefficients a[i] and a[i + N/2].
 */
void unpack(int64_t* g, int32_t* a) {
    for (int i = 0; i < N/2; i++) {
        a[i] = get_first(g[i]);
        a[i + N/2] = get_second(g[i]);
    }
}

/**
 * Unpacks 30-bit coefficients from N/2 long 60-bit word array into length N polynomial, where
 * word i contains coefficients a[2 * i] and a[2 * i + ].
 */
void unpack_after_ntt(int64_t* g, int32_t* a) {
    for (int i = 0; i < N/2; i++) {
        a[2 * i] = get_first(g[i]);
        a[2 * i + 1] = get_second(g[i]);
    }
}



//////////////////////
// HELPER FUNCTIONS //
//////////////////////

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



//////////////////////
// MODULO FUNCTIONS //
//////////////////////

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
