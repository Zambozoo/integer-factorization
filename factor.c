//#define THREADED

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#ifdef THREADED
#include <pthread.h>
#endif
#include <gmp.h>

// There exists a set of numbers S[n]
//   such that for any distinct primes a,b <= sqrt(n),
//   there exists a number m in S[n] such that a | m /\ ~(b | m)

// Let |S[n]| = l,
//   l choose l / 2 >= pi(sqrt(n)), where pi(x) = number of primes <= x
//   Therefore, min(|S[n]|) = k where k choose k/2 >= pi(n) /\ (k-2) choose (k-2)/2 < pi(n)


// There exists S[n] for all n {Is 2 | pi(sqrt(n)) necessary?}
//   such that for every m in S[n], 
//   the prime factors of m all have multiplicity 1
//   and #_prime_factors(m) = pi(sqrt(n)) / 2

// Let T[n] be the set of all S[n] with minimal cardinality
//   |T[n]| = |S[n]| choose pi(n)

#pragma region Random
uint64_t rand_uint64() { //Generate random number
    uint64_t r = 0;
    for (int i = 0; i < 64; i++) {
        r = (r << 1) + (rand() & 1);
    }
    return r;
}
uint32_t rand_uint32() { //Generate random number
    uint32_t r = 0;
    for (int i = 0; i < 32; i++) {
        r = (r << 1) + (rand() & 1);
    }
    return r;
}

long long mulmod(long long a, long long b, long long mod) { //Modular multiplication
    long long x = 0, y = a % mod;
    while (b > 0) {
        if (b % 2 == 1)
            x = (x + y) % mod;
        y = (y * 2) % mod;
        b /= 2;
    }
    return x % mod;
}

long long modulo(long long base, long long exponent, long long mod) { //Modular exponentiation
    long long x = 1;
    long long y = base;
    while (exponent > 0)
    {
        if (exponent % 2 == 1)
            x = (x * y) % mod;
        y = (y * y) % mod;
        exponent = exponent / 2;
    }
    return x % mod;
}

int Miller(long long p,int iteration) { //Miller-Rabin Primality test
    int i;
    long long s;
    if (p < 2 || (p != 2 && p % 2 == 0))
        return 0;

    s = p - 1;
    while (s % 2 == 0)
        s /= 2;

    for (i = 0; i < iteration; i++) {
        long long a = rand() % (p - 1) + 1, temp = s;
        long long mod = modulo(a, temp, p);
        while (temp != p - 1 && mod != 1 && mod != p - 1) {
            mod = mulmod(mod, mod, p);
            temp *= 2;
        }
        if (mod != p - 1 && temp % 2 == 0)
            return 0;
    }
    return 1;
}

uint32_t rand_uint32_prime() { //Return random miller(10) assured uint32 prime
    while(true) {
        int i = rand_uint32();
        if(Miller(i, 10)) return i;
    }
}
#pragma endregion

#pragma region BigInt
typedef mpz_t BigInt;

double mpz_log2(BigInt x) { //Calculate log2(BigInt)
    signed long int ex;
    const double di = mpz_get_d_2exp(&ex, x);
    return log(di) + log(2) * (double) ex;
}

static void askBigInt(const char* msg, BigInt n) { //Request Integer from std in
    mpz_init(n);
    mpz_set_ui(n, 0);
    char inputStr1[1024];
    printf("%s", msg);
    if(scanf("%1023s", inputStr1) == EOF) {
        printf("  ! Error: Scanning error\n\n");
        exit(-1);
    }
    if(!strcmp(inputStr1, "q")) {
        printf("  - Goodbye\n\n");
        exit(0);
    } else if(!strcmp(inputStr1, "r")) {
        uint64_t r = rand_uint64();
        printf("  - Factoring: %lu\n", r);
        mpz_set_ui(n, rand_uint64());
    } else if(!strcmp(inputStr1, "pq")) {
        uint64_t p = rand_uint32_prime();
        uint64_t q = rand_uint32_prime();
        uint64_t pq = p * q;
        printf("  - Factoring: %lu = %lu * %lu\n", pq, p, q);
        mpz_set_ui(n, pq);
    } else if(!strcmp(inputStr1, "p")) {
        uint64_t p = rand_uint32_prime();
        printf("  - Factoring: %lu\n", p);
        mpz_set_ui(n, p);
    } else {
        int flag = mpz_set_str(n, inputStr1, 10);
        if(flag) {
            printf("  ! Error: Invalid String\n\n");
            askBigInt(msg, n);
        } else if(mpz_sgn(n) <= 0 || mpz_log2(n) > 64) {
            printf("  ! Error: Out of range, 0 < N < 2 ^ 64\n\n");
            askBigInt(msg, n);
        }
    }
}

static void printBigInt(BigInt n) { //Print BigInt in decimal to stdout
    mpz_out_str(stdout,10,n);
}
#pragma endregion

#pragma region Table
static void atkins_sieve(uint32_t* ps) { //Compute primes below INT32_MAX using atkins sieve
    const uint64_t MAX = INT32_MAX;
    const uint64_t SQRT_MAX = (uint64_t) sqrt(MAX) + 1;
    bool *array = (bool*)calloc(MAX, sizeof(bool));
    
    for (uint64_t x = 1; x < SQRT_MAX; x++)
        for (uint64_t y = 1; y < SQRT_MAX; y++) {
            uint64_t k = 4 * x * x + y * y;
            if ((k < MAX) && ((k % 12 == 1) || (k % 12 == 5)))
                array[k] = !array[k];
            k = 3 * x * x + y * y;
            if ((k < MAX) && (k % 12 == 7))
                array[k] = !array[k];
            if (x > y) {
                k = 3 * x * x - y * y;
                if ((k < MAX) && (k % 12 == 11))
                    array[k] = !array[k];
            }
        }

    array[2] = true;
    array[3] = true;
    for (uint64_t n = 5; n <= SQRT_MAX; n++)
        if (array[n]) {
            uint64_t n2 = n * n;
            for (uint64_t k = n2; k < MAX; k += n2)
                array[k] = false;
        }

    int i = 0;
    for(uint64_t j = 0; j < MAX; j++)
        if(array[j])
            ps[i++] = j;

    free(array);
}

uint32_t count_bit(uint32_t x) { //Count set 1 bits in uint32_t
  x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
  x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
  x = (x & 0x0000FFFF) + ((x >> 16)& 0x0000FFFF);
  return x;
}

void combin_mask(int bound, uint32_t num_primes, bool **mask) { //Build mask from int i and ~i count_bit(i) == 15
    //(COMBIN(30, 15) = 155117520) > (num_primes = 105097564)
    for (int i = 0; i < bound; i++ ) 
        mask[i] = calloc(num_primes, sizeof(bool));


    int k = 15;
    int n = 30;
    int i = 0;
    int c = (1 << k) - 1;
    while (c < (1 << (n-1))) {
        uint32_t j = c;
        for(int l = 0; l < bound; l++) {
                if(j & 1) {
                    mask[l][i] = true;
                } else { 
                    mask[l][i + 1] = true;
                }
                j >>= 1;
        }
        i+=2;
        if(i >= num_primes) break;
        int a = c & (-c), b  = c + a;
        c = (c ^ b) / 4 / a | b;
    }
}

#ifdef THREADED
typedef struct {
    BigInt* k;
    BigInt* max;
} ThreadArgs;

void* thread_function(void* void_args) { //Compute table in thread
    ThreadArgs* thread_args = (ThreadArgs*)void_args;
    BigInt* start = thread_args->k;
    BigInt* max = thread_args->max;
    for(int j = 1; thread_args->k + j < max; j <<= 1)
        for(BigInt* k = thread_args->k; k + j < max; k += j << 1) {
            mpz_mul(k[0], k[0], k[j]);
            mpz_clear(k[j]);
        }
}
#endif

void build_table(const char* folder) { //Build table from product of masked primes
    #define TIME(start) (double)(clock() - start) / CLOCKS_PER_SEC
    clock_t overall_time = clock();
    uint64_t n_bound = UINT64_MAX;
    uint32_t primes_bound = UINT32_MAX;
    printf("Factoring Bound:%lu\n",n_bound);

    clock_t step_time = clock();
    printf("  - Generating primes < %u\n",primes_bound);
    const uint32_t num_primes = 105097564; //precalculated
    uint32_t *primes = malloc(num_primes * sizeof(uint32_t));
    atkins_sieve(primes);
    printf("    + Primes[%d] generated in %lfs\n", num_primes, TIME(step_time));

    step_time = clock();
    int bound = 30;
    printf("  - Populating mask(%d choose %d)\n", bound, bound / 2);
    bool *mask[bound];
    combin_mask(bound, num_primes, mask);
    printf("    + Masks[%d][%d] generated in %lfs\n", bound, num_primes, TIME(step_time));
    
    
    step_time = clock();
    printf("  - Masking primes, and computing product table[%d]\n", bound);
    const int row_len = num_primes / 2;
    double total_bits = 0;
    BigInt *row = calloc(row_len, sizeof(BigInt));
    int folder_len = strlen(folder);
    for(int i = 0; i < bound; i++) {
        printf("    + Entry(%d): Masking...", i); fflush(stdout); //Populate row with masked primes, free mask[i]
        clock_t row_time = clock();
        int k = 0;
        for(int j = 0; j < num_primes; j++)
            if(mask[i][j])
                mpz_init_set_ui(row[k++], primes[j]);
        free(mask[i]);

        printf("Multiplying..."); fflush(stdout); //Multiply pairs until 1 remains
#ifdef THREADED
        BigInt* ranges[5] = {row, row + 16777216, row + 33554432, row + 50331648, row + row_len };
        ThreadArgs thread_args[4];
        pthread_t threads[4];
        for(int j = 0; j < 4; j++)
                thread_args[j] = (ThreadArgs) {.k=ranges[j], .max = ranges[j+1]};
                
        for(int k = 0; k < 4; k++) {
            pthread_create(&threads[k], NULL, thread_function, (void*) &thread_args[k]);
        }
        for(int k = 0; k < 4; k++)
            pthread_join(threads[k], NULL);
        
        mpz_mul(row[0], row[0], row[16777216]);
        mpz_mul(row[33554432], row[33554432], row[50331648]);
        mpz_mul(row[0], row[0], row[33554432]);
        mpz_clear(row[16777216]);
        mpz_clear(row[33554432]);
        mpz_clear(row[50331648]);
#endif
#ifndef THREADED
        for(int j = 1; j < row_len; j*=2) 
            for(int k = 0; k + j < row_len; k += 2 * j) {
                mpz_mul(row[k], row[k], row[k + j]);
                mpz_clear(row[k + j]);
            }
#endif

        printf("Writing..."); fflush(stdout); //Grab file path, write to file
        int length = snprintf(NULL, 0, "%d", i);
        char filename[folder_len + length + 1]; 
            strncpy(filename, folder, sizeof(filename) - 1);
            snprintf((char*__restrict__)(filename + folder_len), length + 1, "%d", i);
        FILE *file = fopen(filename,"w");
        uint32_t error = !mpz_out_raw(file, row[0]);
        fclose(file);
        total_bits += mpz_log2(row[0]);
        mpz_clear(row[0]);
        if(error) {
            printf("\n    ! Error writing to file: %s\n", "0");
            break;  
        }
        printf("Done in %lfs\n", TIME(row_time)); fflush(stdout);

        if(i == bound - 1) { //Print last iteration message
            printf("    + Wrote to files: 0-%d after %lfs\n", bound - 1, TIME(step_time));
            printf("    + Total Bits:%lf", total_bits);
            printf("    + Average Bits:%lf", total_bits / bound);
        }
    }
    printf("    + Table[%d] populated in %lfs\n", bound, TIME(step_time));

    printf("  - Cleaning up\n"); //Free rows and primes.
    step_time = clock();
    free(row);
    free(primes);
    printf("    + Cleaned up after %lfs", TIME(step_time));

    printf("Finished after %lfs\n", TIME(overall_time));
}
#pragma endregion

#pragma region Factor
void shuffle(int *array, size_t n) {
    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; i++) {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

void factor(const char* folder, int tries) { //Randomly select entry in table and compute gcd up to tries times
    BigInt n, filter;
    askBigInt("Factor what uint64_t? (r:RAND, p:PRIME, pq:P*Q, q:QUIT) : ", n);
    mpz_init(filter);

    int folder_len = strlen(folder);

    int start_tries = tries;
    int indices[30];
    for(int i = 0; i < tries; i++) 
        indices[i] = i;
    shuffle(indices, 30);
    printf("  - Opened files:");
    while(tries-- > 0) {
        int length = snprintf(NULL, 0, "%d", indices[tries]);
        char filename[folder_len + length + 1]; 
            strncpy(filename, folder, sizeof(filename) - 1);
            snprintf((char*__restrict__)(filename + folder_len), length + 1, "%d", indices[tries]);
        printf(" %d", indices[tries]); fflush(stdout);

        FILE *file = fopen(filename,"r");
        int i = mpz_inp_raw(filter, file);
        fclose(file);
        if(i) {
            BigInt gcd;
            mpz_init(gcd);
            mpz_gcd(gcd, filter, n);
            if(mpz_cmp(gcd, n) != 0 && mpz_cmp_ui(gcd, 1) != 0) {
                BigInt factor; mpz_init(factor); mpz_divexact(factor, n, gcd);
                printf("\n  - N is composite: "); printBigInt(n); printf(" = "); 
                    printBigInt(gcd); printf(" * "); printBigInt(factor); printf("\n\n");
                //uint64_t gcd_t = mpz_get_ui(gcd);
                //uint64_t factor_t = mpz_get_ui(factor);
                mpz_clear(gcd);
                mpz_clear(factor);
                mpz_clear(filter);
                mpz_clear(n);
                return;
            } else if(tries != 0) {
                printf(",");
            }
            mpz_clear(gcd);
        } else {
            printf("\n  ! File Error:%d\n",i);
            exit(-1);
        }
    }
    if(start_tries == 30) 
        printf("\n  - N is prime\n\n");
    else 
        printf("\n  - N is possibly prime\n\n");
    mpz_clear(filter);
    mpz_clear(n);
}
void usage() {
    printf("Usage: ./factor [-build [TBL_DIR=\"table/\"] | -factor [TRIES=30 [TBL_DIR=\"table/\"]]] \n");
    exit(-1);
}
#pragma endregion

int main(int argc, char *argv[]) {
    if(argc == 1 || argc > 4) {
        usage();
    } 

    if(!strcmp(argv[1], "-build")) {
        if(argc == 2) {
            build_table("table/");
        } else if(argc == 3) {
            build_table(argv[2]);
        } else {
            usage();
        }
    } else if (!strcmp(argv[1], "-factor")) {
        time_t t; srand((unsigned) time(&t));
        if(argc == 2) {
            while(true) 
                factor("table/", 30);
        } else if(argc == 3 || argc == 4) {
            int tries = atoi(argv[2]);
            if(tries <= 0 || 30 < tries) {
                printf("Invalid TRIES, 0 < TRIES <= 30)");
                exit(-1);
            }
            while(true) 
                factor(argc == 3 ? "table/" : argv[4], tries);
        } else {
            usage();
        }
    } else {
        usage();
    }
    return 0;
}
