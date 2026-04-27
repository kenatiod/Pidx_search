// primorial_minPidx.c
/*
This program goes through intervals of n between primorials and computes the
pi(GPF) of the product n(n+1).
The prime index of the greatest prime factor, pi(GPF), is called the "Pidx"
in this code. Pidx of the product is the larger Pidx value of the
n and n+1 pair. Computing those values for n and n+1 avoids forming the product and
computing over larger numbers. Also, progressing through the numbers, the values for
n+1 can be saved for use in the product (n+1)(n+2), so successive numbers are not
unnecessarily processed twice.

For each interval [primorial_roots[r], primorial_roots[r+1]) the processing is looking
for the minimum Pidx for that interval. The individual values of the
Pidx jump rapidly each n, but the minimum Pidx grows interval to interval as a trend,
but sags back before surging forward. The maximum omega side starts out greater than the
minimum Pidx, but over many intervals, the minimum Pidx catches up, and then exceeds
maximum omega.

A number with a complete set of prime factors, from its GPF down to 2, is called
"prime-complete" and the products of consecutive integers that are prime-complete are
listed in sequence A141399 at OEIS.com and stop at 633,555x633,556. Numbers are
prime-complete if, and only if, their omega exactly equals their Pidx. For the products
of pairs of consecutive integers, omega grows approximately as 2*ln(ln(n)) while Pidx
grows like a positive power of n divided by ln(n).

The significance of looking at the maximum omega and minimum Pidx is that if the minimum
Pidx is greater than the maximum omega over an interval, then it is impossible for any
point in the interval to have had omega equal to Pidx, thus impossible for any value to be
prime-complete. Once minimum Pidx grows to exceed maximum omega in successive intervals,
the closed set of prime-complete products has been enumerated, and the A141399 sequence
is terminated. The maximum omega is always less than or equal to the r of the primorial 
interval.

On Linux, build with:
    gcc -O3 -std=c11 -fopenmp -o primorial_minPidx primorial_minPidx.c -lm

On Mac you might need:
    brew install libomp
    clang -O3 -std=c11 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
      -o primorial_minPidx primorial_minPidx.c \
      -L/opt/homebrew/opt/libomp/lib -lomp -lm

Run with:
    ./primorial_minPidx [start interval number] [how many intervals to search]

Default run:
    ./primorial_minPidx 1 22

By Ken Clements, April 20, 2026      
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

#ifdef _OPENMP
#include <omp.h>
#else
static int    omp_get_max_threads(void) { return 1; }
static int    omp_get_thread_num(void)  { return 0; }
static double omp_get_wtime(void)       { return 0.0; }
#endif

#define CHECKPOINT_FILE  "primorial_minPidx_checkpoint.dat"
#define DEFAULT_BLOCK_N  524288
#define SENTINEL_PIDX    9999999
/* Number of fields written/read by checkpoint I/O — keep in sync. */
#define CHECKPOINT_FIELDS 7

/* =========================================================================
 * Prime table: primes[1..100] = 2, 3, 5, ..., 541
 * primes[0] is unused (set to 0).
 * ========================================================================= */
static uint64_t primes[101] = {
    0,   2,   3,   5,   7,  11,  13,  17,  19,  23,
   29,  31,  37,  41,  43,  47,  53,  59,  61,  67,
   71,  73,  79,  83,  89,  97, 101, 103, 107, 109,
  113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
  173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
  229, 233, 239, 241, 251, 257, 263, 269, 271, 277,
  281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
  349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
  409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
  463, 467, 479, 487, 491, 499, 503, 509, 521, 523,
  541
};
static int max_pidx = 100;

/* =========================================================================
 * Data structures
 * ========================================================================= */
typedef struct {
    int      pidx_small;   /* index of largest small prime factor found      */
    int      omega_small;  /* count of distinct small prime factors found     */
    uint64_t residue;      /* cofactor remaining after small-prime trial div  */
    int      exact_pidx;   /* 1 if pidx_small is the true GPF index           */
    int      exact_omega;  /* 1 if omega_small is the true omega              */
} po_result_t;

typedef struct {
    int      max_omega;
    uint64_t mo_n;
    int      min_pidx;
    uint64_t mp_n;
    /* diagnostic counters (currently unused but kept for future use) */
    uint64_t near_max_pairs;
    uint64_t residue_prime_checks;
    uint64_t square_checks;
    uint64_t rho_calls;
    uint64_t full_factor_calls;
} thread_result_t;

typedef struct {
    int      next_loop_index;
    uint64_t next_n_start;      /* = primorial_roots[next_loop_index+1]       */
    uint64_t next_interval_size;
    uint64_t orig_interval_start;
    int      orig_intervals;
    int      effective_intervals;
    int      max_pidx_saved;
    uint64_t block_n;
} checkpoint_t;

/* =========================================================================
 * splitmix64 PRNG (used for future Pollard-rho seeding)
 * ========================================================================= */
static uint64_t splitmix64_next(uint64_t *state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}



static inline int pidx_exact(uint64_t n, int r) {
    int pidx = 0;
    int limit_idx = r + 3;
    if (limit_idx > max_pidx) limit_idx = max_pidx;

    /* Strip all factors of 2 in one instruction */
    if (n % 2 == 0) {
        n >>= __builtin_ctzll(n);
        pidx = 1;
        if (n == 1) return 1;
    }

    /* Trial divide by odd primes[2..limit_idx]  */
    for (int i = 2; i <= limit_idx; i++) {
        uint64_t p = primes[i];
        if (n % p == 0) {
            pidx = i;
            do { n /= p; } while (n % p == 0);
            if (n == 1) return pidx;
        }
    }
    return 9999;
}


/* =========================================================================
 * isqrt128 — integer square root of a 128-bit unsigned integer.
 * Returns floor(sqrt(n)).  Result always fits in 64 bits.
 * ========================================================================= */
uint64_t isqrt128(__uint128_t n) {
    if (n == 0) return 0;

    /* Phase 1: long-double seed */
    __uint128_t x = (uint64_t)sqrtl((long double)n);
    if (x > (uint64_t)UINT64_MAX) x = (uint64_t)UINT64_MAX;

    /* Phase 2: Newton-Raphson refinement */
    __uint128_t x1;
    while (1) {
        x1 = (x + n / x) / 2;
        if (x1 >= x) break;
        x = x1;
    }

    /* Phase 3: Correct any ceil overshoot (at most 1 step) */
    while (x * x > n) x--;

    return (uint64_t)x;
}

/* =========================================================================
 * Formatting helpers
 * ========================================================================= */
static char *fmt_num_u64(uint64_t n, char *buf, size_t size) {
    char temp[64];
    snprintf(temp, sizeof(temp), "%" PRIu64, n);
    int len = (int)strlen(temp);
    int out_idx = 0;
    int first_group = ((len - 1) % 3) + 1;
    for (int i = 0; i < len; i++) {
        if (i > 0 && i >= first_group && ((i - first_group) % 3) == 0) {
            if ((size_t)out_idx < size - 1) buf[out_idx++] = ',';
        }
        if ((size_t)out_idx < size - 1) buf[out_idx++] = temp[i];
    }
    buf[out_idx] = '\0';
    return buf;
}

static void uint128_to_str(__uint128_t n, char *str) {
    if (n == 0) { str[0] = '0'; str[1] = '\0'; return; }
    char temp[40];
    int i = 0;
    while (n > 0) { temp[i++] = (char)(n % 10 + '0'); n /= 10; }
    for (int j = 0; j < i; j++) str[j] = temp[i - 1 - j];
    str[i] = '\0';
}

/* Convert __uint128_t to comma-formatted decimal string.
 * buf must be at least 64 bytes.  Returns buf, or NULL if too small. */
static char *uint128_to_commastr(__uint128_t n, char *buf, size_t bufsize) {
    char tmp[64], out[64];
    int tlen = 0, olen = 0;

    if (n == 0) {
        tmp[tlen++] = '0';
    } else {
        while (n > 0) { tmp[tlen++] = '0' + (int)(n % 10); n /= 10; }
    }
    for (int i = 0; i < tlen; i++) {
        if (i > 0 && i % 3 == 0) out[olen++] = ',';
        out[olen++] = tmp[i];
    }
    if ((size_t)(olen + 1) > bufsize) return NULL;
    for (int i = 0; i < olen; i++) buf[i] = out[olen - 1 - i];
    buf[olen] = '\0';
    return buf;
}

static void fmt_time(double sec, char *buf, size_t size) {
    if      (sec < 120.0)    snprintf(buf, size, "%6.2fs", sec);
    else if (sec < 7200.0)   snprintf(buf, size, "%6.2fm", sec / 60.0);
    else if (sec < 172800.0) snprintf(buf, size, "%6.2fh", sec / 3600.0);
    else                     snprintf(buf, size, "%6.2fd", sec / 86400.0);
}

/* =========================================================================
 * Checkpoint I/O
 * The struct has CHECKPOINT_FIELDS = 7 fields; keep in sync if you add more.
 * ========================================================================= */
static void write_checkpoint(const checkpoint_t *cp) {
    FILE *f = fopen(CHECKPOINT_FILE, "w");
    if (!f) { fprintf(stderr, "Warning: could not save checkpoint\n"); return; }
    fprintf(f,
        "orig_interval_start=%" PRIu64 "\n"
        "orig_intervals=%d\n"
        "effective_intervals=%d\n"
        "block_n=%" PRIu64 "\n"
        "next_loop_index=%d\n"
        "next_n_start=%" PRIu64 "\n"
        "next_interval_size=%" PRIu64 "\n",
        cp->orig_interval_start,
        cp->orig_intervals,
        cp->effective_intervals,
        cp->block_n,
        cp->next_loop_index,
        cp->next_n_start,
        cp->next_interval_size);
    fclose(f);
}

static int read_checkpoint(uint64_t orig_interval_start, int orig_intervals,
                            uint64_t block_n, checkpoint_t *cp) {
    FILE *f = fopen(CHECKPOINT_FILE, "r");
    if (!f) return 0;
    memset(cp, 0, sizeof(*cp));
    char line[256];
    int fields = 0;
    unsigned long long ull;
    int iv;
    while (fgets(line, sizeof(line), f)) {
        if      (sscanf(line, "orig_interval_start=%" SCNu64, &cp->orig_interval_start) == 1) fields++;
        else if (sscanf(line, "orig_intervals=%d",  &iv) == 1) { cp->orig_intervals = iv; fields++; }
        else if (sscanf(line, "effective_intervals=%d", &iv) == 1) { cp->effective_intervals = iv; fields++; }
        else if (sscanf(line, "block_n=%" SCNu64, &cp->block_n) == 1) fields++;
        else if (sscanf(line, "next_loop_index=%d", &iv) == 1) { cp->next_loop_index = iv; fields++; }
        else if (sscanf(line, "next_n_start=%" SCNu64, &cp->next_n_start) == 1) fields++;
        else if (sscanf(line, "next_interval_size=%" SCNu64, &cp->next_interval_size) == 1) fields++;
    }
    fclose(f);
    /* BUG FIX: was "!= 10"; struct has exactly CHECKPOINT_FIELDS = 7 fields */
    if (fields != CHECKPOINT_FIELDS)            return 0;
    if (cp->orig_interval_start != orig_interval_start) return 0;
    if (cp->orig_intervals      != orig_intervals)      return 0;
    if (cp->block_n             != block_n)             return 0;
    return 1;
}

/* =========================================================================
 * main
 * ========================================================================= */
int main(int argc, char *argv[]) {
    uint64_t interval_start = 1;
    uint64_t intervals      = 22;
    uint64_t block_n        = DEFAULT_BLOCK_N;

    if (argc > 1) interval_start = strtoull(argv[1], NULL, 10);
    if (argc > 2) intervals      = strtoull(argv[2], NULL, 10);

    uint64_t orig_interval_start = interval_start;
    int      orig_intervals      = (int)intervals;

    /* --- Build primorials and their integer square roots --- */
    __uint128_t primorials[30];
    uint64_t    primorial_roots[30];
    primorials[0]      = 1;
    primorial_roots[0] = 1;
    for (int idx = 1; idx < 30; idx++) {
        primorials[idx]      = primorials[idx-1] * (__uint128_t)primes[idx];
        primorial_roots[idx] = isqrt128(primorials[idx]);
    }

    /* --- Checkpoint resume --- */
    int      resuming           = 0;
    int      first_interval     = 0;
    int      effective_intervals = (int)intervals;
    checkpoint_t cp;

    if (read_checkpoint(orig_interval_start, orig_intervals, block_n, &cp)) {
        resuming            = 1;
        first_interval      = cp.next_loop_index;
        effective_intervals = cp.effective_intervals;
        /* n boundaries are re-derived from primorial_roots, not from checkpoint */
    }

    int nthreads = omp_get_max_threads();

    printf("\n primorial_minPidx: Find minPidx over Primorial Intervals\n");
    printf("============================================================\n\n");
    if (resuming) {
        printf("*** Resuming from checkpoint: starting at loop interval %d ***\n\n",
               first_interval);
    }
    printf("Tracking: min(Pidx) per interval.\n");
    printf("Intervals: %d between primorial boundaries of n(n+1).\n", orig_intervals);
    printf("Threads: %d,  block_n: %" PRIu64 "\n\n", nthreads, block_n);
    printf("  r  %26s  minPidx %18s  %8s\n",
            "Interval Start n,", "n@minPidx", "time");
    printf("-----------------------------------------------------------------------\n");
    fflush(stdout);

    /* --- Handle the trivial r=1 row when starting from the beginning --- */
    if (!resuming && orig_interval_start == 1) {
        printf("  1 %26s  %6d %18s   %8s\n",
                "1", 1, "1", "  0.00s");
        first_interval = 1;
        effective_intervals = (int)intervals - 1;
    }

    /* --- CSV output --- */
    FILE *csv = fopen("primorial_minPidx_data.csv", resuming ? "a" : "w");
    if (csv) {
        if (!resuming) {
            fprintf(csv, "interval,n_start,n_stop,max_omega,mo_n,min_pidx,mp_n\n");
            if (orig_interval_start == 1)
                fprintf(csv, "1,1,1,1,1,1,1\n");
        }
    }

    /* =====================================================================
     * Main interval loop.
     *
     * BUG FIX (primary): intervals are now anchored directly to the
     * primorial roots, not computed via cumulative n_stop arithmetic.
     *
     * Interval r (loop variable "interval") covers:
     *   n in [ primorial_roots[interval+1], primorial_roots[interval+2] )
     * so that n*(n+1) falls between p_{r}# and p_{r+1}#.
     * ===================================================================== */
    for (int interval = first_interval; interval <= effective_intervals; interval++) {

        /* BUG FIX: anchor start and stop directly to primorial_roots */
        uint64_t n_start = primorial_roots[interval + 1];
        uint64_t n_stop  = primorial_roots[interval + 2];  /* exclusive */

        double t0 = omp_get_wtime();
        thread_result_t *thr = (thread_result_t *)calloc((size_t)nthreads, sizeof(thread_result_t));
        if (!thr) { fprintf(stderr, "Allocation failure for thread results\n"); return 1; }
        for (int t = 0; t < nthreads; t++) {
            thr[t].max_omega = -1;
            thr[t].min_pidx  = SENTINEL_PIDX;
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int             tid        = omp_get_thread_num();
            thread_result_t *my        = &thr[tid];
            /* seed_state kept for future Pollard-rho use */
            uint64_t seed_state = 0x123456789abcdef0ULL
                                  ^ ((uint64_t)(tid + 1) * 0x9E3779B97F4A7C15ULL);
            (void)seed_state; /* suppress unused-variable warning for now */

            uint64_t total_numbers = n_stop - n_start;
            uint64_t num_blocks    = (total_numbers + block_n - 1) / block_n;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (uint64_t blk = 0; blk < num_blocks; blk++) {
                uint64_t block_start     = n_start + blk * block_n;
                /* BUG FIX: clamp to n_stop (exclusive), not n_stop+1 */
                uint64_t block_end_excl  = block_start + block_n;
                if (block_end_excl > n_stop) block_end_excl = n_stop;

                size_t len = (size_t)(block_end_excl - block_start);
                /* Allocate len+1 so we can always form the pair (n, n+1)
                 * for the last n = block_end_excl - 1, whose n+1 may be
                 * the first element of the next block. */
                po_result_t *arr = (po_result_t *)malloc((len + 1) * sizeof(po_result_t));
                if (!arr) { fprintf(stderr, "Allocation failure for block array\n"); exit(1); }

                for (size_t i = 0; i <= len; i++) {
//                     arr[i] = find_po_residue(block_start + (uint64_t)i, max_pidx);
//                     arr[i] = find_po_residue(block_start + (uint64_t)i, interval + 10);

                    arr[i].pidx_small = pidx_exact(block_start + (uint64_t) i, interval);
                    arr[i].omega_small = 0;
                    arr[i].residue    = 0;
                    arr[i].exact_pidx = (arr[i].pidx_small < 999);
                    arr[i].exact_omega = 0;
                }

                /* Pairs: (block_start+i) and (block_start+i+1), for i in [0, len-1] */
                for (size_t i = 0; i < len; i++) {
                    uint64_t      n = block_start + (uint64_t)i;
                    po_result_t  *a = &arr[i];
                    po_result_t  *b = &arr[i + 1];

                    /* --- min Pidx: only valid when BOTH sides are fully factored --- */
                    if (a->exact_pidx && b->exact_pidx) {
                        int pidx = (a->pidx_small > b->pidx_small)
                                    ? a->pidx_small : b->pidx_small;
                        if (pidx < my->min_pidx ||
                            (pidx == my->min_pidx && n > my->mp_n)) {
                            my->min_pidx = pidx;
                            my->mp_n     = n;
                        }
                    }

                    /* --- max omega: only valid when BOTH sides give exact omega --- */
                    if (a->exact_omega && b->exact_omega) {
                        int omega = a->omega_small + b->omega_small;
                        if (omega > my->max_omega ||
                            (omega == my->max_omega && n > my->mo_n)) {
                            my->max_omega = omega;
                            my->mo_n      = n;
                        }
                    }
                }

                free(arr);
            } /* blk loop */
        } /* omp parallel */

        /* --- Reduce thread results --- */
        thread_result_t best;
        memset(&best, 0, sizeof(best));
        best.max_omega = -1;
        best.min_pidx  = SENTINEL_PIDX;
        for (int t = 0; t < nthreads; t++) {
            if (thr[t].max_omega > best.max_omega ||
                (thr[t].max_omega == best.max_omega && thr[t].mo_n > best.mo_n)) {
                best.max_omega = thr[t].max_omega;
                best.mo_n      = thr[t].mo_n;
            }
            if (thr[t].min_pidx < best.min_pidx ||
                (thr[t].min_pidx == best.min_pidx && thr[t].mp_n > best.mp_n)) {
                best.min_pidx = thr[t].min_pidx;
                best.mp_n     = thr[t].mp_n;
            }
        }
        free(thr);

        /* --- Print row --- */
        double dt = omp_get_wtime() - t0;
        char primorial_str[128], start_str[64], nxnp1_str[128],
             mo_str[64], mp_str[64], time_str[32];

        uint128_to_commastr(primorials[interval + 1], primorial_str, sizeof(primorial_str));
        fmt_num_u64(n_start, start_str, sizeof(start_str));
        /* n(n+1) displayed uses __uint128_t to avoid overflow */
        uint128_to_commastr((__uint128_t)n_start * (n_start + 1),
                            nxnp1_str, sizeof(nxnp1_str));
        fmt_num_u64(best.mo_n, mo_str, sizeof(mo_str));
        fmt_num_u64(best.mp_n, mp_str, sizeof(mp_str));
        fmt_time(dt, time_str, sizeof(time_str));

        printf("%3d %26s  %6d %18s   %8s\n",
               interval + 1,
               start_str,
               best.min_pidx == SENTINEL_PIDX ? 9999 : best.min_pidx,
               mp_str,
               time_str);
        fflush(stdout);

        if (csv) {
            fprintf(csv, "%d,%" PRIu64 ",%" PRIu64 ",%d,%" PRIu64 ",%d,%" PRIu64 "\n",
                    interval + 1,
                    n_start,
                    n_stop - 1,
                    best.max_omega,
                    best.mo_n,
                    best.min_pidx == SENTINEL_PIDX ? 9999 : best.min_pidx,
                    best.mp_n);
            fflush(csv);
        }

        /* --- Save checkpoint --- */
        checkpoint_t outcp;
        outcp.next_loop_index      = interval + 1;
        outcp.next_n_start         = n_stop;        /* = primorial_roots[interval+2] */
        outcp.next_interval_size   = n_stop - n_start;
        outcp.orig_interval_start  = orig_interval_start;
        outcp.orig_intervals       = orig_intervals;
        outcp.effective_intervals  = effective_intervals;
        outcp.max_pidx_saved       = max_pidx;
        outcp.block_n              = block_n;
        write_checkpoint(&outcp);

    } /* interval loop */

    if (csv) fclose(csv);
    remove(CHECKPOINT_FILE);
    printf("\nDone.\n");
    return 0;
}
