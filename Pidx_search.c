/*
 * Pidx_search.c
 *
 * Search products m(m+1) over intervals defined by successive square roots of
 * primorials.  For each r, the program scans
 *
 *   ceil(sqrt(p_r#)) <= m < ceil(sqrt(p_{r+1}#))
 *
 * and reports interval statistics in the same general format as Delta_min.c.
 *
 * IMPORTANT LIMITATION:
 *   This program is optimized to advance min_Pidx quickly.  It factors each m
 *   and m+1 only up to the current best interval min_Pidx threshold.  As a
 *   result, min_Pidx is exact, while minDelta and maxOmega are only recorded on
 *   pairs that factor completely within the active threshold at the time they
 *   are examined.  Those two columns are therefore auxiliary / opportunistic.
 *
 * Compile (Linux):
 *   gcc -O3 -std=c11 -march=native -fopenmp -o Pidx_search Pidx_search.c -lm
 *
 * Compile (macOS / Homebrew llvm):
 *   clang -O3 -std=c11 -march=native -Xpreprocessor -fopenmp -o Pidx_search Pidx_search.c -lomp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_max_threads(void) { return 1; }
static int omp_get_thread_num(void) { return 0; }
static double omp_get_wtime(void) { return 0.0; }
#endif

static const uint64_t primes[101] = {0, 2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
    113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
    173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277,
    281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
    349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
    409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
    463, 467, 479, 487, 491, 499, 503, 509, 521, 523,
    541};
static const int primes_length = 100;

static volatile int g_best_min_pidx = INT_MAX;

typedef struct {
    int min_delta;
    int md_pidx;
    int md_omega;
    uint64_t md_n;

    int max_omega;
    uint64_t mo_n;

    int min_pidx;
    uint64_t mp_n;

    uint64_t exact_pairs;
    uint64_t pruned_left;
    uint64_t pruned_right;
} thread_result_t;

typedef struct {
    int pidx;
    int omega;
    uint64_t residue;
} po_result_t;

static char* fmt_num_u64(uint64_t n, char *buf, size_t size) {
    char temp[64];
    snprintf(temp, sizeof(temp), "%llu", (unsigned long long)n);
    int len = (int)strlen(temp);
    int out_idx = 0;
    int first_group = ((len - 1) % 3) + 1;
    for (int i = 0; i < len; i++) {
        if (i > 0 && ((i - first_group) % 3) == 0 && i >= first_group) {
            if (out_idx < (int)size - 1) buf[out_idx++] = ',';
        }
        if (out_idx < (int)size - 1) buf[out_idx++] = temp[i];
    }
    buf[out_idx] = '\0';
    return buf;
}

static void fmt_time(double sec, char *buf, size_t size) {
    if (sec < 120.0) snprintf(buf, size, "%6.2fs", sec);
    else if (sec < 7200.0) snprintf(buf, size, "%6.2fm", sec / 60.0);
    else if (sec < 172800.0) snprintf(buf, size, "%6.2fh", sec / 3600.0);
    else snprintf(buf, size, "%6.2fd", sec / 86400.0);
}

static uint64_t ceil_sqrt_u128(unsigned __int128 x) {
    long double xd = (long double)x;
    long double rd = sqrtl(xd);
    uint64_t r = (uint64_t)rd;
    while ((unsigned __int128)r * (unsigned __int128)r < x) r++;
    while (r > 0 && (unsigned __int128)(r - 1) * (unsigned __int128)(r - 1) >= x) r--;
    return r;
}

static inline int factor_to_limit(uint64_t n, int limit_idx, po_result_t *out) {
    int pidx = 0, omega = 0;
    uint64_t residue = n;

    for (int i = 1; i <= limit_idx; i++) {
        uint64_t p = primes[i];
        if (p > residue) break;
        if (residue % p == 0) {
            omega++;
            pidx = i;
            do { residue /= p; } while (residue % p == 0);
            if (residue == 1) break;
        }
    }

    out->pidx = pidx;
    out->omega = omega;
    out->residue = residue;
    return residue == 1;
}

int main(int argc, char *argv[]) {
    int r_start = 5;
    int r_end = 22;
    int max_pidx = primes_length;

    if (argc > 1) r_start = atoi(argv[1]);
    if (argc > 2) r_end = atoi(argv[2]);
    if (argc > 3) max_pidx = atoi(argv[3]);

    if (r_start < 1) r_start = 1;
    if (r_end > primes_length - 1) r_end = primes_length - 1;
    if (r_start > r_end) {
        fprintf(stderr, "Invalid r range.\n");
        return 1;
    }
    if (max_pidx < 1) max_pidx = 1;
    if (max_pidx > primes_length) max_pidx = primes_length;
    if (r_end + 1 > max_pidx) max_pidx = r_end + 1;
    if (max_pidx > primes_length) max_pidx = primes_length;

    unsigned __int128 primorials[101] = {0};
    uint64_t sqrt_bounds[101] = {0};
    primorials[0] = 1;
    for (int i = 1; i <= max_pidx; i++) {
        unsigned __int128 next = primorials[i - 1] * (unsigned __int128)primes[i];
        primorials[i] = next;
        sqrt_bounds[i] = ceil_sqrt_u128(next);
        if (sqrt_bounds[i] == 0) {
            fprintf(stderr, "Overflow while building sqrt(primorial) bounds at r=%d.\n", i);
            return 1;
        }
    }

    FILE *csv = fopen("Pidx_search_data.csv", "w");
    if (csv) {
        fprintf(csv, "r,n_start,n_stop,min_delta,md_pidx,md_omega,md_n,max_omega,mo_n,min_pidx,mp_n,exact_pairs,pruned_left,pruned_right\n");
    }

    printf("\nPidx_search: sqrt(primorial)-band search for advancing min_Pidx\n");
    printf("===============================================================================\n\n");
    printf("Note: min_Pidx is exact. minDelta and maxOmega are opportunistic exact-only\n");
    printf("      values from pairs that fully factor inside the active min_Pidx cutoff.\n\n");
    printf("  r        Interval Start  minDelta = Pidx - omega      n@minDelta  maxOmega      n@maxOmega  minPidx        n@minPidx     time\n");
    printf("----------------------------------------------------------------------------------------------------------------------------------------\n");

    int nthreads = omp_get_max_threads();

    for (int r = r_start; r <= r_end; r++) {
        uint64_t n_start = sqrt_bounds[r];
        uint64_t n_stop_excl = sqrt_bounds[r + 1];
        if (n_stop_excl <= n_start) n_stop_excl = n_start + 1;

        double t0 = omp_get_wtime();
        g_best_min_pidx = INT_MAX;

        thread_result_t *thr = calloc((size_t)nthreads, sizeof(thread_result_t));
        if (!thr) {
            fprintf(stderr, "Allocation failed.\n");
            if (csv) fclose(csv);
            return 1;
        }
        for (int t = 0; t < nthreads; t++) {
            thr[t].min_delta = INT_MAX;
            thr[t].max_omega = -1;
            thr[t].min_pidx = INT_MAX;
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_result_t *my = &thr[tid];

            #pragma omp for schedule(dynamic, 65536)
            for (uint64_t n = n_start; n < n_stop_excl; n++) {
                int best_now;
                #pragma omp atomic read
                best_now = g_best_min_pidx;

                int cutoff = max_pidx;
                if (best_now != INT_MAX) {
                    cutoff = best_now - 1;
                    if (cutoff < 1) cutoff = 1;
                    if (cutoff > max_pidx) cutoff = max_pidx;
                }

                po_result_t a, b;
                int exact_a = factor_to_limit(n, cutoff, &a);
                if (!exact_a) {
                    my->pruned_left++;
                    continue;
                }
                int exact_b = factor_to_limit(n + 1, cutoff, &b);
                if (!exact_b) {
                    my->pruned_right++;
                    continue;
                }

                my->exact_pairs++;
                int pidx = (a.pidx > b.pidx) ? a.pidx : b.pidx;
                int omega = a.omega + b.omega;
                int delta = pidx - omega;

                if (pidx < my->min_pidx || (pidx == my->min_pidx && n > my->mp_n)) {
                    my->min_pidx = pidx;
                    my->mp_n = n;
                }
                if (delta < my->min_delta || (delta == my->min_delta && n > my->md_n)) {
                    my->min_delta = delta;
                    my->md_pidx = pidx;
                    my->md_omega = omega;
                    my->md_n = n;
                }
                if (omega > my->max_omega || (omega == my->max_omega && n > my->mo_n)) {
                    my->max_omega = omega;
                    my->mo_n = n;
                }

                int updated = 0;
                while (!updated) {
                    int cur;
                    #pragma omp atomic read
                    cur = g_best_min_pidx;
                    if (pidx >= cur) break;
                    #pragma omp critical
                    {
                        if (pidx < g_best_min_pidx) g_best_min_pidx = pidx;
                    }
                    updated = 1;
                }
            }
        }

        thread_result_t best;
        memset(&best, 0, sizeof(best));
        best.min_delta = INT_MAX;
        best.max_omega = -1;
        best.min_pidx = INT_MAX;

        for (int t = 0; t < nthreads; t++) {
            if (thr[t].min_delta < best.min_delta ||
                (thr[t].min_delta == best.min_delta && thr[t].md_n > best.md_n)) {
                best.min_delta = thr[t].min_delta;
                best.md_pidx = thr[t].md_pidx;
                best.md_omega = thr[t].md_omega;
                best.md_n = thr[t].md_n;
            }
            if (thr[t].max_omega > best.max_omega ||
                (thr[t].max_omega == best.max_omega && thr[t].mo_n > best.mo_n)) {
                best.max_omega = thr[t].max_omega;
                best.mo_n = thr[t].mo_n;
            }
            if (thr[t].min_pidx < best.min_pidx ||
                (thr[t].min_pidx == best.min_pidx && thr[t].mp_n > best.mp_n)) {
                best.min_pidx = thr[t].min_pidx;
                best.mp_n = thr[t].mp_n;
            }
            best.exact_pairs += thr[t].exact_pairs;
            best.pruned_left += thr[t].pruned_left;
            best.pruned_right += thr[t].pruned_right;
        }
        free(thr);

        double dt = omp_get_wtime() - t0;

        char start_str[64], md_n_str[64], mo_n_str[64], mp_n_str[64], time_str[32];
        fmt_num_u64(n_start, start_str, sizeof(start_str));
        fmt_num_u64(best.md_n, md_n_str, sizeof(md_n_str));
        fmt_num_u64(best.mo_n, mo_n_str, sizeof(mo_n_str));
        fmt_num_u64(best.mp_n, mp_n_str, sizeof(mp_n_str));
        fmt_time(dt, time_str, sizeof(time_str));

        printf("%3d %24s %9d %5d %5d %18s %6d %18s %6d %18s %8s\n",
               r,
               start_str,
               best.min_delta == INT_MAX ? 9999 : best.min_delta,
               best.md_pidx,
               best.md_omega,
               md_n_str,
               best.max_omega < 0 ? -1 : best.max_omega,
               mo_n_str,
               best.min_pidx == INT_MAX ? 9999 : best.min_pidx,
               mp_n_str,
               time_str);
        fflush(stdout);

        if (csv) {
            fprintf(csv, "%d,%llu,%llu,%d,%d,%d,%llu,%d,%llu,%d,%llu,%llu,%llu,%llu\n",
                    r,
                    (unsigned long long)n_start,
                    (unsigned long long)(n_stop_excl - 1),
                    best.min_delta == INT_MAX ? 9999 : best.min_delta,
                    best.md_pidx,
                    best.md_omega,
                    (unsigned long long)best.md_n,
                    best.max_omega < 0 ? -1 : best.max_omega,
                    (unsigned long long)best.mo_n,
                    best.min_pidx == INT_MAX ? 9999 : best.min_pidx,
                    (unsigned long long)best.mp_n,
                    (unsigned long long)best.exact_pairs,
                    (unsigned long long)best.pruned_left,
                    (unsigned long long)best.pruned_right);
            fflush(csv);
        }
    }

    if (csv) fclose(csv);
    printf("\nDone.\n");
    return 0;
}
