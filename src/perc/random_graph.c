#include "hp/perc/random_graph.h"

#include <math.h>
#include <string.h>

#include "hp/core/memory.h"

double hp_perc_er_giant_theory(double mean_degree) {
    const double c = mean_degree;
    if (!(c > 1.0) || !isfinite(c)) return 0.0;
    /* g(S) = S - 1 + exp(-c S); g(0)=0 is the trivial root, g(1)>0.
     * For c > 1 there is a root in (0,1) with g < 0 just above zero. */
    double lo = 1e-12, hi = 1.0;
    double glo = lo - 1.0 + exp(-c * lo);
    if (!(glo < 0.0)) return 0.0;
    for (int it = 0; it < 200; ++it) {
        double mid = 0.5 * (lo + hi);
        double g = mid - 1.0 + exp(-c * mid);
        if (g < 0.0)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

hp_status hp_perc_config_giant(const double *pk, size_t k_max,
                               double transmissibility, double *giant,
                               double *r0, double *mean_degree) {
    if (!pk || !(transmissibility >= 0.0 && transmissibility <= 1.0))
        return HP_ERR_INVALID;
    double sum = 0.0, mean = 0.0, second = 0.0;
    for (size_t k = 0; k <= k_max; ++k) {
        double p = pk[k];
        if (!isfinite(p) || p < 0.0) return HP_ERR_DOMAIN;
        sum += p;
        mean += (double)k * p;
        second += (double)k * (double)(k ? k - 1 : 0) * p;
    }
    if (fabs(sum - 1.0) > 1e-6) return HP_ERR_DOMAIN;
    double r = mean > 0.0 ? transmissibility * second / mean : 0.0;
    const double T = transmissibility;
    double S = 0.0;
    if (r > 1.0 && mean > 0.0 && T > 0.0) {
        /* Solve the branching-process extinction equation
         *   q = sum_k (k p_k / mean) (1 - T + T q)^(k-1)
         * for its smallest root, then S = 1 - sum_k p_k (1 - T + T q)^k.
         * This is the tree approximation used by the configuration model;
         * the giant-component theorem is cited in docs/mathematics.md. */
        double q = 0.0;
        double hi = 1.0 - 1e-9;
        double ghi = 0.0;
        for (size_t k = 1; k <= k_max; ++k)
            ghi += ((double)k * pk[k] / mean) *
                   pow(1.0 - T + T * hi, (double)(k - 1));
        ghi -= hi;
        if (ghi < 0.0) {
            double lo = 0.0;
            for (int it = 0; it < 200; ++it) {
                double mid = 0.5 * (lo + hi);
                double h = 0.0;
                for (size_t k = 1; k <= k_max; ++k)
                    h += ((double)k * pk[k] / mean) *
                         pow(1.0 - T + T * mid, (double)(k - 1));
                if (h - mid < 0.0)
                    hi = mid;
                else
                    lo = mid;
            }
            q = 0.5 * (lo + hi);
        }
        double z = 1.0 - T + T * q;
        S = 1.0;
        for (size_t k = 0; k <= k_max; ++k)
            S -= pk[k] * pow(z, (double)k);
        if (S < 0.0) S = 0.0;
        if (S > 1.0) S = 1.0;
    }
    if (giant) *giant = S;
    if (r0) *r0 = r;
    if (mean_degree) *mean_degree = mean;
    return HP_OK;
}

hp_graph *hp_perc_er_generate(size_t n, double p, hp_rng *rng) {
    if (n == 0 || !rng || !(p >= 0.0 && p <= 1.0)) return NULL;
    if (p >= 1.0 && n > 4000) return NULL; /* would be O(n^2) memory */
    size_t cap = 1024, m = 0;
    uint32_t *u = (uint32_t *)hp_xmalloc(cap * sizeof(uint32_t));
    uint32_t *v = (uint32_t *)hp_xmalloc(cap * sizeof(uint32_t));
    const double logq = p < 1.0 ? log1p(-p) : -INFINITY;
    for (size_t a = 0; a + 1 < n; ++a) {
        size_t b = a + 1;
        while (b < n) {
            if (p <= 0.0) break;
            size_t gap = 0;
            if (p < 1.0) {
                double q = hp_rng_uniform(rng);
                while (q <= 0.0) q = hp_rng_uniform(rng);
                gap = (size_t)(log(q) / logq);
            }
            if (gap > n) break;
            b += gap;
            if (b >= n) break;
            if (m == cap) {
                cap *= 2;
                u = (uint32_t *)hp_xrealloc(u, cap * sizeof(uint32_t));
                v = (uint32_t *)hp_xrealloc(v, cap * sizeof(uint32_t));
            }
            u[m] = (uint32_t)a;
            v[m] = (uint32_t)b;
            ++m;
            ++b;
        }
    }
    hp_graph *g = hp_graph_from_edges(n, u, v, NULL, m, 1);
    hp_free(u);
    hp_free(v);
    return g;
}

double hp_perc_bond_threshold(double mean_degree) {
    return mean_degree > 0.0 ? 1.0 / mean_degree : NAN;
}
