#include "hp/perc/corr.h"

#include <math.h>
#include <string.h>

#include "hp/core/memory.h"

hp_graph *hp_graph_from_corr(const double *corr, size_t n, double theta) {
    if (!corr || n == 0 || !isfinite(theta)) return NULL;
    size_t cap = n > 1 ? n * (n - 1) / 2 : 1;
    uint32_t *u = (uint32_t *)hp_xmalloc(cap * sizeof(uint32_t));
    uint32_t *v = (uint32_t *)hp_xmalloc(cap * sizeof(uint32_t));
    size_t m = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double r = corr[i * n + j];
            if (!isfinite(r)) continue;
            if (fabs(r) >= theta) {
                u[m] = (uint32_t)i;
                v[m] = (uint32_t)j;
                ++m;
            }
        }
    }
    hp_graph *g = hp_graph_from_edges(n, u, v, NULL, m, 1);
    hp_free(u);
    hp_free(v);
    return g;
}

hp_status hp_corr_sweep(const double *corr, size_t n, const double *thetas,
                        size_t n_theta, hp_perc_point *out) {
    if (!corr || !thetas || !out || n < 2 || n_theta == 0)
        return HP_ERR_INVALID;
    for (size_t k = 0; k < n_theta; ++k) {
        if (!isfinite(thetas[k])) return HP_ERR_INVALID;
        hp_graph *g = hp_graph_from_corr(corr, n, thetas[k]);
        if (!g) return HP_ERR_NOMEM;
        hp_components c;
        hp_status st = hp_graph_components(g, -INFINITY, &c);
        if (st != HP_OK) {
            hp_graph_free(g);
            return st;
        }
        out[k].theta = thetas[k];
        out[k].n_edges = hp_graph_edge_count(g);
        out[k].n_components = c.n_components;
        out[k].giant_fraction = c.giant_fraction;
        out[k].susceptibility = c.susceptibility;
        out[k].mean_degree = 2.0 * (double)out[k].n_edges / (double)n;
        hp_components_free(&c);
        hp_graph_free(g);
    }
    return HP_OK;
}

hp_status hp_corr_critical(const double *corr, size_t n, double theta_lo,
                           double theta_hi, size_t steps,
                           hp_corr_critical_result *out) {
    if (!corr || !out || n < 2 || steps < 2 || !(theta_lo < theta_hi))
        return HP_ERR_INVALID;
    out->theta_c = theta_lo;
    out->theta_half = theta_lo;
    out->chi_max = 0.0;
    out->giant_at_c = 0.0;
    out->censored_high = 0;
    out->censored_low = 0;
    double step = (theta_hi - theta_lo) / (double)(steps - 1);
    hp_perc_point *pts = (hp_perc_point *)hp_xcalloc(steps, sizeof(*pts));
    double *thetas = hp_doubles(steps);
    for (size_t k = 0; k < steps; ++k) thetas[k] = theta_lo + step * (double)k;
    hp_status st = hp_corr_sweep(corr, n, thetas, steps, pts);
    if (st != HP_OK) {
        hp_free(pts);
        hp_free(thetas);
        return st;
    }
    size_t kbest = 0;
    /* Use the last maximum: on a plateau chi cannot locate the transition,
     * and the last maximum sits at the point where the giant is about to
     * disappear. */
    for (size_t k = 1; k < steps; ++k)
        if (pts[k].susceptibility >= pts[kbest].susceptibility) kbest = k;
    out->chi_max = pts[kbest].susceptibility;
    out->giant_at_c = pts[kbest].giant_fraction;
    double theta_c = pts[kbest].theta;
    if (kbest > 0 && kbest + 1 < steps) {
        double y0 = pts[kbest - 1].susceptibility;
        double y1 = pts[kbest].susceptibility;
        double y2 = pts[kbest + 1].susceptibility;
        double denom = y0 - 2.0 * y1 + y2;
        if (denom < 0.0) {
            double off = 0.5 * (y0 - y2) / denom;
            if (off > 1.0) off = 1.0;
            if (off < -1.0) off = -1.0;
            theta_c += off * step;
        }
    }
    if (kbest == 0) out->censored_low = 1;
    if (kbest + 1 == steps) out->censored_high = 1;
    out->theta_c = theta_c;

    for (size_t k = 1; k < steps; ++k) {
        double g0 = pts[k - 1].giant_fraction;
        double g1 = pts[k].giant_fraction;
        if (g0 >= 0.5 && g1 < 0.5) {
            double frac = (g1 - g0) != 0.0 ? (0.5 - g0) / (g1 - g0) : 0.5;
            out->theta_half = pts[k - 1].theta + frac * step;
            hp_free(pts);
            hp_free(thetas);
            return HP_OK;
        }
    }
    if (pts[0].giant_fraction >= 0.5)
        out->theta_half = theta_hi;
    else
        out->theta_half = theta_lo;
    hp_free(pts);
    hp_free(thetas);
    return HP_OK;
}
