/* HP · correlation-threshold graphs and their percolation sweep.
 *
 * For N variables we build a graph whose edge (i,j) exists when
 * |corr(i,j)| >= theta.  Lowering theta adds edges, so the giant component
 * appears at a critical threshold theta_c.  Two readings are produced:
 *
 *   theta_c    argmax of the percolation susceptibility chi(theta), which is
 *              the finite-size estimate of the transition;
 *   theta_half threshold where the giant fraction crosses 0.5.
 *
 * A high theta_c means even strongly correlated pairs keep the graph
 * connected.  This is a structural diagnostic, not a forecast of direction.
 */
#ifndef HP_PERC_CORR_H
#define HP_PERC_CORR_H

#include <stddef.h>

#include "hp/core/status.h"
#include "hp/perc/graph.h"

#ifdef __cplusplus
extern "C" {
#endif

hp_graph *hp_graph_from_corr(const double *corr, size_t n, double theta);

typedef struct {
    double theta;
    size_t n_edges;
    size_t n_components;
    double giant_fraction;
    double susceptibility;
    double mean_degree;
} hp_perc_point;

hp_status hp_corr_sweep(const double *corr, size_t n, const double *thetas,
                        size_t n_theta, hp_perc_point *out);

typedef struct {
    double theta_c;
    double theta_half;
    double chi_max;
    double giant_at_c;
    int censored_high;
    int censored_low;
} hp_corr_critical_result;

/* Scans [theta_lo, theta_hi] on a uniform grid.  theta_c is the susceptibility
 * peak (parabolically interpolated); theta_half is the 0.5 giant crossing
 * (linearly interpolated). */
hp_status hp_corr_critical(const double *corr, size_t n, double theta_lo,
                           double theta_hi, size_t steps,
                           hp_corr_critical_result *out);

#ifdef __cplusplus
}
#endif

#endif /* HP_PERC_CORR_H */
