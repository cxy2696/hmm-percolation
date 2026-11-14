/* HP · mean-field percolation: Erdos-Renyi and configuration models.
 *
 * For a configuration model with degree distribution p_k and bond
 * transmissibility T, the tree approximation solves the neighbour extinction
 * probability q from
 *
 *     q = sum_k (k p_k / <k>) (1 - T + T q)^(k-1),
 *
 * and then
 *
 *     S = 1 - sum_k p_k (1 - T + T q)^k.
 *
 * A non-trivial root exists when the mean excess degree
 *
 *     R0 = T * <k(k-1)> / <k>  >  1.
 *
 * For a Poisson degree distribution q = 1-S and this reduces to
 * S = 1 - exp(-c T S).  The identification of S with a giant component of a
 * finite random graph is a theorem under conditions; see docs/mathematics.md.
 */
#ifndef HP_PERC_RANDOM_GRAPH_H
#define HP_PERC_RANDOM_GRAPH_H

#include <stddef.h>

#include "hp/core/rng.h"
#include "hp/core/status.h"
#include "hp/perc/graph.h"

#ifdef __cplusplus
extern "C" {
#endif

/* S solves S = 1 - exp(-c S); returns 0 for c <= 1. */
double hp_perc_er_giant_theory(double mean_degree);

/* pk[k], k = 0..k_max, must sum to 1.  transmissibility in [0,1]. */
hp_status hp_perc_config_giant(const double *pk, size_t k_max,
                               double transmissibility, double *giant,
                               double *r0, double *mean_degree);

/* Sparse Erdos-Renyi G(n,p) with geometric skip sampling: O(n + m). */
hp_graph *hp_perc_er_generate(size_t n, double p, hp_rng *rng);

/* Bond percolation threshold on a graph with mean degree mu: p_c = 1/mu. */
double hp_perc_bond_threshold(double mean_degree);

#ifdef __cplusplus
}
#endif

#endif /* HP_PERC_RANDOM_GRAPH_H */
