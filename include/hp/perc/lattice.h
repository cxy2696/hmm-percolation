/* HP · two-dimensional site and bond percolation on the square lattice.
 *
 * Site percolation: each vertex open with probability p, edges always present
 * between open neighbours.  Bond percolation: every vertex is present, each of
 * the 2*L*L bonds open with probability p.  Open boundary is the default; with
 * periodic=1 the occupied graph is wrapped on a torus.
 *
 * The order parameter is S = largest cluster / L^2 and the susceptibility is
 * chi = sum_{clusters except the largest} s^2 / L^2, following the standard
 * finite-size definitions.
 */
#ifndef HP_PERC_LATTICE_H
#define HP_PERC_LATTICE_H

#include <stddef.h>
#include <stdint.h>

#include "hp/core/rng.h"
#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    HP_PERC_SITE = 0,
    HP_PERC_BOND = 1
} hp_perc_kind;

typedef struct {
    int L;
    int periodic;
    hp_perc_kind kind;
    uint8_t *site; /* L*L */
    uint8_t *h;    /* L*L: (i,j)-(i,j+1) */
    uint8_t *v;    /* L*L: (i,j)-(i+1,j) */
} hp_lattice;

hp_lattice *hp_lattice_new(int L, hp_perc_kind kind, int periodic);
void hp_lattice_free(hp_lattice *lat);
void hp_lattice_clear(hp_lattice *lat);
/* Fills occupancy from independent Bernoulli(p) draws. */
hp_status hp_lattice_bernoulli(hp_lattice *lat, double p, hp_rng *rng);
/* Site occupancy for site percolation, always 1 for bond percolation. */
int hp_lattice_site_open(const hp_lattice *lat, int i, int j);

typedef struct {
    size_t n_sites;
    size_t n_occupied;
    size_t n_clusters;
    size_t max_cluster;
    double s_max;          /* max_cluster / n_sites */
    double susceptibility; /* non-largest cluster mass / n_sites */
    double occupancy;
    int spanning_x;
    int spanning_y;
    int spanning;
} hp_lattice_stats;

hp_status hp_lattice_analyze(const hp_lattice *lat, hp_lattice_stats *out);

/* Cluster sizes in descending order; *n_out receives the count.  Caller frees
 * the returned array with hp_free(). */
size_t *hp_lattice_cluster_sizes(const hp_lattice *lat, size_t *n_out);

const char *hp_perc_kind_name(hp_perc_kind kind);

/* ---- finite-size scaling helpers ---- */

typedef struct {
    int L;
    hp_perc_kind kind;
    int periodic;
    size_t reps;
    uint64_t seed;
} hp_perc_exp;

/* Fraction of reps with a cluster touching both opposite boundaries. */
double hp_perc_spanning_prob(const hp_perc_exp *exp, double p);

/* Bisection for the p where the spanning probability is 0.5.  The same seed
 * is reused, so the result is deterministic. */
hp_status hp_perc_crossing(const hp_perc_exp *exp, double p_lo, double p_hi,
                           unsigned iters, double *p_cross,
                           double *span_prob);

#ifdef __cplusplus
}
#endif

#endif /* HP_PERC_LATTICE_H */
