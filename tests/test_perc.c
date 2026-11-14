#include <math.h>
#include <string.h>

#include "hp/core/memory.h"
#include "hp/core/rng.h"
#include "hp/perc/corr.h"
#include "hp/perc/graph.h"
#include "hp/perc/lattice.h"
#include "hp/perc/random_graph.h"
#include "hp/perc/uf.h"
#include "hp_test.h"

HP_TEST(test_uf_basics) {
    hp_uf *uf = hp_uf_new(6);
    HP_CHECK(uf != NULL);
    HP_CHECK(hp_uf_union(uf, 0, 1) == 1);
    HP_CHECK(hp_uf_union(uf, 1, 2) == 1);
    HP_CHECK(hp_uf_union(uf, 2, 0) == 0);
    HP_CHECK(hp_uf_union(uf, 3, 4) == 1);
    HP_CHECK(hp_uf_size(uf, 0) == 3);
    HP_CHECK(hp_uf_size(uf, 4) == 2);
    HP_CHECK(uf->count == 3);
    hp_uf_reset(uf);
    HP_CHECK(uf->count == 6);
    HP_CHECK(hp_uf_size(uf, 0) == 1);
    hp_uf_free(uf);
    return 0;
}

HP_TEST(test_graph_components) {
    uint32_t u[5] = {0, 1, 2, 3, 4};
    uint32_t v[5] = {1, 2, 3, 4, 5};
    double w[5] = {0.9, 0.8, 0.2, 0.7, 0.6};
    hp_graph *g = hp_graph_from_edges(6, u, v, w, 5, 1);
    HP_CHECK(g != NULL);
    HP_CHECK(hp_graph_edge_count(g) == 5);
    hp_components c;
    HP_CHECK(hp_graph_components(g, 0.5, &c) == HP_OK);
    HP_CHECK(c.n_components == 2);
    HP_CHECK(c.giant == 3);
    HP_CHECK_NEAR(c.giant_fraction, 0.5, 1e-15);
    HP_CHECK_NEAR(c.susceptibility, 1.5, 1e-15);
    HP_CHECK(c.labels[0] == c.labels[2]);
    HP_CHECK(c.labels[0] != c.labels[3]);
    hp_components_free(&c);

    HP_CHECK(hp_graph_components(g, 0.75, &c) == HP_OK);
    HP_CHECK(c.n_components == 4);
    HP_CHECK(c.giant == 3);
    HP_CHECK_NEAR(c.susceptibility, 0.5, 1e-15);
    hp_components_free(&c);

    HP_CHECK(hp_graph_components(g, 0.95, &c) == HP_OK);
    HP_CHECK(c.n_components == 6);
    HP_CHECK(c.giant == 1);
    HP_CHECK_NEAR(c.susceptibility, 5.0 / 6.0, 1e-15);
    hp_components_free(&c);
    hp_graph_free(g);
    return 0;
}

static void set_site(hp_lattice *lat, int i, int j, uint8_t v) {
    lat->site[(size_t)i * (size_t)lat->L + (size_t)j] = v;
}

HP_TEST(test_lattice_hand) {
    hp_lattice *lat = hp_lattice_new(4, HP_PERC_SITE, 0);
    HP_CHECK(lat != NULL);
    hp_lattice_clear(lat);
    int open[5][2] = {{0, 0}, {0, 1}, {1, 1}, {2, 1}, {2, 2}};
    for (int k = 0; k < 5; ++k) set_site(lat, open[k][0], open[k][1], 1);
    hp_lattice_stats st;
    HP_CHECK(hp_lattice_analyze(lat, &st) == HP_OK);
    HP_CHECK(st.n_occupied == 5);
    HP_CHECK(st.max_cluster == 5);
    HP_CHECK(st.n_clusters == 1);
    HP_CHECK_NEAR(st.s_max, 5.0 / 16.0, 1e-15);
    HP_CHECK_NEAR(st.susceptibility, 0.0, 1e-15);
    HP_CHECK(st.spanning == 0);

    set_site(lat, 3, 0, 1);
    HP_CHECK(hp_lattice_analyze(lat, &st) == HP_OK);
    HP_CHECK(st.n_occupied == 6);
    HP_CHECK(st.n_clusters == 2);
    HP_CHECK(st.max_cluster == 5);
    HP_CHECK_NEAR(st.susceptibility, 1.0 / 16.0, 1e-15);

    size_t ns = 0;
    size_t *sizes = hp_lattice_cluster_sizes(lat, &ns);
    HP_CHECK(ns == 2);
    HP_CHECK(sizes[0] == 5 && sizes[1] == 1);
    hp_free(sizes);
    hp_lattice_free(lat);
    return 0;
}

HP_TEST(test_lattice_pc_sanity) {
    hp_perc_exp e;
    e.L = 128;
    e.kind = HP_PERC_BOND;
    e.periodic = 0;
    e.reps = 160;
    e.seed = 2024;
    double below = hp_perc_spanning_prob(&e, 0.40);
    double above = hp_perc_spanning_prob(&e, 0.55);
    HP_CHECK(below < 0.05);
    HP_CHECK(above > 0.95);

    e.kind = HP_PERC_SITE;
    e.seed = 2025;
    below = hp_perc_spanning_prob(&e, 0.50);
    above = hp_perc_spanning_prob(&e, 0.66);
    HP_CHECK(below < 0.05);
    HP_CHECK(above > 0.85);
    return 0;
}

HP_TEST(test_lattice_crossing) {
    hp_perc_exp e;
    e.L = 64;
    e.kind = HP_PERC_BOND;
    e.periodic = 0;
    e.reps = 80;
    e.seed = 31;
    double pc = 0, sp = 0;
    HP_CHECK(hp_perc_crossing(&e, 0.42, 0.58, 14, &pc, &sp) == HP_OK);
    HP_CHECK(fabs(pc - 0.5) < 0.04);

    e.kind = HP_PERC_SITE;
    e.seed = 32;
    HP_CHECK(hp_perc_crossing(&e, 0.52, 0.68, 14, &pc, &sp) == HP_OK);
    HP_CHECK(fabs(pc - 0.5927) < 0.04);
    return 0;
}

static void poisson_pk(double *pk, size_t kmax, double c) {
    double fact = 1.0, term = exp(-c);
    pk[0] = term;
    for (size_t k = 1; k <= kmax; ++k) {
        fact *= (double)k;
        term = exp(-c) * pow(c, (double)k) / fact;
        pk[k] = term;
    }
}

HP_TEST(test_meanfield_theory) {
    const size_t kmax = 30;
    double pk[31];
    poisson_pk(pk, kmax, 2.0);
    double sum = 0;
    for (size_t k = 0; k <= kmax; ++k) sum += pk[k];
    for (size_t k = 0; k <= kmax; ++k) pk[k] /= sum;

    double S = 0, r0 = 0, mean = 0;
    HP_CHECK(hp_perc_config_giant(pk, kmax, 1.0, &S, &r0, &mean) == HP_OK);
    HP_CHECK_NEAR(mean, 2.0, 1e-3);
    HP_CHECK_NEAR(r0, 2.0, 1e-3);
    HP_CHECK_NEAR(S, hp_perc_er_giant_theory(2.0), 1e-6);
    HP_CHECK_NEAR(S, 0.797, 2e-3);

    HP_CHECK(hp_perc_config_giant(pk, kmax, 0.4, &S, &r0, &mean) == HP_OK);
    HP_CHECK_NEAR(r0, 0.8, 1e-3);
    HP_CHECK(S == 0.0);

    HP_CHECK(hp_perc_config_giant(pk, kmax, 0.6, &S, &r0, &mean) == HP_OK);
    HP_CHECK_NEAR(r0, 1.2, 1e-3);
    HP_CHECK(S > 0.0 && S < 0.5);
    HP_CHECK_NEAR(hp_perc_bond_threshold(2.0), 0.5, 1e-15);

    /* A 2-regular graph is the critical cycle case: P(2)=1 gives R0=1 at
     * T=1, so the tree approximation must return S=0, not S=1. */
    double p2[3] = {0.0, 0.0, 1.0};
    HP_CHECK(hp_perc_config_giant(p2, 2, 1.0, &S, &r0, &mean) == HP_OK);
    HP_CHECK_NEAR(mean, 2.0, 1e-12);
    HP_CHECK_NEAR(r0, 1.0, 1e-12);
    HP_CHECK(S == 0.0);

    /* 3-regular: R0=2, the tree approximation has S=1 (a connected
     * 3-regular graph is a giant). */
    double p3[4] = {0.0, 0.0, 0.0, 1.0};
    HP_CHECK(hp_perc_config_giant(p3, 3, 1.0, &S, &r0, &mean) == HP_OK);
    HP_CHECK_NEAR(r0, 2.0, 1e-12);
    HP_CHECK(S > 0.99);

    /* Half degree-1, half degree-3, T=1.  The old S-equation gives S=1;
     * the neighbour-extinction equation gives q=1/3 and S=22/27. */
    double pmix[4] = {0.0, 0.5, 0.0, 0.5};
    HP_CHECK(hp_perc_config_giant(pmix, 3, 1.0, &S, &r0, &mean) == HP_OK);
    HP_CHECK_NEAR(mean, 2.0, 1e-12);
    HP_CHECK_NEAR(r0, 1.5, 1e-12);
    HP_CHECK_NEAR(S, 22.0 / 27.0, 1e-12);
    return 0;
}

HP_TEST(test_er_simulation) {
    hp_rng rng;
    hp_rng_seed(&rng, 4242);
    hp_graph *g = hp_perc_er_generate(20000, 0.0002, &rng);
    HP_CHECK(g != NULL);
    hp_components c;
    HP_CHECK(hp_graph_components(g, -INFINITY, &c) == HP_OK);
    double theory = hp_perc_er_giant_theory(4.0);
    HP_CHECK_NEAR(c.giant_fraction, theory, 0.02);
    HP_CHECK(c.n_components < 20000);
    hp_components_free(&c);
    hp_graph_free(g);

    hp_graph *g2 = hp_perc_er_generate(20000, 0.00004, &rng);
    HP_CHECK(g2 != NULL);
    HP_CHECK(hp_graph_components(g2, -INFINITY, &c) == HP_OK);
    HP_CHECK(c.giant_fraction < 0.01);
    hp_components_free(&c);
    hp_graph_free(g2);
    return 0;
}

static void block_corr(double *corr, size_t n, size_t half, double within) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j)
                corr[i * n + j] = 1.0;
            else if ((i < half) == (j < half))
                corr[i * n + j] = within;
            else
                corr[i * n + j] = 0.0;
        }
    }
}

HP_TEST(test_corr_graph) {
    const size_t n = 20;
    double *corr = hp_doubles(n * n);
    block_corr(corr, n, 10, 0.8);
    hp_graph *g = hp_graph_from_corr(corr, n, 0.5);
    HP_CHECK(g != NULL);
    HP_CHECK(hp_graph_edge_count(g) == 90);
    hp_components c;
    HP_CHECK(hp_graph_components(g, -INFINITY, &c) == HP_OK);
    HP_CHECK(c.n_components == 2);
    HP_CHECK(c.giant == 10);
    HP_CHECK_NEAR(c.giant_fraction, 0.5, 1e-15);
    hp_components_free(&c);
    hp_graph_free(g);

    double thetas[6] = {0.3, 0.5, 0.75, 0.8, 0.81, 0.9};
    hp_perc_point pts[6];
    HP_CHECK(hp_corr_sweep(corr, n, thetas, 6, pts) == HP_OK);
    HP_CHECK_NEAR(pts[0].giant_fraction, 0.5, 1e-15);
    HP_CHECK_NEAR(pts[1].giant_fraction, 0.5, 1e-15);
    HP_CHECK_NEAR(pts[3].giant_fraction, 0.5, 1e-15);
    HP_CHECK_NEAR(pts[4].giant_fraction, 0.05, 1e-15);
    HP_CHECK(pts[1].mean_degree == 9.0);

    hp_corr_critical_result cr;
    HP_CHECK(hp_corr_critical(corr, n, 0.5, 0.9, 41, &cr) == HP_OK);
    HP_CHECK(cr.theta_half > 0.79 && cr.theta_half < 0.81);
    HP_CHECK(cr.theta_c >= 0.5 && cr.theta_c <= 0.81);
    HP_CHECK(cr.chi_max >= 5.0 - 1e-12);
    hp_free(corr);
    return 0;
}

void hp_test_perc_register(void) {
    HP_REG(test_uf_basics);
    HP_REG(test_graph_components);
    HP_REG(test_lattice_hand);
    HP_REG(test_lattice_pc_sanity);
    HP_REG(test_lattice_crossing);
    HP_REG(test_meanfield_theory);
    HP_REG(test_er_simulation);
    HP_REG(test_corr_graph);
}
