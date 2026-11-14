/* HP · undirected CSR graph and component statistics. */
#ifndef HP_PERC_GRAPH_H
#define HP_PERC_GRAPH_H

#include <stddef.h>
#include <stdint.h>

#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Undirected adjacency in CSR form.  Every undirected edge appears twice
 * (u->v and v->u); m is the number of adjacency entries, not edges. */
typedef struct {
    size_t n;
    size_t m;
    size_t *offset;   /* n+1 */
    uint32_t *dst;    /* m */
    double *weight;   /* m, or NULL when unweighted */
} hp_graph;

hp_graph *hp_graph_new(size_t n, size_t m);
void hp_graph_free(hp_graph *g);
/* Builds a graph from an edge list.  Self loops are skipped.  `weight` may be
 * NULL.  Duplicate edges are kept; they do not change connectivity. */
hp_graph *hp_graph_from_edges(size_t n, const uint32_t *u, const uint32_t *v,
                              const double *weight, size_t n_edges,
                              int undirected);
size_t hp_graph_edge_count(const hp_graph *g);

typedef struct {
    size_t n;             /* number of nodes */
    size_t n_components;
    size_t giant;         /* largest component size */
    double giant_fraction;
    double susceptibility; /* sum of s^2 over non-largest components / n */
    size_t *labels;       /* n, component id 0..n_components-1 */
    size_t *sizes;        /* n_components, descending */
} hp_components;

/* Components using only edges with weight >= threshold (all edges when the
 * graph has no weights).  `threshold` is ignored for unweighted graphs. */
hp_status hp_graph_components(const hp_graph *g, double threshold,
                              hp_components *out);
void hp_components_free(hp_components *c);

/* Convenience: giant fraction only. */
hp_status hp_graph_giant(const hp_graph *g, double threshold,
                         double *giant_fraction);

#ifdef __cplusplus
}
#endif

#endif /* HP_PERC_GRAPH_H */
