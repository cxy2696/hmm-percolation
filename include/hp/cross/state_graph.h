/* HP · the HMM transition matrix as a sparse graph.
 *
 * The link between the two theories in the source note: if a_ij is read as an
 * edge weight, thresholding A yields a random-graph-like object whose giant
 * component tells us whether the state space is effectively communicating.
 * A well-sized model has one giant component; an oversized one fragments into
 * weakly linked cliques (and states that are rarely visited).
 */
#ifndef HP_CROSS_STATE_GRAPH_H
#define HP_CROSS_STATE_GRAPH_H

#include <stddef.h>

#include "hp/core/status.h"
#include "hp/perc/graph.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Undirected projection: edge (i,j) when a_ij >= theta or a_ji >= theta. */
hp_graph *hp_graph_from_transition(const double *a, size_t n, double theta);
hp_status hp_transition_components(const double *a, size_t n, double theta,
                                   hp_components *out);

#ifdef __cplusplus
}
#endif

#endif /* HP_CROSS_STATE_GRAPH_H */
