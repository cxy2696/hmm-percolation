#include "hp/cross/state_graph.h"

#include <math.h>
#include <string.h>

#include "hp/core/memory.h"

hp_graph *hp_graph_from_transition(const double *a, size_t n, double theta) {
    if (!a || n < 2 || !isfinite(theta) || theta < 0.0) return NULL;
    size_t cap = n * (n - 1) / 2;
    uint32_t *u = (uint32_t *)hp_xmalloc((cap ? cap : 1) * sizeof(uint32_t));
    uint32_t *v = (uint32_t *)hp_xmalloc((cap ? cap : 1) * sizeof(uint32_t));
    size_t m = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double aij = a[i * n + j];
            double aji = a[j * n + i];
            if (!isfinite(aij) || !isfinite(aji)) continue;
            if (aij >= theta || aji >= theta) {
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

hp_status hp_transition_components(const double *a, size_t n, double theta,
                                   hp_components *out) {
    if (!a || !out || n < 2) return HP_ERR_INVALID;
    hp_graph *g = hp_graph_from_transition(a, n, theta);
    if (!g) return HP_ERR_INVALID;
    hp_status st = hp_graph_components(g, -INFINITY, out);
    hp_graph_free(g);
    return st;
}
