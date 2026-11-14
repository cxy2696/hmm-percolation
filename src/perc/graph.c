#include "hp/perc/graph.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/memory.h"
#include "hp/perc/uf.h"

hp_graph *hp_graph_new(size_t n, size_t m) {
    if (n == 0) return NULL;
    hp_graph *g = (hp_graph *)hp_xcalloc(1, sizeof(*g));
    g->n = n;
    g->m = m;
    g->offset = hp_sizes(n + 1);
    g->dst = m ? hp_u32s(m) : NULL;
    return g;
}

void hp_graph_free(hp_graph *g) {
    if (!g) return;
    hp_free(g->offset);
    hp_free(g->dst);
    hp_free(g->weight);
    hp_free(g);
}

hp_graph *hp_graph_from_edges(size_t n, const uint32_t *u, const uint32_t *v,
                              const double *weight, size_t n_edges,
                              int undirected) {
    if (n == 0) return NULL;
    if (n_edges && (!u || !v)) return NULL;
    size_t valid = 0;
    for (size_t e = 0; e < n_edges; ++e) {
        if (u[e] >= n || v[e] >= n || u[e] == v[e]) continue;
        ++valid;
    }
    size_t entries = undirected ? valid * 2 : valid;
    hp_graph *g = hp_graph_new(n, entries);
    g->weight = weight ? hp_doubles(entries) : NULL;
    size_t *deg = hp_sizes(n);
    for (size_t e = 0; e < n_edges; ++e) {
        if (u[e] >= n || v[e] >= n || u[e] == v[e]) continue;
        ++deg[u[e]];
        if (undirected) ++deg[v[e]];
    }
    for (size_t i = 0; i < n; ++i) g->offset[i + 1] = g->offset[i] + deg[i];
    size_t *cursor = hp_sizes(n);
    memcpy(cursor, g->offset, n * sizeof(size_t));
    for (size_t e = 0; e < n_edges; ++e) {
        if (u[e] >= n || v[e] >= n || u[e] == v[e]) continue;
        size_t pu = cursor[u[e]]++;
        g->dst[pu] = v[e];
        if (g->weight) g->weight[pu] = weight ? weight[e] : 1.0;
        if (undirected) {
            size_t pv = cursor[v[e]]++;
            g->dst[pv] = u[e];
            if (g->weight) g->weight[pv] = weight ? weight[e] : 1.0;
        }
    }
    hp_free(deg);
    hp_free(cursor);
    return g;
}

size_t hp_graph_edge_count(const hp_graph *g) { return g ? g->m / 2 : 0; }

typedef struct {
    size_t root;
    size_t size;
} comp_pair;

static int cmp_comp(const void *pa, const void *pb) {
    const comp_pair *a = (const comp_pair *)pa;
    const comp_pair *b = (const comp_pair *)pb;
    if (a->size != b->size) return a->size > b->size ? -1 : 1;
    return a->root < b->root ? -1 : (a->root > b->root);
}

hp_status hp_graph_components(const hp_graph *g, double threshold,
                              hp_components *out) {
    if (!g || !out || g->n == 0) return HP_ERR_INVALID;
    memset(out, 0, sizeof(*out));
    hp_uf *uf = hp_uf_new(g->n);
    for (size_t u = 0; u < g->n; ++u) {
        for (size_t p = g->offset[u]; p < g->offset[u + 1]; ++p) {
            if (g->weight && !(g->weight[p] >= threshold)) continue;
            hp_uf_union(uf, u, g->dst[p]);
        }
    }
    size_t *cnt = hp_sizes(g->n);
    for (size_t i = 0; i < g->n; ++i) cnt[hp_uf_find(uf, i)]++;
    comp_pair *list = (comp_pair *)hp_xmalloc(g->n * sizeof(comp_pair));
    size_t nc = 0;
    for (size_t i = 0; i < g->n; ++i) {
        if (cnt[i] > 0) {
            list[nc].root = i;
            list[nc].size = cnt[i];
            ++nc;
        }
    }
    qsort(list, nc, sizeof(comp_pair), cmp_comp);

    size_t *root_id = hp_sizes(g->n);
    for (size_t i = 0; i < g->n; ++i) root_id[i] = SIZE_MAX;
    out->labels = hp_sizes(g->n);
    out->sizes = hp_sizes(nc ? nc : 1);
    for (size_t k = 0; k < nc; ++k) {
        root_id[list[k].root] = k;
        out->sizes[k] = list[k].size;
    }
    for (size_t i = 0; i < g->n; ++i)
        out->labels[i] = root_id[hp_uf_find(uf, i)];
    out->n = g->n;
    out->n_components = nc;
    out->giant = nc ? list[0].size : 0;
    out->giant_fraction = (double)out->giant / (double)g->n;
    double chi = 0.0;
    for (size_t k = 1; k < nc; ++k)
        chi += (double)list[k].size * (double)list[k].size;
    out->susceptibility = chi / (double)g->n;

    hp_free(cnt);
    hp_free(list);
    hp_free(root_id);
    hp_uf_free(uf);
    return HP_OK;
}

void hp_components_free(hp_components *c) {
    if (!c) return;
    hp_free(c->labels);
    hp_free(c->sizes);
    memset(c, 0, sizeof(*c));
}

hp_status hp_graph_giant(const hp_graph *g, double threshold,
                         double *giant_fraction) {
    if (!giant_fraction) return HP_ERR_INVALID;
    hp_components c;
    HP_TRY(hp_graph_components(g, threshold, &c));
    *giant_fraction = c.giant_fraction;
    hp_components_free(&c);
    return HP_OK;
}
