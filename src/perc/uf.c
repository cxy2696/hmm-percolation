#include "hp/perc/uf.h"

#include "hp/core/memory.h"

hp_uf *hp_uf_new(size_t n) {
    if (n == 0) return NULL;
    hp_uf *uf = (hp_uf *)hp_xcalloc(1, sizeof(*uf));
    uf->n = n;
    uf->count = n;
    uf->parent = hp_sizes(n);
    uf->size = hp_sizes(n);
    for (size_t i = 0; i < n; ++i) {
        uf->parent[i] = i;
        uf->size[i] = 1;
    }
    return uf;
}

void hp_uf_free(hp_uf *uf) {
    if (!uf) return;
    hp_free(uf->parent);
    hp_free(uf->size);
    hp_free(uf);
}

void hp_uf_reset(hp_uf *uf) {
    if (!uf) return;
    uf->count = uf->n;
    for (size_t i = 0; i < uf->n; ++i) {
        uf->parent[i] = i;
        uf->size[i] = 1;
    }
}

size_t hp_uf_find(hp_uf *uf, size_t x) {
    size_t root = x;
    while (uf->parent[root] != root) root = uf->parent[root];
    while (uf->parent[x] != root) {
        size_t next = uf->parent[x];
        uf->parent[x] = root;
        x = next;
    }
    return root;
}

int hp_uf_union(hp_uf *uf, size_t a, size_t b) {
    size_t ra = hp_uf_find(uf, a);
    size_t rb = hp_uf_find(uf, b);
    if (ra == rb) return 0;
    if (uf->size[ra] < uf->size[rb]) {
        size_t t = ra;
        ra = rb;
        rb = t;
    }
    uf->parent[rb] = ra;
    uf->size[ra] += uf->size[rb];
    --uf->count;
    return 1;
}

size_t hp_uf_size(hp_uf *uf, size_t x) { return uf->size[hp_uf_find(uf, x)]; }
