/* HP · union-find with union by size and path compression. */
#ifndef HP_PERC_UF_H
#define HP_PERC_UF_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t n;
    size_t count;   /* number of components, including untouched nodes */
    size_t *parent;
    size_t *size;
} hp_uf;

hp_uf *hp_uf_new(size_t n);
void hp_uf_free(hp_uf *uf);
void hp_uf_reset(hp_uf *uf);
size_t hp_uf_find(hp_uf *uf, size_t x);
/* Returns 1 if a and b were in different sets. */
int hp_uf_union(hp_uf *uf, size_t a, size_t b);
size_t hp_uf_size(hp_uf *uf, size_t x);

#ifdef __cplusplus
}
#endif

#endif /* HP_PERC_UF_H */
