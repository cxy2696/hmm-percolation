#include "hp/core/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void *hp_malloc(size_t n) { return malloc(n ? n : 1); }

void *hp_calloc(size_t n, size_t size) {
    if (n && size > SIZE_MAX / n) return NULL;
    return calloc(n ? n : 1, size ? size : 1);
}

void *hp_realloc(void *p, size_t n) { return realloc(p, n ? n : 1); }

void hp_free(void *p) { free(p); }

static void hp_oom(size_t n) {
    fprintf(stderr, "hp: out of memory allocating %zu bytes\n", n);
    abort();
}

void *hp_xmalloc(size_t n) {
    void *p = hp_malloc(n);
    if (!p) hp_oom(n);
    return p;
}

void *hp_xcalloc(size_t n, size_t size) {
    void *p = hp_calloc(n, size);
    if (!p) hp_oom(n * (size ? size : 1));
    return p;
}

void *hp_xrealloc(void *p, size_t n) {
    void *q = hp_realloc(p, n);
    if (!q) hp_oom(n);
    return q;
}

double *hp_doubles(size_t n) { return (double *)hp_xcalloc(n ? n : 1, sizeof(double)); }
size_t *hp_sizes(size_t n) { return (size_t *)hp_xcalloc(n ? n : 1, sizeof(size_t)); }
uint32_t *hp_u32s(size_t n) { return (uint32_t *)hp_xcalloc(n ? n : 1, sizeof(uint32_t)); }
uint8_t *hp_u8s(size_t n) { return (uint8_t *)hp_xcalloc(n ? n : 1, sizeof(uint8_t)); }
int32_t *hp_i32s(size_t n) { return (int32_t *)hp_xcalloc(n ? n : 1, sizeof(int32_t)); }

void hp_mem_zero(void *p, size_t n) { memset(p, 0, n); }

void hp_copy_doubles(double *dst, const double *src, size_t n) {
    if (n) memcpy(dst, src, n * sizeof(double));
}

void hp_fill_doubles(double *dst, double v, size_t n) {
    for (size_t i = 0; i < n; ++i) dst[i] = v;
}
