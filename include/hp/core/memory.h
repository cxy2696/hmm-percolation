/* HP · allocation helpers.
 *
 * The numeric code uses the checked `hp_x*` family: on allocation failure it
 * prints a diagnostic and aborts.  Rationale: every algorithm here allocates a
 * handful of arrays whose sizes are known up front; threading NOMEM through
 * every inner loop would add noise without making the process recoverable.
 * The raw `hp_*` family returns NULL and is available to callers that need to
 * handle failure themselves.
 */
#ifndef HP_CORE_MEMORY_H
#define HP_CORE_MEMORY_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void *hp_malloc(size_t n);
void *hp_calloc(size_t n, size_t size);
void *hp_realloc(void *p, size_t n);
void hp_free(void *p);

void *hp_xmalloc(size_t n);
void *hp_xcalloc(size_t n, size_t size);
void *hp_xrealloc(void *p, size_t n);

/* Typed zero-initialised arrays (abort on OOM). */
double *hp_doubles(size_t n);
size_t *hp_sizes(size_t n);
uint32_t *hp_u32s(size_t n);
uint8_t *hp_u8s(size_t n);
int32_t *hp_i32s(size_t n);

void hp_mem_zero(void *p, size_t n);
/* Copy n doubles; src and dst must not overlap. */
void hp_copy_doubles(double *dst, const double *src, size_t n);
/* Fill n doubles with v. */
void hp_fill_doubles(double *dst, double v, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_MEMORY_H */
