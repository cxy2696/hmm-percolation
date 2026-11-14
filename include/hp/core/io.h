/* HP · file helpers. */
#ifndef HP_CORE_IO_H
#define HP_CORE_IO_H

#include <stddef.h>
#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Reads the complete file into a NUL-terminated malloc'ed buffer. */
hp_status hp_read_file(const char *path, char **out, size_t *out_len);
hp_status hp_write_file(const char *path, const char *data, size_t len);
hp_status hp_write_text(const char *path, const char *text);

/* strdup that aborts on OOM. */
char *hp_strdup(const char *s);
/* Trims leading/trailing ASCII whitespace in place, returns s. */
char *hp_trim(char *s);
/* Case-insensitive ASCII compare. */
int hp_str_ieq(const char *a, const char *b);

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_IO_H */
