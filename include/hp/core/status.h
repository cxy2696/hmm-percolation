/* HP · error model.
 *
 * Every public function that can fail returns hp_status (or NULL for a
 * constructor).  Status codes are stable and are part of the ABI: add new
 * codes at the end, never renumber.  HP_OK is 0 so `if (st)` is wrong on
 * purpose -- always compare against HP_OK.
 */
#ifndef HP_CORE_STATUS_H
#define HP_CORE_STATUS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    HP_OK = 0,
    HP_ERR_INVALID = -1,   /* caller passed a malformed argument        */
    HP_ERR_NOMEM = -2,     /* allocation failed (only from raw allocs)  */
    HP_ERR_IO = -3,        /* open/read/write failure                   */
    HP_ERR_FORMAT = -4,    /* file content is not the expected format   */
    HP_ERR_DOMAIN = -5,    /* probability zero, non-finite value, etc.  */
    HP_ERR_CONVERGE = -6,  /* iterative algorithm hit max_iter          */
    HP_ERR_UNSUPPORTED = -7,
    HP_ERR_RANGE = -8,     /* index out of range                        */
    HP_ERR_STATE = -9,     /* object in the wrong state                 */
    HP_ERR_LAST = -10
} hp_status;

/* Short machine-readable name, e.g. "HP_ERR_INVALID". */
const char *hp_status_name(hp_status st);
/* One-line human readable description (English). */
const char *hp_status_message(hp_status st);

/* Run `expr`; on failure return its status from the enclosing function. */
#define HP_TRY(expr)            \
    do {                        \
        hp_status hp__st = (expr); \
        if (hp__st != HP_OK) return hp__st; \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif /* HP_CORE_STATUS_H */
