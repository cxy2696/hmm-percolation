/* HP · price panel: dates x symbols, close prices and log returns.
 *
 * CSV layout (first row is the header):
 *     date,A1,A2,...
 *     2012-01-02,100.0,100.0,...
 * Empty/NA cells become NAN; a return is only defined when both closes are
 * finite and positive.
 */
#ifndef HP_MARKET_PANEL_H
#define HP_MARKET_PANEL_H

#include <stddef.h>

#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t T;
    size_t N;
    char **dates;    /* T */
    char **symbols;  /* N */
    double *close;   /* T*N, row-major */
    double *ret;     /* T*N, row-major, ret[0]=NAN */
} hp_panel;

hp_panel *hp_panel_new(size_t T, size_t N);
void hp_panel_free(hp_panel *p);
hp_panel *hp_panel_load_csv(const char *path, hp_status *status);
hp_status hp_panel_save_csv(const hp_panel *p, const char *path);
hp_status hp_panel_write_returns_csv(const hp_panel *p, const char *path);
hp_status hp_panel_compute_returns(hp_panel *p);
int hp_panel_symbol_index(const hp_panel *p, const char *symbol);
/* Copies column j of the return matrix into out[T]. */
hp_status hp_panel_return_series(const hp_panel *p, size_t j, double *out);
/* Cross-sectional equal-weight return for each date (NaN if no asset). */
hp_status hp_panel_ew_returns(const hp_panel *p, double *out);

#ifdef __cplusplus
}
#endif

#endif /* HP_MARKET_PANEL_H */
