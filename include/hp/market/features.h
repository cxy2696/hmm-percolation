/* HP · market features for the regime HMM.
 *
 * Two columns: daily log return in percent, and log annualised rolling
 * volatility in percent.  The volatility column is a slow state variable that
 * separates calm from stressed regimes; the return column carries direction.
 * Rows before `start` are NAN and the app trims them.
 */
#ifndef HP_MARKET_FEATURES_H
#define HP_MARKET_FEATURES_H

#include <stddef.h>

#include "hp/core/status.h"
#include "hp/market/panel.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t T;
    size_t dim;
    size_t start;
    double *ret; /* T */
    double *vol; /* T, annualised */
    double *x;   /* T*dim */
} hp_features;

hp_features *hp_features_from_returns(const double *ret, size_t T,
                                      size_t vol_window);
hp_features *hp_features_from_panel(const hp_panel *p, size_t asset_idx,
                                    size_t vol_window);
void hp_features_free(hp_features *f);

#ifdef __cplusplus
}
#endif

#endif /* HP_MARKET_FEATURES_H */
