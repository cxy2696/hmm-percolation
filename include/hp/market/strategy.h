/* HP · small, causal signal rules.
 *
 * These functions only look backwards.  A signal at index t is computed from
 * data at or before t and is applied to the return at t+1 by hp_backtest().
 * NaN marks "not enough history"; the application layer decides whether to
 * hold cash or wait.
 */
#ifndef HP_MARKET_STRATEGY_H
#define HP_MARKET_STRATEGY_H

#include <stddef.h>

#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t lookback;
    double target_vol;      /* annualised, e.g. 0.10 */
    double min_position;    /* usually 0 */
    double max_position;    /* usually 1 */
    double periods_per_year;
} hp_vol_target_opts;

void hp_vol_target_opts_default(hp_vol_target_opts *opts);

/* pos[t] = clip(target_vol / trailing annualised vol, min, max). */
hp_status hp_signal_vol_target(const double *ret, size_t T,
                               const hp_vol_target_opts *opts, double *pos);

typedef struct {
    size_t lookback;
    double threshold;     /* cumulative log return threshold */
    double long_position;
    double flat_position;
} hp_momentum_opts;

void hp_momentum_opts_default(hp_momentum_opts *opts);

/* Long if the trailing cumulative log return exceeds `threshold`, else flat. */
hp_status hp_signal_momentum(const double *ret, size_t T,
                             const hp_momentum_opts *opts, double *pos);

/* pos[t] = sum_j filtered[t*n+j] * exposure[j]. */
hp_status hp_signal_state_exposure(const double *filtered, size_t T, size_t n,
                                   const double *exposure, double *pos);

/* out[t] = a[t] * b[t]; NaN propagates. */
hp_status hp_signal_multiply(const double *a, const double *b, size_t T,
                             double *out);

hp_status hp_signal_constant(size_t T, double position, double *pos);

#ifdef __cplusplus
}
#endif

#endif /* HP_MARKET_STRATEGY_H */
