/* HP · event-driven backtest.
 *
 * Timing (no look-ahead): signal[t] is decided after observing the close of t
 * and is applied to the return from t to t+1.  A change of position is charged
 * cost_bps on the traded notional, including the initial entry.
 */
#ifndef HP_MARKET_BACKTEST_H
#define HP_MARKET_BACKTEST_H

#include <stddef.h>

#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double cost_bps;
    int apply_costs;
    /* Periods per year used for CAGR/vol/Sharpe.  Default 252 (daily). */
    double periods_per_year;
} hp_bt_opts;

void hp_bt_opts_default(hp_bt_opts *opts);

typedef struct {
    double total_return;
    double cagr;
    double ann_vol;
    double sharpe;
    double max_drawdown;
    double calmar;
    double avg_turnover;
    double exposure;
    double hit_rate;
    double final_equity;
    size_t n_days;
} hp_bt_metrics;

typedef struct {
    double *equity;       /* T */
    double *strategy_ret; /* T */
    double *position;     /* T, same as signal */
} hp_bt_result;

void hp_bt_result_free(hp_bt_result *r);

hp_status hp_backtest(const double *asset_ret, size_t T, const double *signal,
                      const hp_bt_opts *opts, hp_bt_result *result,
                      hp_bt_metrics *metrics);

#ifdef __cplusplus
}
#endif

#endif /* HP_MARKET_BACKTEST_H */
