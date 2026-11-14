#include "hp/market/backtest.h"

#include <math.h>

#include "hp/core/memory.h"
#include "hp/core/stats.h"

void hp_bt_opts_default(hp_bt_opts *opts) {
    if (!opts) return;
    opts->cost_bps = 5.0;
    opts->apply_costs = 1;
}

void hp_bt_result_free(hp_bt_result *r) {
    if (!r) return;
    hp_free(r->equity);
    hp_free(r->strategy_ret);
    hp_free(r->position);
    r->equity = NULL;
    r->strategy_ret = NULL;
    r->position = NULL;
}

hp_status hp_backtest(const double *asset_ret, size_t T, const double *signal,
                      const hp_bt_opts *opts, hp_bt_result *result,
                      hp_bt_metrics *metrics) {
    if (!asset_ret || !signal || !result || T < 2) return HP_ERR_INVALID;
    hp_bt_opts def;
    if (!opts) {
        hp_bt_opts_default(&def);
        opts = &def;
    }
    if (!(opts->cost_bps >= 0.0)) return HP_ERR_INVALID;
    for (size_t t = 0; t < T; ++t)
        if (!isfinite(signal[t])) return HP_ERR_DOMAIN;
    for (size_t t = 1; t < T; ++t)
        if (!isfinite(asset_ret[t])) return HP_ERR_DOMAIN;

    result->equity = hp_doubles(T);
    result->strategy_ret = hp_doubles(T);
    result->position = hp_doubles(T);
    result->equity[0] = 1.0;
    result->strategy_ret[0] = 0.0;
    result->position[0] = hp_clamp(signal[0], -1.0, 1.0);
    double peak = 1.0, max_dd = 0.0;
    double cost_rate = opts->apply_costs ? opts->cost_bps * 1e-4 : 0.0;

    for (size_t t = 1; t < T; ++t) {
        double pos = hp_clamp(signal[t - 1], -1.0, 1.0);
        double prev = t >= 2 ? hp_clamp(signal[t - 2], -1.0, 1.0) : 0.0;
        double turnover = fabs(pos - prev);
        double r = pos * asset_ret[t] - turnover * cost_rate;
        result->strategy_ret[t] = r;
        result->equity[t] = result->equity[t - 1] * (1.0 + r);
        result->position[t] = hp_clamp(signal[t], -1.0, 1.0);
        if (result->equity[t] > peak) peak = result->equity[t];
        double dd = result->equity[t] / peak - 1.0;
        if (dd < max_dd) max_dd = dd;
    }

    if (metrics) {
        const size_t n = T - 1;
        double *sr = result->strategy_ret + 1;
        double mean = hp_mean(sr, n);
        double sd = hp_sd(sr, n);
        double total = result->equity[T - 1] - 1.0;
        double cagr = result->equity[T - 1] > 0.0
                          ? pow(result->equity[T - 1], 252.0 / (double)n) - 1.0
                          : -1.0;
        metrics->total_return = total;
        metrics->cagr = cagr;
        metrics->ann_vol = sd * sqrt(252.0);
        metrics->sharpe = (sd > 0.0) ? mean / sd * sqrt(252.0) : NAN;
        metrics->max_drawdown = max_dd;
        metrics->calmar = max_dd < 0.0 ? cagr / fabs(max_dd) : NAN;
        double turn = 0.0, expo = 0.0;
        for (size_t t = 1; t < T; ++t) {
            double pos = result->position[t - 1];
            double prev = t >= 2 ? result->position[t - 2] : 0.0;
            turn += fabs(pos - prev);
            expo += fabs(pos);
        }
        metrics->avg_turnover = turn / (double)n;
        metrics->exposure = expo / (double)n;
        size_t hits = 0;
        for (size_t t = 1; t < T; ++t)
            if (result->strategy_ret[t] > 0.0) ++hits;
        metrics->hit_rate = (double)hits / (double)n;
        metrics->final_equity = result->equity[T - 1];
        metrics->n_days = n;
    }
    return HP_OK;
}
