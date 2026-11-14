#include "hp/market/strategy.h"

#include <math.h>

#include "hp/core/stats.h"

void hp_vol_target_opts_default(hp_vol_target_opts *opts) {
    if (!opts) return;
    opts->lookback = 20;
    opts->target_vol = 0.10;
    opts->min_position = 0.0;
    opts->max_position = 1.0;
    opts->periods_per_year = 252.0;
}

hp_status hp_signal_vol_target(const double *ret, size_t T,
                               const hp_vol_target_opts *opts, double *pos) {
    if (!ret || !opts || !pos || T == 0) return HP_ERR_INVALID;
    if (opts->lookback < 2 || !isfinite(opts->target_vol) ||
        !(opts->target_vol > 0.0) || !isfinite(opts->periods_per_year) ||
        !(opts->periods_per_year > 0.0) ||
        !(opts->min_position <= opts->max_position) ||
        !isfinite(opts->min_position) || !isfinite(opts->max_position))
        return HP_ERR_INVALID;
    for (size_t t = 0; t < T; ++t) {
        pos[t] = NAN;
        if (t + 1 < opts->lookback) continue;
        double mean = 0.0, var = 0.0;
        size_t count = 0;
        for (size_t k = t + 1 - opts->lookback; k <= t; ++k) {
            if (!isfinite(ret[k])) {
                count = 0;
                break;
            }
            double d = ret[k] - mean;
            mean += d / (double)(count + 1);
            var += d * (ret[k] - mean);
            ++count;
        }
        if (count != opts->lookback) continue;
        var /= (double)(count - 1);
        double vol = sqrt(var * opts->periods_per_year);
        double p = vol > 1e-12 ? opts->target_vol / vol : opts->max_position;
        pos[t] = hp_clamp(p, opts->min_position, opts->max_position);
    }
    return HP_OK;
}

void hp_momentum_opts_default(hp_momentum_opts *opts) {
    if (!opts) return;
    opts->lookback = 60;
    opts->threshold = 0.0;
    opts->long_position = 1.0;
    opts->flat_position = 0.0;
}

hp_status hp_signal_momentum(const double *ret, size_t T,
                             const hp_momentum_opts *opts, double *pos) {
    if (!ret || !opts || !pos || T == 0) return HP_ERR_INVALID;
    if (opts->lookback == 0 || !isfinite(opts->threshold) ||
        !isfinite(opts->long_position) || !isfinite(opts->flat_position))
        return HP_ERR_INVALID;
    for (size_t t = 0; t < T; ++t) {
        pos[t] = NAN;
        if (t + 1 < opts->lookback) continue;
        double cum = 0.0;
        int ok = 1;
        for (size_t k = t + 1 - opts->lookback; k <= t; ++k) {
            if (!isfinite(ret[k])) {
                ok = 0;
                break;
            }
            cum += ret[k];
        }
        if (ok)
            pos[t] = cum > opts->threshold ? opts->long_position
                                           : opts->flat_position;
    }
    return HP_OK;
}

hp_status hp_signal_state_exposure(const double *filtered, size_t T, size_t n,
                                   const double *exposure, double *pos) {
    if (!filtered || !exposure || !pos || T == 0 || n == 0)
        return HP_ERR_INVALID;
    for (size_t t = 0; t < T; ++t) {
        double p = 0.0;
        for (size_t j = 0; j < n; ++j) {
            double g = filtered[t * n + j];
            if (!isfinite(g) || !isfinite(exposure[j])) {
                p = NAN;
                break;
            }
            p += g * exposure[j];
        }
        pos[t] = p;
    }
    return HP_OK;
}

hp_status hp_signal_multiply(const double *a, const double *b, size_t T,
                             double *out) {
    if (!a || !b || !out || T == 0) return HP_ERR_INVALID;
    for (size_t t = 0; t < T; ++t) {
        if (!isfinite(a[t]) || !isfinite(b[t]))
            out[t] = NAN;
        else
            out[t] = a[t] * b[t];
    }
    return HP_OK;
}

hp_status hp_signal_constant(size_t T, double position, double *pos) {
    if (!pos || T == 0 || !isfinite(position)) return HP_ERR_INVALID;
    for (size_t t = 0; t < T; ++t) pos[t] = position;
    return HP_OK;
}
