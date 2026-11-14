#include "hp/market/features.h"

#include <math.h>
#include <string.h>

#include "hp/core/memory.h"

hp_features *hp_features_from_returns(const double *ret, size_t T,
                                      size_t vol_window) {
    if (!ret || T == 0 || vol_window < 2) return NULL;
    hp_features *f = (hp_features *)hp_xcalloc(1, sizeof(*f));
    f->T = T;
    f->dim = 2;
    f->start = T;
    f->ret = hp_doubles(T);
    f->vol = hp_doubles(T);
    f->x = hp_doubles(T * 2);
    memcpy(f->ret, ret, T * sizeof(double));
    for (size_t t = 0; t < T; ++t) {
        f->vol[t] = NAN;
        f->x[t * 2 + 0] = NAN;
        f->x[t * 2 + 1] = NAN;
    }
    for (size_t t = vol_window; t < T; ++t) {
        size_t cnt = 0;
        double s = 0.0, s2 = 0.0;
        for (size_t k = t - vol_window + 1; k <= t; ++k) {
            double r = ret[k];
            if (!isfinite(r)) continue;
            ++cnt;
            s += r;
            s2 += r * r;
        }
        if (cnt != vol_window) continue;
        double mean = s / (double)vol_window;
        double var = (s2 - (double)vol_window * mean * mean) /
                     (double)(vol_window - 1);
        if (var < 0.0) var = 0.0;
        double vol = sqrt(var * 252.0);
        if (vol < 1e-4) vol = 1e-4;
        f->vol[t] = vol;
        f->x[t * 2 + 0] = ret[t] * 100.0;
        f->x[t * 2 + 1] = log(vol * 100.0);
        if (f->start == T) f->start = t;
    }
    return f;
}

hp_features *hp_features_from_panel(const hp_panel *p, size_t asset_idx,
                                    size_t vol_window) {
    if (!p || asset_idx >= p->N) return NULL;
    double *ret = hp_doubles(p->T);
    hp_status st = hp_panel_return_series(p, asset_idx, ret);
    if (st != HP_OK) {
        hp_free(ret);
        return NULL;
    }
    hp_features *f = hp_features_from_returns(ret, p->T, vol_window);
    hp_free(ret);
    return f;
}

void hp_features_free(hp_features *f) {
    if (!f) return;
    hp_free(f->ret);
    hp_free(f->vol);
    hp_free(f->x);
    hp_free(f);
}
