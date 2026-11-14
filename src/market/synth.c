#include "hp/market/synth.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "hp/core/csv.h"
#include "hp/core/io.h"
#include "hp/core/memory.h"
#include "hp/core/rng.h"
#include "hp/core/stats.h"

static int days_in_month(int y, int m) {
    static const int d[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (m == 2 && ((y % 4 == 0 && y % 100 != 0) || y % 400 == 0)) return 29;
    return d[m - 1];
}

static void next_weekday(int *y, int *m, int *d, int *dow) {
    do {
        ++(*d);
        if (*d > days_in_month(*y, *m)) {
            *d = 1;
            ++(*m);
            if (*m > 12) {
                *m = 1;
                ++(*y);
            }
        }
        *dow = (*dow + 1) % 7;
    } while (*dow == 0 || *dow == 6);
}

void hp_synth_opts_default(hp_synth_opts *opts) {
    if (!opts) return;
    opts->n_days = 3000;
    opts->n_assets = 10;
    opts->n_regimes = 3;
    opts->seed = 20240101;
}

static void stationary(const double *a, size_t n, double *pi) {
    for (size_t i = 0; i < n; ++i) pi[i] = 1.0 / (double)n;
    double *next = hp_doubles(n);
    for (int it = 0; it < 5000; ++it) {
        for (size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (size_t i = 0; i < n; ++i) s += pi[i] * a[i * n + j];
            next[j] = s;
        }
        memcpy(pi, next, n * sizeof(double));
    }
    hp_free(next);
}

hp_synth_data *hp_synth_generate(const hp_synth_opts *opts) {
    if (!opts || opts->n_days < 100 || opts->n_assets < 2 ||
        (opts->n_regimes != 2 && opts->n_regimes != 3))
        return NULL;
    const size_t T = opts->n_days, N = opts->n_assets;
    const size_t R = opts->n_regimes;

    double A[9] = {0}, mu[3] = {0}, sig[3] = {0}, idio_mult[3] = {0};
    if (R == 3) {
        double a[9] = {0.990, 0.009, 0.001,
                       0.018, 0.965, 0.017,
                       0.006, 0.018, 0.976};
        memcpy(A, a, sizeof(a));
        mu[0] = 0.15;
        mu[1] = 0.02;
        mu[2] = -0.25;
        sig[0] = 0.14;
        sig[1] = 0.20;
        sig[2] = 0.38;
        idio_mult[0] = 1.00;
        idio_mult[1] = 0.80;
        idio_mult[2] = 0.45;
    } else {
        double a[4] = {0.990, 0.010, 0.040, 0.960};
        memcpy(A, a, sizeof(a));
        mu[0] = 0.10;
        mu[1] = -0.20;
        sig[0] = 0.16;
        sig[1] = 0.34;
        idio_mult[0] = 1.0;
        idio_mult[1] = 0.45;
    }

    hp_synth_data *out = (hp_synth_data *)hp_xcalloc(1, sizeof(*out));
    out->n_regimes = R;
    out->panel = hp_panel_new(T, N);
    out->regime = hp_u8s(T);
    out->market_ret = hp_doubles(T);
    out->factor_beta = hp_doubles(N);

    for (size_t j = 0; j < N; ++j) {
        char name[16];
        snprintf(name, sizeof(name), "A%zu", j + 1);
        out->panel->symbols[j] = hp_strdup(name);
    }

    int y = 2012, m = 1, d = 2, dow = 1; /* Monday */
    for (size_t t = 0; t < T; ++t) {
        char buf[16];
        snprintf(buf, sizeof(buf), "%04d-%02d-%02d", y, m, d);
        out->panel->dates[t] = hp_strdup(buf);
        next_weekday(&y, &m, &d, &dow);
    }

    hp_rng rng;
    hp_rng_seed_stream(&rng, opts->seed, 1);
    hp_rng rng_assets;
    hp_rng_seed_stream(&rng_assets, opts->seed, 2);

    double pi[3] = {0};
    stationary(A, R, pi);

    /* Draw the latent regime path first; asset correlations then depend on
     * the current regime (idiosyncratic risk falls in the bear state). */
    double u0 = hp_rng_uniform(&rng);
    uint8_t state = 0;
    double cum = 0.0;
    for (size_t r = 0; r < R; ++r) {
        cum += pi[r];
        if (u0 <= cum) {
            state = (uint8_t)r;
            break;
        }
        state = (uint8_t)r;
    }
    for (size_t t = 0; t < T; ++t) {
        if (t > 0) {
            double u = hp_rng_uniform(&rng);
            double c = 0.0;
            uint8_t next = state;
            for (size_t r = 0; r < R; ++r) {
                c += A[state * R + r];
                if (u <= c) {
                    next = (uint8_t)r;
                    break;
                }
            }
            state = next;
        }
        out->regime[t] = state;
    }

    for (size_t j = 0; j < N; ++j) {
        out->factor_beta[j] = 0.55 + 0.90 * hp_rng_uniform(&rng_assets);
    }
    double *idio_ratio = hp_doubles(N);
    for (size_t j = 0; j < N; ++j)
        idio_ratio[j] = 0.35 + 0.45 * hp_rng_uniform(&rng_assets);

    const double sqrt252 = sqrt(252.0);
    const double tscale = sqrt(5.0 / 3.0); /* unit variance for df=5 */

    double *price = hp_doubles(N);
    for (size_t j = 0; j < N; ++j) price[j] = 100.0;

    for (size_t t = 0; t < T; ++t) {
        state = out->regime[t];
        double mret = mu[state] / 252.0;
        double msig = sig[state] / sqrt252;
        double shock = hp_rng_student_t(&rng, 5.0) / tscale;
        shock = hp_clamp(shock, -8.0, 8.0);
        double market = mret + msig * shock;
        out->market_ret[t] = market;

        for (size_t j = 0; j < N; ++j) {
            double beta = out->factor_beta[j];
            double idio_sig = msig * idio_ratio[j] * idio_mult[state];
            double z = hp_rng_normal(&rng_assets);
            double r = beta * market + idio_sig * z;
            r = hp_clamp(r, -0.25, 0.25);
            price[j] *= exp(r);
            out->panel->close[t * N + j] = price[j];
        }
        out->panel->ret[t * N + 0] = NAN;
        if (t > 0) {
            for (size_t j = 0; j < N; ++j) {
                double a = out->panel->close[(t - 1) * N + j];
                double b = out->panel->close[t * N + j];
                out->panel->ret[t * N + j] =
                    (a > 0.0 && b > 0.0) ? log(b / a) : NAN;
            }
        }
    }
    for (size_t t = 0; t < T; ++t)
        for (size_t j = 0; j < N; ++j)
            if (!isfinite(out->panel->ret[t * N + j]))
                out->panel->ret[t * N + j] = NAN;

    hp_free(idio_ratio);
    hp_free(price);
    return out;
}

void hp_synth_free(hp_synth_data *d) {
    if (!d) return;
    hp_panel_free(d->panel);
    hp_free(d->regime);
    hp_free(d->market_ret);
    hp_free(d->factor_beta);
    hp_free(d);
}

hp_status hp_synth_write_truth(const hp_synth_data *d, const char *path) {
    if (!d || !d->panel || !path) return HP_ERR_INVALID;
    const hp_panel *p = d->panel;
    size_t n_names = 3 + p->N;
    const char **names = (const char **)hp_xcalloc(n_names, sizeof(char *));
    names[0] = "date";
    names[1] = "regime";
    names[2] = "market_ret";
    for (size_t j = 0; j < p->N; ++j) names[3 + j] = p->symbols[j];
    hp_csv_writer *w = hp_csv_writer_open(path, names, n_names, 1);
    hp_free(names);
    if (!w) return HP_ERR_IO;
    double *values = hp_doubles(n_names - 1);
    for (size_t t = 0; t < p->T; ++t) {
        values[0] = (double)d->regime[t];
        values[1] = d->market_ret[t];
        for (size_t j = 0; j < p->N; ++j) values[2 + j] = p->ret[t * p->N + j];
        hp_status st = hp_csv_writer_row_label(w, p->dates[t], values, n_names - 1);
        if (st != HP_OK) {
            hp_csv_writer_abort(w);
            hp_free(values);
            return st;
        }
    }
    hp_free(values);
    return hp_csv_writer_close(w);
}
