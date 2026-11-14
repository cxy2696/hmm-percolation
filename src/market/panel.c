#include "hp/market/panel.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "hp/core/csv.h"
#include "hp/core/io.h"
#include "hp/core/memory.h"

hp_panel *hp_panel_new(size_t T, size_t N) {
    if (T == 0 || N == 0) return NULL;
    hp_panel *p = (hp_panel *)hp_xcalloc(1, sizeof(*p));
    p->T = T;
    p->N = N;
    p->dates = (char **)hp_xcalloc(T, sizeof(char *));
    p->symbols = (char **)hp_xcalloc(N, sizeof(char *));
    p->close = hp_doubles(T * N);
    p->ret = hp_doubles(T * N);
    for (size_t i = 0; i < T * N; ++i) {
        p->close[i] = NAN;
        p->ret[i] = NAN;
    }
    return p;
}

void hp_panel_free(hp_panel *p) {
    if (!p) return;
    if (p->dates) {
        for (size_t t = 0; t < p->T; ++t) hp_free(p->dates[t]);
        hp_free(p->dates);
    }
    if (p->symbols) {
        for (size_t j = 0; j < p->N; ++j) hp_free(p->symbols[j]);
        hp_free(p->symbols);
    }
    hp_free(p->close);
    hp_free(p->ret);
    hp_free(p);
}

hp_panel *hp_panel_load_csv(const char *path, hp_status *status) {
    if (status) *status = HP_ERR_INVALID;
    if (!path) return NULL;
    hp_csv csv;
    hp_status st = hp_csv_read(path, &csv);
    if (st != HP_OK) {
        if (status) *status = st;
        return NULL;
    }
    if (csv.n_cols < 2 || csv.n_rows < 2) {
        hp_csv_free(&csv);
        if (status) *status = HP_ERR_FORMAT;
        return NULL;
    }
    const size_t T = csv.n_rows, N = csv.n_cols - 1;
    hp_panel *p = hp_panel_new(T, N);
    for (size_t j = 0; j < N; ++j) p->symbols[j] = hp_strdup(csv.header[j + 1]);
    for (size_t t = 0; t < T; ++t) {
        p->dates[t] = hp_strdup(hp_csv_cell(&csv, t, 0));
        for (size_t j = 0; j < N; ++j) {
            const char *cell = hp_csv_cell(&csv, t, j + 1);
            double v = NAN;
            st = hp_csv_parse_double(cell, &v);
            if (st != HP_OK) {
                hp_csv_free(&csv);
                hp_panel_free(p);
                if (status) *status = HP_ERR_FORMAT;
                return NULL;
            }
            p->close[t * N + j] = v;
        }
    }
    hp_csv_free(&csv);
    st = hp_panel_compute_returns(p);
    if (st != HP_OK) {
        hp_panel_free(p);
        if (status) *status = st;
        return NULL;
    }
    if (status) *status = HP_OK;
    return p;
}

hp_status hp_panel_compute_returns(hp_panel *p) {
    if (!p || !p->close || !p->ret) return HP_ERR_INVALID;
    const size_t T = p->T, N = p->N;
    for (size_t j = 0; j < N; ++j) p->ret[j] = NAN;
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < N; ++j) {
            double a = p->close[(t - 1) * N + j];
            double b = p->close[t * N + j];
            p->ret[t * N + j] =
                (isfinite(a) && isfinite(b) && a > 0.0 && b > 0.0)
                    ? log(b / a)
                    : NAN;
        }
    }
    return HP_OK;
}

int hp_panel_symbol_index(const hp_panel *p, const char *symbol) {
    if (!p || !symbol) return -1;
    for (size_t j = 0; j < p->N; ++j)
        if (strcmp(p->symbols[j], symbol) == 0) return (int)j;
    for (size_t j = 0; j < p->N; ++j)
        if (hp_str_ieq(p->symbols[j], symbol)) return (int)j;
    return -1;
}

hp_status hp_panel_return_series(const hp_panel *p, size_t j, double *out) {
    if (!p || !out || j >= p->N) return HP_ERR_INVALID;
    for (size_t t = 0; t < p->T; ++t) out[t] = p->ret[t * p->N + j];
    return HP_OK;
}

hp_status hp_panel_ew_returns(const hp_panel *p, double *out) {
    if (!p || !out) return HP_ERR_INVALID;
    for (size_t t = 0; t < p->T; ++t) {
        double s = 0.0;
        size_t cnt = 0;
        for (size_t j = 0; j < p->N; ++j) {
            double r = p->ret[t * p->N + j];
            if (isfinite(r)) {
                s += r;
                ++cnt;
            }
        }
        out[t] = cnt ? s / (double)cnt : NAN;
    }
    return HP_OK;
}

hp_status hp_panel_save_csv(const hp_panel *p, const char *path) {
    if (!p || !path) return HP_ERR_INVALID;
    const char **names = (const char **)hp_xcalloc(p->N + 1, sizeof(char *));
    names[0] = "date";
    for (size_t j = 0; j < p->N; ++j) names[j + 1] = p->symbols[j];
    hp_csv_writer *w = hp_csv_writer_open(path, names, p->N + 1, 1);
    hp_free(names);
    if (!w) return HP_ERR_IO;
    for (size_t t = 0; t < p->T; ++t) {
        hp_status st = hp_csv_writer_row_label(w, p->dates[t], p->close + t * p->N,
                                               p->N);
        if (st != HP_OK) {
            hp_csv_writer_abort(w);
            return st;
        }
    }
    return hp_csv_writer_close(w);
}

hp_status hp_panel_write_returns_csv(const hp_panel *p, const char *path) {
    if (!p || !path) return HP_ERR_INVALID;
    const char **names = (const char **)hp_xcalloc(p->N + 1, sizeof(char *));
    names[0] = "date";
    for (size_t j = 0; j < p->N; ++j) names[j + 1] = p->symbols[j];
    hp_csv_writer *w = hp_csv_writer_open(path, names, p->N + 1, 1);
    hp_free(names);
    if (!w) return HP_ERR_IO;
    for (size_t t = 0; t < p->T; ++t) {
        hp_status st = hp_csv_writer_row_label(w, p->dates[t], p->ret + t * p->N,
                                               p->N);
        if (st != HP_OK) {
            hp_csv_writer_abort(w);
            return st;
        }
    }
    return hp_csv_writer_close(w);
}
