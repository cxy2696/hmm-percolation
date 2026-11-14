#include "hp/hmm/ghmm.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/io.h"
#include "hp/core/log.h"
#include "hp/core/memory.h"
#include "hp/core/stats.h"

#define HP_GHMM_FLOOR 1e-8
#define HP_GHMM_LOG_FLOOR (-745.0)

static double safe_log(double x) { return x > 0.0 ? log(x) : -INFINITY; }

static double lse2(double a, double b) {
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    if (a > b) return a + log1p(exp(b - a));
    return b + log1p(exp(a - b));
}

hp_ghmm *hp_ghmm_new(size_t n_states, size_t dim, size_t n_mix) {
    if (n_states == 0 || dim == 0 || n_mix == 0) return NULL;
    hp_ghmm *m = (hp_ghmm *)hp_xcalloc(1, sizeof(*m));
    m->n_states = n_states;
    m->dim = dim;
    m->n_mix = n_mix;
    m->pi = hp_doubles(n_states);
    m->a = hp_doubles(n_states * n_states);
    m->mix = hp_doubles(n_states * n_mix);
    m->mean = hp_doubles(n_states * n_mix * dim);
    m->var = hp_doubles(n_states * n_mix * dim);
    hp_ghmm_set_uniform(m);
    return m;
}

hp_ghmm *hp_ghmm_clone(const hp_ghmm *m) {
    if (!m) return NULL;
    hp_ghmm *c = hp_ghmm_new(m->n_states, m->dim, m->n_mix);
    if (!c) return NULL;
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    memcpy(c->pi, m->pi, n * sizeof(double));
    memcpy(c->a, m->a, n * n * sizeof(double));
    memcpy(c->mix, m->mix, n * k * sizeof(double));
    memcpy(c->mean, m->mean, n * k * d * sizeof(double));
    memcpy(c->var, m->var, n * k * d * sizeof(double));
    return c;
}

void hp_ghmm_free(hp_ghmm *m) {
    if (!m) return;
    hp_free(m->pi);
    hp_free(m->a);
    hp_free(m->mix);
    hp_free(m->mean);
    hp_free(m->var);
    hp_free(m);
}

hp_status hp_ghmm_resize(hp_ghmm *m, size_t n_states, size_t dim, size_t n_mix) {
    if (!m || n_states == 0 || dim == 0 || n_mix == 0) return HP_ERR_INVALID;
    hp_free(m->pi);
    hp_free(m->a);
    hp_free(m->mix);
    hp_free(m->mean);
    hp_free(m->var);
    m->n_states = n_states;
    m->dim = dim;
    m->n_mix = n_mix;
    m->pi = hp_doubles(n_states);
    m->a = hp_doubles(n_states * n_states);
    m->mix = hp_doubles(n_states * n_mix);
    m->mean = hp_doubles(n_states * n_mix * dim);
    m->var = hp_doubles(n_states * n_mix * dim);
    return hp_ghmm_set_uniform(m);
}

hp_status hp_ghmm_set_uniform(hp_ghmm *m) {
    if (!m || !m->pi || !m->a || !m->mix || !m->mean || !m->var)
        return HP_ERR_INVALID;
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    for (size_t i = 0; i < n; ++i) {
        m->pi[i] = 1.0 / (double)n;
        for (size_t j = 0; j < n; ++j) m->a[i * n + j] = 1.0 / (double)n;
        for (size_t c = 0; c < k; ++c) {
            m->mix[i * k + c] = 1.0 / (double)k;
            for (size_t dd = 0; dd < d; ++dd) {
                m->mean[(i * k + c) * d + dd] = 0.0;
                m->var[(i * k + c) * d + dd] = 1.0;
            }
        }
    }
    return HP_OK;
}

hp_status hp_ghmm_validate(const hp_ghmm *m) {
    if (!m || m->n_states == 0 || m->dim == 0 || m->n_mix == 0)
        return HP_ERR_INVALID;
    if (!m->pi || !m->a || !m->mix || !m->mean || !m->var)
        return HP_ERR_INVALID;
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double v = m->pi[i];
        if (!isfinite(v) || v < -1e-12 || v > 1.0 + 1e-9) return HP_ERR_DOMAIN;
        s += v;
    }
    if (fabs(s - 1.0) > 1e-6) return HP_ERR_DOMAIN;
    for (size_t i = 0; i < n; ++i) {
        s = 0.0;
        for (size_t j = 0; j < n; ++j) {
            double v = m->a[i * n + j];
            if (!isfinite(v) || v < -1e-12 || v > 1.0 + 1e-9)
                return HP_ERR_DOMAIN;
            s += v;
        }
        if (fabs(s - 1.0) > 1e-6) return HP_ERR_DOMAIN;
    }
    for (size_t j = 0; j < n; ++j) {
        s = 0.0;
        for (size_t c = 0; c < k; ++c) {
            double v = m->mix[j * k + c];
            if (!isfinite(v) || v < -1e-12 || v > 1.0 + 1e-9)
                return HP_ERR_DOMAIN;
            s += v;
        }
        if (fabs(s - 1.0) > 1e-6) return HP_ERR_DOMAIN;
        for (size_t c = 0; c < k; ++c) {
            for (size_t dd = 0; dd < d; ++dd) {
                double mu = m->mean[(j * k + c) * d + dd];
                double var = m->var[(j * k + c) * d + dd];
                if (!isfinite(mu) || !isfinite(var) || !(var > 0.0))
                    return HP_ERR_DOMAIN;
            }
        }
    }
    return HP_OK;
}

size_t hp_ghmm_n_params(size_t n_states, size_t dim, size_t n_mix) {
    if (n_states == 0 || dim == 0 || n_mix == 0) return 0;
    size_t p = (n_states - 1) + n_states * (n_states - 1);
    p += n_states * (n_mix - 1);
    p += n_states * n_mix * dim; /* means */
    p += n_states * n_mix * dim; /* diagonal variances */
    return p;
}

double hp_ghmm_state_mean(const hp_ghmm *m, size_t state, size_t dim) {
    if (!m || state >= m->n_states || dim >= m->dim) return NAN;
    double mu = 0.0;
    for (size_t c = 0; c < m->n_mix; ++c)
        mu += m->mix[state * m->n_mix + c] *
              m->mean[(state * m->n_mix + c) * m->dim + dim];
    return mu;
}

static void emission_row(const hp_ghmm *m, const double *x, double *logc,
                         double *logb, double *terms) {
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    for (size_t j = 0; j < n; ++j) {
        for (size_t c = 0; c < k; ++c) {
            double lp = safe_log(m->mix[j * k + c]);
            for (size_t dd = 0; dd < d; ++dd) {
                if (lp == -INFINITY) break;
                double var = m->var[(j * k + c) * d + dd];
                double mu = m->mean[(j * k + c) * d + dd];
                if (!(var > 0.0)) {
                    lp = -INFINITY;
                    break;
                }
                double diff = x[dd] - mu;
                lp += -0.5 * (HP_LOG_2PI + log(var) + diff * diff / var);
            }
            if (logc) logc[j * k + c] = lp;
            terms[c] = lp;
        }
        logb[j] = hp_logsumexp(terms, k);
    }
}

static void compute_tables(const hp_ghmm *m, const double *x, size_t T,
                           double *logb, double *logc) {
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    double *terms = hp_doubles(k);
    for (size_t t = 0; t < T; ++t) {
        emission_row(m, x + t * d, logc ? logc + t * n * k : NULL,
                     logb + t * n, terms);
    }
    hp_free(terms);
}

static double forward_from_logb(const hp_ghmm *m, const double *logb, size_t T,
                                double *la) {
    const size_t n = m->n_states;
    for (size_t j = 0; j < n; ++j)
        la[j] = safe_log(m->pi[j]) + logb[j];
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < n; ++j) {
            double acc = -INFINITY;
            for (size_t i = 0; i < n; ++i) {
                double aij = m->a[i * n + j];
                if (aij <= 0.0) continue;
                acc = lse2(acc, la[(t - 1) * n + i] + log(aij));
            }
            la[t * n + j] = acc + logb[t * n + j];
        }
    }
    return hp_logsumexp(la + (T - 1) * n, n);
}

static void backward_from_logb(const hp_ghmm *m, const double *logb, size_t T,
                               double *lb) {
    const size_t n = m->n_states;
    for (size_t i = 0; i < n; ++i) lb[(T - 1) * n + i] = 0.0;
    for (size_t t = T - 1; t-- > 0;) {
        for (size_t i = 0; i < n; ++i) {
            double acc = -INFINITY;
            for (size_t j = 0; j < n; ++j) {
                double aij = m->a[i * n + j];
                if (aij <= 0.0) continue;
                double v = log(aij) + logb[(t + 1) * n + j] + lb[(t + 1) * n + j];
                acc = lse2(acc, v);
            }
            lb[t * n + i] = acc;
        }
    }
}

static int x_valid(const double *x, size_t T, size_t d) {
    if (!x) return 0;
    for (size_t i = 0; i < T * d; ++i)
        if (!isfinite(x[i])) return 0;
    return 1;
}

hp_status hp_ghmm_log_emission(const hp_ghmm *m, const double *x, double *logb) {
    if (!m || !x || !logb) return HP_ERR_INVALID;
    HP_TRY(hp_ghmm_validate(m));
    for (size_t i = 0; i < m->dim; ++i)
        if (!isfinite(x[i])) return HP_ERR_DOMAIN;
    double *terms = hp_doubles(m->n_mix);
    emission_row(m, x, NULL, logb, terms);
    hp_free(terms);
    return HP_OK;
}

hp_status hp_ghmm_forward_log(const hp_ghmm *m, const double *x, size_t T,
                              double *logalpha, double *loglik) {
    if (!m || !x || !logalpha || !loglik || T == 0) return HP_ERR_INVALID;
    HP_TRY(hp_ghmm_validate(m));
    if (!x_valid(x, T, m->dim)) return HP_ERR_DOMAIN;
    const size_t n = m->n_states;
    double *logb = hp_doubles(T * n);
    compute_tables(m, x, T, logb, NULL);
    double ll = forward_from_logb(m, logb, T, logalpha);
    hp_free(logb);
    *loglik = ll;
    return isfinite(ll) ? HP_OK : HP_ERR_DOMAIN;
}

hp_status hp_ghmm_backward_log(const hp_ghmm *m, const double *x, size_t T,
                               double *logbeta) {
    if (!m || !x || !logbeta || T == 0) return HP_ERR_INVALID;
    HP_TRY(hp_ghmm_validate(m));
    if (!x_valid(x, T, m->dim)) return HP_ERR_DOMAIN;
    const size_t n = m->n_states;
    double *logb = hp_doubles(T * n);
    compute_tables(m, x, T, logb, NULL);
    backward_from_logb(m, logb, T, logbeta);
    hp_free(logb);
    return HP_OK;
}

hp_status hp_ghmm_posterior(const hp_ghmm *m, const double *x, size_t T,
                            double *gamma, double *loglik) {
    if (!m || !gamma || T == 0) return HP_ERR_INVALID;
    const size_t n = m->n_states;
    double *la = hp_doubles(T * n);
    double *lb = hp_doubles(T * n);
    double ll = 0;
    hp_status st = hp_ghmm_forward_log(m, x, T, la, &ll);
    if (st == HP_OK) st = hp_ghmm_backward_log(m, x, T, lb);
    if (st == HP_OK) {
        for (size_t t = 0; t < T; ++t) {
            double *g = gamma + t * n;
            double s = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double v = exp(la[t * n + i] + lb[t * n + i] - ll);
                if (!isfinite(v) || v < 0.0) v = 0.0;
                g[i] = v;
                s += v;
            }
            if (!(s > 0.0)) {
                st = HP_ERR_DOMAIN;
                break;
            }
            for (size_t i = 0; i < n; ++i) g[i] /= s;
        }
    }
    if (loglik) *loglik = ll;
    hp_free(la);
    hp_free(lb);
    return st;
}

hp_status hp_ghmm_filter(const hp_ghmm *m, const double *x, size_t T,
                         double *filt, double *loglik) {
    if (!m || !filt || T == 0) return HP_ERR_INVALID;
    const size_t n = m->n_states;
    double *la = hp_doubles(T * n);
    double ll = 0;
    hp_status st = hp_ghmm_forward_log(m, x, T, la, &ll);
    if (st == HP_OK) {
        for (size_t t = 0; t < T; ++t) {
            double lse = hp_logsumexp(la + t * n, n);
            if (!isfinite(lse)) {
                st = HP_ERR_DOMAIN;
                break;
            }
            for (size_t i = 0; i < n; ++i)
                filt[t * n + i] = exp(la[t * n + i] - lse);
        }
    }
    if (loglik) *loglik = ll;
    hp_free(la);
    return st;
}

hp_status hp_ghmm_filter_step(const hp_ghmm *m, const double *x, double *filt,
                              double *loglik_step) {
    if (!m || !x || !filt) return HP_ERR_INVALID;
    HP_TRY(hp_ghmm_validate(m));
    for (size_t dd = 0; dd < m->dim; ++dd)
        if (!isfinite(x[dd])) return HP_ERR_DOMAIN;
    const size_t n = m->n_states;
    double *logb = hp_doubles(n);
    double *out = hp_doubles(n);
    double *terms = hp_doubles(m->n_mix);
    emission_row(m, x, NULL, logb, terms);
    for (size_t j = 0; j < n; ++j) {
        double acc = -INFINITY;
        for (size_t i = 0; i < n; ++i) {
            if (filt[i] <= 0.0) continue;
            double aij = m->a[i * n + j];
            if (aij <= 0.0) continue;
            acc = lse2(acc, log(filt[i]) + log(aij));
        }
        out[j] = acc + logb[j];
    }
    double lse = hp_logsumexp(out, n);
    if (!isfinite(lse)) {
        hp_free(logb);
        hp_free(out);
        hp_free(terms);
        return HP_ERR_DOMAIN;
    }
    for (size_t j = 0; j < n; ++j) filt[j] = exp(out[j] - lse);
    if (loglik_step) *loglik_step = lse;
    hp_free(logb);
    hp_free(out);
    hp_free(terms);
    return HP_OK;
}

hp_status hp_ghmm_viterbi(const hp_ghmm *m, const double *x, size_t T,
                          uint32_t *path, double *logp) {
    if (!m || !x || !path || T == 0) return HP_ERR_INVALID;
    HP_TRY(hp_ghmm_validate(m));
    if (!x_valid(x, T, m->dim)) return HP_ERR_DOMAIN;
    const size_t n = m->n_states;
    double *logb = hp_doubles(T * n);
    compute_tables(m, x, T, logb, NULL);
    double *delta = hp_doubles(n);
    double *next = hp_doubles(n);
    uint32_t *psi = hp_u32s(T * n);
    for (size_t j = 0; j < n; ++j)
        delta[j] = safe_log(m->pi[j]) + logb[j];
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < n; ++j) {
            double best = -INFINITY;
            uint32_t arg = 0;
            for (size_t i = 0; i < n; ++i) {
                double aij = m->a[i * n + j];
                if (aij <= 0.0 || delta[i] == -INFINITY) continue;
                double v = delta[i] + log(aij);
                if (v > best) {
                    best = v;
                    arg = (uint32_t)i;
                }
            }
            next[j] = best + logb[t * n + j];
            psi[t * n + j] = arg;
        }
        double *tmp = delta;
        delta = next;
        next = tmp;
    }
    double best = -INFINITY;
    uint32_t arg = 0;
    for (size_t j = 0; j < n; ++j) {
        if (delta[j] > best) {
            best = delta[j];
            arg = (uint32_t)j;
        }
    }
    if (!isfinite(best)) {
        hp_free(logb);
        hp_free(delta);
        hp_free(next);
        hp_free(psi);
        return HP_ERR_DOMAIN;
    }
    path[T - 1] = arg;
    for (size_t t = T - 1; t-- > 0;) path[t] = psi[(t + 1) * n + path[t + 1]];
    if (logp) *logp = best;
    hp_free(logb);
    hp_free(delta);
    hp_free(next);
    hp_free(psi);
    return HP_OK;
}

void hp_ghmm_train_opts_default(hp_ghmm_train_opts *opts) {
    if (!opts) return;
    opts->max_iter = 200;
    opts->tol = 1e-6;
    opts->n_restarts = 4;
    opts->seed = 42;
    opts->alpha = 1.0;
    opts->var_floor = 1e-4;
    opts->mix_floor = 1e-8;
    opts->init = HP_GHMM_INIT_QUANTILE;
    opts->verbose = 0;
}

hp_status hp_ghmm_randomize(hp_ghmm *m, hp_rng *rng, double alpha,
                            const double *global_var) {
    if (!m || !rng || !(alpha > 0.0)) return HP_ERR_INVALID;
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    const size_t nmax = n > k ? n : k;
    double *dir = hp_doubles(nmax > d ? nmax : d);
    for (size_t i = 0; i < n; ++i) dir[i] = alpha;
    hp_rng_dirichlet(rng, dir, n, m->pi);
    for (size_t i = 0; i < n; ++i) hp_rng_dirichlet(rng, dir, n, m->a + i * n);
    for (size_t j = 0; j < n; ++j) {
        for (size_t c = 0; c < k; ++c) dir[c] = alpha;
        hp_rng_dirichlet(rng, dir, k, m->mix + j * k);
        for (size_t c = 0; c < k; ++c) {
            for (size_t dd = 0; dd < d; ++dd) {
                double sd = global_var ? sqrt(global_var[dd]) : 1.0;
                m->mean[(j * k + c) * d + dd] = hp_rng_normal(rng) * sd;
                double v = global_var ? global_var[dd] * (0.5 + hp_rng_uniform(rng))
                                      : (0.5 + hp_rng_uniform(rng));
                m->var[(j * k + c) * d + dd] = v > 1e-6 ? v : 1e-6;
            }
        }
    }
    hp_free(dir);
    return HP_OK;
}

typedef struct {
    double key;
    size_t idx;
} ghmm_pair;

static int cmp_pair(const void *a, const void *b) {
    double x = ((const ghmm_pair *)a)->key;
    double y = ((const ghmm_pair *)b)->key;
    return (x > y) - (x < y);
}

/* Quantile initialisation: sort by the first feature, cut into N contiguous
 * blocks, seed each state with one block.  Deterministic and usually puts
 * regime states on the right side of the return distribution immediately. */
static void init_quantile(hp_ghmm *m, const double *x, size_t T,
                          double var_floor) {
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    ghmm_pair *order = (ghmm_pair *)hp_xmalloc(T * sizeof(ghmm_pair));
    for (size_t t = 0; t < T; ++t) {
        order[t].key = x[t * d];
        order[t].idx = t;
    }
    qsort(order, T, sizeof(ghmm_pair), cmp_pair);

    for (size_t i = 0; i < n; ++i) m->pi[i] = 1.0 / (double)n;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            m->a[i * n + j] = (i == j) ? 0.9
                                       : (n > 1 ? 0.1 / (double)(n - 1) : 0.0);
        if (n == 1) m->a[i] = 1.0;
    }

    for (size_t j = 0; j < n; ++j) {
        size_t start = (j * T) / n;
        size_t end = ((j + 1) * T) / n;
        if (end <= start) end = start + 1;
        if (end > T) end = T;
        size_t cnt = end - start;
        double *mu = hp_doubles(d), *var = hp_doubles(d);
        for (size_t dd = 0; dd < d; ++dd) {
            double s = 0.0, s2 = 0.0;
            for (size_t q = start; q < end; ++q) {
                double v = x[order[q].idx * d + dd];
                s += v;
                s2 += v * v;
            }
            double mean = s / (double)cnt;
            double vv = s2 / (double)cnt - mean * mean;
            if (vv < 0.0) vv = 0.0;
            mu[dd] = mean;
            var[dd] = vv > var_floor ? vv : var_floor;
        }
        for (size_t c = 0; c < k; ++c) {
            m->mix[j * k + c] = 1.0 / (double)k;
            double shift = k > 1 ? ((double)c - 0.5 * (double)(k - 1)) * 0.35 : 0.0;
            for (size_t dd = 0; dd < d; ++dd) {
                m->mean[(j * k + c) * d + dd] =
                    mu[dd] + shift * sqrt(var[dd]);
                m->var[(j * k + c) * d + dd] = var[dd];
            }
        }
        hp_free(mu);
        hp_free(var);
    }
    hp_free(order);
}

static void ghmm_copy_params(hp_ghmm *dst, const hp_ghmm *src) {
    const size_t n = src->n_states, k = src->n_mix, d = src->dim;
    memcpy(dst->pi, src->pi, n * sizeof(double));
    memcpy(dst->a, src->a, n * n * sizeof(double));
    memcpy(dst->mix, src->mix, n * k * sizeof(double));
    memcpy(dst->mean, src->mean, n * k * d * sizeof(double));
    memcpy(dst->var, src->var, n * k * d * sizeof(double));
}

static void ghmm_canonicalize(hp_ghmm *m) {
    const size_t n = m->n_states, k = m->n_mix, d = m->dim;
    double *key = hp_doubles(n);
    size_t *order = hp_sizes(n);
    for (size_t j = 0; j < n; ++j) {
        key[j] = hp_ghmm_state_mean(m, j, 0);
        order[j] = j;
    }
    for (size_t i = 1; i < n; ++i) {
        size_t v = order[i];
        size_t q = i;
        while (q > 0) {
            double a = key[v], b = key[order[q - 1]];
            if (a > b || (a == b && v > order[q - 1])) break;
            order[q] = order[q - 1];
            --q;
        }
        order[q] = v;
    }
    double *pi = hp_doubles(n), *a = hp_doubles(n * n);
    double *mix = hp_doubles(n * k), *mean = hp_doubles(n * k * d),
           *var = hp_doubles(n * k * d);
    for (size_t j = 0; j < n; ++j) {
        pi[j] = m->pi[order[j]];
        for (size_t c = 0; c < k; ++c) {
            mix[j * k + c] = m->mix[order[j] * k + c];
            for (size_t dd = 0; dd < d; ++dd) {
                mean[(j * k + c) * d + dd] =
                    m->mean[(order[j] * k + c) * d + dd];
                var[(j * k + c) * d + dd] =
                    m->var[(order[j] * k + c) * d + dd];
            }
        }
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            a[i * n + j] = m->a[order[i] * n + order[j]];
    memcpy(m->pi, pi, n * sizeof(double));
    memcpy(m->a, a, n * n * sizeof(double));
    memcpy(m->mix, mix, n * k * sizeof(double));
    memcpy(m->mean, mean, n * k * d * sizeof(double));
    memcpy(m->var, var, n * k * d * sizeof(double));
    hp_free(key);
    hp_free(order);
    hp_free(pi);
    hp_free(a);
    hp_free(mix);
    hp_free(mean);
    hp_free(var);
}

hp_status hp_ghmm_report_fill(const hp_ghmm *m, size_t T, double loglik,
                              unsigned iterations, int converged,
                              unsigned restarts, hp_ghmm_report *report) {
    if (!m || !report || T == 0) return HP_ERR_INVALID;
    size_t p = hp_ghmm_n_params(m->n_states, m->dim, m->n_mix);
    report->loglik = loglik;
    report->n_params = p;
    report->bic = -2.0 * loglik + (double)p * log((double)T);
    report->aic = -2.0 * loglik + 2.0 * (double)p;
    report->iterations = iterations;
    report->converged = converged;
    report->restarts = restarts;
    return HP_OK;
}

static hp_status ghmm_em(const double *x, size_t T, const hp_ghmm *init,
                         const hp_ghmm_train_opts *opts, hp_ghmm *result,
                         unsigned *iterations, int *converged,
                         double *loglik_out) {
    const size_t n = init->n_states, k = init->n_mix, d = init->dim;
    hp_ghmm *cur = hp_ghmm_clone(init);
    hp_ghmm *best = hp_ghmm_clone(init);
    double *logb = hp_doubles(T * n);
    double *logc = hp_doubles(T * n * k);
    double *la = hp_doubles(T * n);
    double *lb = hp_doubles(T * n);
    double *gamma = hp_doubles(T * n);
    double *xi = hp_doubles(n * n);
    double *gsum = hp_doubles(n);
    double *gsumk = hp_doubles(n * k);
    double *sumx = hp_doubles(n * k * d);
    double *sumxx = hp_doubles(n * k * d);
    double *a_new = hp_doubles(n * n);
    double *mix_new = hp_doubles(n * k);
    double *mean_new = hp_doubles(n * k * d);
    double *var_new = hp_doubles(n * k * d);

    compute_tables(cur, x, T, logb, logc);
    double ll = forward_from_logb(cur, logb, T, la);
    double best_ll = ll;
    unsigned it = 0;
    int conv = 0;
    hp_status st = isfinite(ll) ? HP_OK : HP_ERR_DOMAIN;
    if (st != HP_OK) goto done;

    for (it = 1; it <= opts->max_iter; ++it) {
        backward_from_logb(cur, logb, T, lb);
        memset(xi, 0, n * n * sizeof(double));
        memset(gsum, 0, n * sizeof(double));
        memset(gsumk, 0, n * k * sizeof(double));
        memset(sumx, 0, n * k * d * sizeof(double));
        memset(sumxx, 0, n * k * d * sizeof(double));

        for (size_t t = 0; t < T; ++t) {
            double *g = gamma + t * n;
            double s = 0.0;
            for (size_t j = 0; j < n; ++j) {
                double v = exp(la[t * n + j] + lb[t * n + j] - ll);
                if (!isfinite(v) || v < 0.0) v = 0.0;
                g[j] = v;
                s += v;
            }
            if (!(s > 0.0)) {
                st = HP_ERR_DOMAIN;
                goto done;
            }
            for (size_t j = 0; j < n; ++j) {
                g[j] /= s;
                gsum[j] += g[j];
                double lbj = logb[t * n + j];
                if (g[j] <= 0.0 || lbj == -INFINITY) continue;
                for (size_t c = 0; c < k; ++c) {
                    double lr = logc[(t * n + j) * k + c] - lbj;
                    double gk = lr > HP_GHMM_LOG_FLOOR ? g[j] * exp(lr) : 0.0;
                    size_t idx = j * k + c;
                    gsumk[idx] += gk;
                    for (size_t dd = 0; dd < d; ++dd) {
                        double xv = x[t * d + dd];
                        sumx[idx * d + dd] += gk * xv;
                        sumxx[idx * d + dd] += gk * xv * xv;
                    }
                }
            }
            if (t + 1 < T) {
                for (size_t i = 0; i < n; ++i) {
                    for (size_t j = 0; j < n; ++j) {
                        double lbj = logb[(t + 1) * n + j];
                        if (lbj == -INFINITY) continue;
                        double lv = la[t * n + i] + safe_log(cur->a[i * n + j]) +
                                    lbj + lb[(t + 1) * n + j] - ll;
                        if (lv > HP_GHMM_LOG_FLOOR) xi[i * n + j] += exp(lv);
                    }
                }
            }
        }

        for (size_t j = 0; j < n; ++j) cur->pi[j] = gamma[j];
        double pis = 0.0;
        for (size_t j = 0; j < n; ++j) pis += cur->pi[j];
        for (size_t j = 0; j < n; ++j)
            cur->pi[j] = pis > 0.0 ? cur->pi[j] / pis : 1.0 / (double)n;

        for (size_t i = 0; i < n; ++i) {
            double den = gsum[i] - gamma[(T - 1) * n + i];
            double rowsum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                double v = den > 1e-300 ? xi[i * n + j] / den
                                        : cur->a[i * n + j];
                if (!(v >= 0.0) || !isfinite(v)) v = 0.0;
                if (v < HP_GHMM_FLOOR) v = HP_GHMM_FLOOR;
                a_new[i * n + j] = v;
                rowsum += v;
            }
            for (size_t j = 0; j < n; ++j) a_new[i * n + j] /= rowsum;
        }

        for (size_t j = 0; j < n; ++j) {
            double rowsum = 0.0;
            for (size_t c = 0; c < k; ++c) {
                double v = gsum[j] > 1e-300 ? gsumk[j * k + c] / gsum[j]
                                            : cur->mix[j * k + c];
                if (!(v >= 0.0) || !isfinite(v)) v = 0.0;
                if (v < opts->mix_floor) v = opts->mix_floor;
                mix_new[j * k + c] = v;
                rowsum += v;
            }
            for (size_t c = 0; c < k; ++c) mix_new[j * k + c] /= rowsum;
            for (size_t c = 0; c < k; ++c) {
                size_t idx = j * k + c;
                double den = gsumk[idx];
                for (size_t dd = 0; dd < d; ++dd) {
                    double mu, vv;
                    if (den > 1e-300) {
                        mu = sumx[idx * d + dd] / den;
                        vv = sumxx[idx * d + dd] / den - mu * mu;
                    } else {
                        mu = cur->mean[idx * d + dd];
                        vv = cur->var[idx * d + dd];
                    }
                    if (!isfinite(mu)) mu = cur->mean[idx * d + dd];
                    if (!isfinite(vv) || vv < opts->var_floor)
                        vv = opts->var_floor;
                    mean_new[idx * d + dd] = mu;
                    var_new[idx * d + dd] = vv;
                }
            }
        }
        memcpy(cur->a, a_new, n * n * sizeof(double));
        memcpy(cur->mix, mix_new, n * k * sizeof(double));
        memcpy(cur->mean, mean_new, n * k * d * sizeof(double));
        memcpy(cur->var, var_new, n * k * d * sizeof(double));

        compute_tables(cur, x, T, logb, logc);
        double ll_new = forward_from_logb(cur, logb, T, la);
        if (!isfinite(ll_new) || ll_new < best_ll - 1e-9) {
            ghmm_copy_params(cur, best);
            st = HP_OK;
            break;
        }
        if (ll_new > best_ll + 1e-12) {
            best_ll = ll_new;
            ghmm_copy_params(best, cur);
        }
        double improvement = ll_new - ll;
        ll = ll_new;
        if (improvement >= 0.0 &&
            improvement < opts->tol * (1.0 + fabs(ll_new))) {
            conv = 1;
            break;
        }
    }
    if (it > opts->max_iter) it = opts->max_iter;
    ghmm_canonicalize(best);
    hp_status s2 = hp_ghmm_resize(result, n, d, k);
    if (s2 == HP_OK) {
        ghmm_copy_params(result, best);
    } else {
        st = s2;
    }
    *iterations = it;
    *converged = conv;
    *loglik_out = best_ll;
done:
    hp_ghmm_free(cur);
    hp_ghmm_free(best);
    hp_free(logb);
    hp_free(logc);
    hp_free(la);
    hp_free(lb);
    hp_free(gamma);
    hp_free(xi);
    hp_free(gsum);
    hp_free(gsumk);
    hp_free(sumx);
    hp_free(sumxx);
    hp_free(a_new);
    hp_free(mix_new);
    hp_free(mean_new);
    hp_free(var_new);
    return st;
}

hp_status hp_ghmm_train(const double *x, size_t T, size_t dim,
                        size_t n_states, size_t n_mix,
                        const hp_ghmm_train_opts *opts, hp_ghmm *out,
                        hp_ghmm_report *report) {
    if (!x || !out || T < 2 || dim == 0 || n_states == 0 || n_mix == 0)
        return HP_ERR_INVALID;
    if (n_states > T) return HP_ERR_INVALID;
    if (!x_valid(x, T, dim)) return HP_ERR_DOMAIN;
    hp_ghmm_train_opts def;
    if (!opts) {
        hp_ghmm_train_opts_default(&def);
        opts = &def;
    }
    if (opts->max_iter == 0 || !(opts->tol > 0.0) || !(opts->alpha > 0.0) ||
        !(opts->var_floor > 0.0) || !(opts->mix_floor > 0.0))
        return HP_ERR_INVALID;
    unsigned restarts = opts->n_restarts ? opts->n_restarts : 1;

    double *gmean = hp_doubles(dim), *gvar = hp_doubles(dim);
    for (size_t dd = 0; dd < dim; ++dd) {
        double s = 0.0, s2 = 0.0;
        for (size_t t = 0; t < T; ++t) {
            double v = x[t * dim + dd];
            s += v;
            s2 += v * v;
        }
        double mu = s / (double)T;
        gmean[dd] = mu;
        double vv = s2 / (double)T - mu * mu;
        gvar[dd] = vv > 1e-6 ? vv : 1.0;
    }

    hp_ghmm *init = hp_ghmm_new(n_states, dim, n_mix);
    hp_ghmm *best = hp_ghmm_new(n_states, dim, n_mix);
    hp_ghmm *run = hp_ghmm_new(n_states, dim, n_mix);
    double best_ll = -INFINITY;
    unsigned best_it = 0;
    int best_conv = 0;
    hp_status st = HP_OK;

    for (unsigned r = 0; r < restarts; ++r) {
        if (r == 0 && opts->init == HP_GHMM_INIT_QUANTILE) {
            init_quantile(init, x, T, opts->var_floor);
        } else {
            hp_rng rng;
            hp_rng_seed_stream(&rng, opts->seed, 0x2000u + r);
            hp_ghmm_randomize(init, &rng, opts->alpha, gvar);
            for (size_t dd = 0; dd < dim; ++dd) {
                /* Centre random means on the data instead of on zero. */
                for (size_t j = 0; j < n_states * n_mix; ++j)
                    init->mean[j * dim + dd] += gmean[dd];
            }
        }
        unsigned it = 0;
        int conv = 0;
        double ll = -INFINITY;
        hp_status s = ghmm_em(x, T, init, opts, run, &it, &conv, &ll);
        if (s != HP_OK) {
            st = s;
            continue;
        }
        if (ll > best_ll) {
            best_ll = ll;
            best_it = it;
            best_conv = conv;
            ghmm_copy_params(best, run);
        }
        if (opts->verbose)
            HP_LOGI("ghmm n=%zu k=%zu restart %u: loglik=%.10f it=%u conv=%d",
                    n_states, n_mix, r, ll, it, conv);
    }
    if (!isfinite(best_ll)) {
        hp_ghmm_free(init);
        hp_ghmm_free(best);
        hp_ghmm_free(run);
        hp_free(gmean);
        hp_free(gvar);
        return st != HP_OK ? st : HP_ERR_DOMAIN;
    }
    hp_status s2 = hp_ghmm_resize(out, n_states, dim, n_mix);
    if (s2 != HP_OK) {
        hp_ghmm_free(init);
        hp_ghmm_free(best);
        hp_ghmm_free(run);
        hp_free(gmean);
        hp_free(gvar);
        return s2;
    }
    ghmm_copy_params(out, best);
    if (report)
        hp_ghmm_report_fill(out, T, best_ll, best_it, best_conv, restarts, report);
    hp_ghmm_free(init);
    hp_ghmm_free(best);
    hp_ghmm_free(run);
    hp_free(gmean);
    hp_free(gvar);
    return HP_OK;
}

hp_status hp_ghmm_select(const double *x, size_t T, size_t dim, size_t n_min,
                         size_t n_max, size_t n_mix,
                         const hp_ghmm_train_opts *opts, hp_ghmm *best,
                         size_t *best_n, hp_ghmm_model_info *table) {
    if (!x || !best || !table || n_min < 1 || n_max < n_min || n_mix == 0)
        return HP_ERR_INVALID;
    hp_ghmm *cur = hp_ghmm_new(n_min, dim, n_mix);
    double best_bic = INFINITY;
    size_t best_states = 0;
    hp_status st = HP_OK;
    for (size_t n = n_min; n <= n_max; ++n) {
        hp_ghmm_report rep;
        memset(&rep, 0, sizeof(rep));
        hp_status s = hp_ghmm_train(x, T, dim, n, n_mix, opts, cur, &rep);
        hp_ghmm_model_info *info = &table[n - n_min];
        info->n_states = n;
        info->loglik = rep.loglik;
        info->bic = rep.bic;
        info->aic = rep.aic;
        info->n_params = rep.n_params;
        info->converged = rep.converged;
        if (s != HP_OK) {
            st = s;
            continue;
        }
        if (rep.bic < best_bic) {
            hp_status s2 = hp_ghmm_resize(best, n, dim, n_mix);
            if (s2 != HP_OK) {
                hp_ghmm_free(cur);
                return s2;
            }
            best_bic = rep.bic;
            best_states = n;
            ghmm_copy_params(best, cur);
        }
    }
    hp_ghmm_free(cur);
    if (best_states == 0) return st != HP_OK ? st : HP_ERR_DOMAIN;
    if (best_n) *best_n = best_states;
    return HP_OK;
}

/* ---------- persistence ---------- */

static int ghmm_next_token(char **p, char **tok) {
    while (**p && isspace((unsigned char)**p)) ++(*p);
    if (!**p) return 0;
    *tok = *p;
    while (**p && !isspace((unsigned char)**p)) ++(*p);
    if (**p) *(*p)++ = '\0';
    return 1;
}

static int ghmm_expect(char **p, const char *want) {
    char *tok = NULL;
    if (!ghmm_next_token(p, &tok)) return 0;
    return strcmp(tok, want) == 0;
}

static int ghmm_read_size(char **p, size_t *out) {
    char *tok = NULL;
    if (!ghmm_next_token(p, &tok)) return 0;
    char *end = NULL;
    unsigned long long v = strtoull(tok, &end, 10);
    if (end == tok || (end && *end)) return 0;
    *out = (size_t)v;
    return 1;
}

static int ghmm_read_double(char **p, double *out) {
    char *tok = NULL;
    if (!ghmm_next_token(p, &tok)) return 0;
    char *end = NULL;
    double v = strtod(tok, &end);
    if (end == tok || (end && *end)) return 0;
    *out = v;
    return 1;
}

static void ghmm_write_array(FILE *f, const double *v, size_t n) {
    for (size_t i = 0; i < n; ++i)
        fprintf(f, "%.17g%c", v[i], (i % 10 == 9 || i + 1 == n) ? '\n' : ' ');
}

hp_status hp_ghmm_save(const hp_ghmm *m, const char *path) {
    if (!m || !path) return HP_ERR_INVALID;
    HP_TRY(hp_ghmm_validate(m));
    FILE *f = fopen(path, "wb");
    if (!f) return HP_ERR_IO;
    fprintf(f, "hp-ghmm 1\nstates %zu\ndim %zu\nnmix %zu\npi\n", m->n_states,
            m->dim, m->n_mix);
    ghmm_write_array(f, m->pi, m->n_states);
    fprintf(f, "a\n");
    ghmm_write_array(f, m->a, m->n_states * m->n_states);
    fprintf(f, "mix\n");
    ghmm_write_array(f, m->mix, m->n_states * m->n_mix);
    fprintf(f, "mean\n");
    ghmm_write_array(f, m->mean, m->n_states * m->n_mix * m->dim);
    fprintf(f, "var\n");
    ghmm_write_array(f, m->var, m->n_states * m->n_mix * m->dim);
    int bad = ferror(f);
    if (fclose(f) != 0 || bad) return HP_ERR_IO;
    return HP_OK;
}

hp_ghmm *hp_ghmm_load(const char *path) {
    if (!path) return NULL;
    char *buf = NULL;
    if (hp_read_file(path, &buf, NULL) != HP_OK) return NULL;
    char *p = buf;
    size_t n = 0, d = 0, k = 0, version = 0;
    hp_ghmm *m = NULL;
    int ok = 1;
    ok = ghmm_expect(&p, "hp-ghmm");
    ok = ok && ghmm_read_size(&p, &version) && version == 1;
    ok = ok && ghmm_expect(&p, "states") && ghmm_read_size(&p, &n);
    ok = ok && ghmm_expect(&p, "dim") && ghmm_read_size(&p, &d);
    ok = ok && ghmm_expect(&p, "nmix") && ghmm_read_size(&p, &k);
    if (!ok || n == 0 || d == 0 || k == 0) goto fail;
    m = hp_ghmm_new(n, d, k);
    ok = ghmm_expect(&p, "pi");
    for (size_t i = 0; ok && i < n; ++i) ok = ghmm_read_double(&p, &m->pi[i]);
    ok = ok && ghmm_expect(&p, "a");
    for (size_t i = 0; ok && i < n * n; ++i) ok = ghmm_read_double(&p, &m->a[i]);
    ok = ok && ghmm_expect(&p, "mix");
    for (size_t i = 0; ok && i < n * k; ++i) ok = ghmm_read_double(&p, &m->mix[i]);
    ok = ok && ghmm_expect(&p, "mean");
    for (size_t i = 0; ok && i < n * k * d; ++i)
        ok = ghmm_read_double(&p, &m->mean[i]);
    ok = ok && ghmm_expect(&p, "var");
    for (size_t i = 0; ok && i < n * k * d; ++i)
        ok = ghmm_read_double(&p, &m->var[i]);
    if (!ok || hp_ghmm_validate(m) != HP_OK) goto fail;
    hp_free(buf);
    return m;
fail:
    hp_ghmm_free(m);
    hp_free(buf);
    return NULL;
}
