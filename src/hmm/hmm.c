#include "hp/hmm/hmm.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hp/core/io.h"
#include "hp/core/log.h"
#include "hp/core/memory.h"
#include "hp/core/stats.h"

#define HP_HMM_FLOOR 1e-12

static double safe_log(double x) { return x > 0.0 ? log(x) : -INFINITY; }

static int obs_valid(const uint32_t *obs, size_t T, size_t m) {
    if (!obs) return 0;
    for (size_t t = 0; t < T; ++t)
        if ((size_t)obs[t] >= m) return 0;
    return 1;
}

hp_hmm *hp_hmm_new(size_t n_states, size_t n_symbols) {
    if (n_states == 0 || n_symbols == 0) return NULL;
    hp_hmm *h = (hp_hmm *)hp_xcalloc(1, sizeof(*h));
    h->n_states = n_states;
    h->n_symbols = n_symbols;
    h->pi = hp_doubles(n_states);
    h->a = hp_doubles(n_states * n_states);
    h->b = hp_doubles(n_states * n_symbols);
    hp_hmm_set_uniform(h);
    return h;
}

hp_hmm *hp_hmm_clone(const hp_hmm *hmm) {
    if (!hmm) return NULL;
    hp_hmm *h = hp_hmm_new(hmm->n_states, hmm->n_symbols);
    if (!h) return NULL;
    memcpy(h->pi, hmm->pi, h->n_states * sizeof(double));
    memcpy(h->a, hmm->a, h->n_states * h->n_states * sizeof(double));
    memcpy(h->b, hmm->b, h->n_states * h->n_symbols * sizeof(double));
    return h;
}

void hp_hmm_free(hp_hmm *hmm) {
    if (!hmm) return;
    hp_free(hmm->pi);
    hp_free(hmm->a);
    hp_free(hmm->b);
    hp_free(hmm);
}

hp_status hp_hmm_resize(hp_hmm *hmm, size_t n_states, size_t n_symbols) {
    if (!hmm || n_states == 0 || n_symbols == 0) return HP_ERR_INVALID;
    hp_free(hmm->pi);
    hp_free(hmm->a);
    hp_free(hmm->b);
    hmm->n_states = n_states;
    hmm->n_symbols = n_symbols;
    hmm->pi = hp_doubles(n_states);
    hmm->a = hp_doubles(n_states * n_states);
    hmm->b = hp_doubles(n_states * n_symbols);
    hp_hmm_set_uniform(hmm);
    return HP_OK;
}

hp_status hp_hmm_validate(const hp_hmm *hmm) {
    if (!hmm || hmm->n_states == 0 || hmm->n_symbols == 0) return HP_ERR_INVALID;
    if (!hmm->pi || !hmm->a || !hmm->b) return HP_ERR_INVALID;
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double v = hmm->pi[i];
        if (!isfinite(v) || v < -1e-12 || v > 1.0 + 1e-9) return HP_ERR_DOMAIN;
        s += v;
    }
    if (fabs(s - 1.0) > 1e-6) return HP_ERR_DOMAIN;
    for (size_t i = 0; i < n; ++i) {
        s = 0.0;
        for (size_t j = 0; j < n; ++j) {
            double v = hmm->a[i * n + j];
            if (!isfinite(v) || v < -1e-12 || v > 1.0 + 1e-9) return HP_ERR_DOMAIN;
            s += v;
        }
        if (fabs(s - 1.0) > 1e-6) return HP_ERR_DOMAIN;
    }
    for (size_t j = 0; j < n; ++j) {
        s = 0.0;
        for (size_t k = 0; k < m; ++k) {
            double v = hmm->b[j * m + k];
            if (!isfinite(v) || v < -1e-12 || v > 1.0 + 1e-9) return HP_ERR_DOMAIN;
            s += v;
        }
        if (fabs(s - 1.0) > 1e-6) return HP_ERR_DOMAIN;
    }
    return HP_OK;
}

hp_status hp_hmm_set_uniform(hp_hmm *hmm) {
    if (!hmm || !hmm->pi || !hmm->a || !hmm->b) return HP_ERR_INVALID;
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    for (size_t i = 0; i < n; ++i) {
        hmm->pi[i] = 1.0 / (double)n;
        for (size_t j = 0; j < n; ++j) hmm->a[i * n + j] = 1.0 / (double)n;
        for (size_t k = 0; k < m; ++k) hmm->b[i * m + k] = 1.0 / (double)m;
    }
    return HP_OK;
}

hp_status hp_hmm_randomize(hp_hmm *hmm, hp_rng *rng, double alpha) {
    if (!hmm || !rng || !(alpha > 0.0)) return HP_ERR_INVALID;
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    const size_t nmax = n > m ? n : m;
    double *a = hp_doubles(nmax);
    for (size_t i = 0; i < n; ++i) a[i] = alpha;
    hp_rng_dirichlet(rng, a, n, hmm->pi);
    for (size_t i = 0; i < n; ++i) hp_rng_dirichlet(rng, a, n, hmm->a + i * n);
    for (size_t j = 0; j < n; ++j) {
        for (size_t k = 0; k < m; ++k) a[k] = alpha;
        hp_rng_dirichlet(rng, a, m, hmm->b + j * m);
    }
    hp_free(a);
    return HP_OK;
}

/* Fills log(pi), log(A), log(B). */
static void hmm_build_logs(const hp_hmm *hmm, double *loga, double *logb) {
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    for (size_t i = 0; i < n * n; ++i) loga[i] = safe_log(hmm->a[i]);
    for (size_t i = 0; i < n * m; ++i) logb[i] = safe_log(hmm->b[i]);
}

hp_status hp_hmm_forward_log(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                             double *logalpha, double *loglik) {
    if (!hmm || !obs || !logalpha || !loglik || T == 0) return HP_ERR_INVALID;
    HP_TRY(hp_hmm_validate(hmm));
    if (!obs_valid(obs, T, hmm->n_symbols)) return HP_ERR_RANGE;
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    const double *pi = hmm->pi, *A = hmm->a, *B = hmm->b;

    for (size_t j = 0; j < n; ++j) {
        double lpj = safe_log(pi[j]);
        double lbj = safe_log(B[j * m + obs[0]]);
        logalpha[j] = lpj + lbj;
    }
    for (size_t t = 1; t < T; ++t) {
        const double *prev = logalpha + (t - 1) * n;
        double *cur = logalpha + t * n;
        for (size_t j = 0; j < n; ++j) {
            double acc = -INFINITY;
            for (size_t i = 0; i < n; ++i) {
                double aij = A[i * n + j];
                if (aij <= 0.0) continue;
                double v = prev[i] + log(aij);
                if (v > acc) {
                    if (acc == -INFINITY) {
                        acc = v;
                    } else {
                        acc = acc + log1p(exp(v - acc));
                    }
                } else if (v > -INFINITY) {
                    acc = acc + log1p(exp(v - acc));
                }
            }
            cur[j] = acc + safe_log(B[j * m + obs[t]]);
        }
    }
    *loglik = hp_logsumexp(logalpha + (T - 1) * n, n);
    return isfinite(*loglik) ? HP_OK : HP_ERR_DOMAIN;
}

hp_status hp_hmm_backward_log(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                              double *logbeta) {
    if (!hmm || !obs || !logbeta || T == 0) return HP_ERR_INVALID;
    HP_TRY(hp_hmm_validate(hmm));
    if (!obs_valid(obs, T, hmm->n_symbols)) return HP_ERR_RANGE;
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    const double *A = hmm->a, *B = hmm->b;
    for (size_t i = 0; i < n; ++i) logbeta[(T - 1) * n + i] = 0.0;
    for (size_t t = T - 1; t-- > 0;) {
        double *cur = logbeta + t * n;
        const double *next = logbeta + (t + 1) * n;
        uint32_t o = obs[t + 1];
        for (size_t i = 0; i < n; ++i) {
            double acc = -INFINITY;
            for (size_t j = 0; j < n; ++j) {
                double aij = A[i * n + j];
                double bj = B[j * m + o];
                if (aij <= 0.0 || bj <= 0.0) continue;
                double v = log(aij) + log(bj) + next[j];
                if (v > acc) {
                    if (acc == -INFINITY) {
                        acc = v;
                    } else {
                        acc = acc + log1p(exp(v - acc));
                    }
                } else if (v > -INFINITY) {
                    acc = acc + log1p(exp(v - acc));
                }
            }
            cur[i] = acc;
        }
    }
    return HP_OK;
}

hp_status hp_hmm_posterior(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                           double *gamma, double *loglik) {
    if (!gamma) return HP_ERR_INVALID;
    const size_t n = hmm ? hmm->n_states : 0;
    double *la = hp_doubles(T * n);
    double *lb = hp_doubles(T * n);
    double ll = 0.0;
    hp_status st = hp_hmm_forward_log(hmm, obs, T, la, &ll);
    if (st == HP_OK) st = hp_hmm_backward_log(hmm, obs, T, lb);
    if (st == HP_OK) {
        for (size_t t = 0; t < T; ++t) {
            double *g = gamma + t * n;
            double s = 0.0;
            for (size_t i = 0; i < n; ++i) {
                g[i] = exp(la[t * n + i] + lb[t * n + i] - ll);
                if (!isfinite(g[i]) || g[i] < 0.0) g[i] = 0.0;
                s += g[i];
            }
            if (s <= 0.0) {
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

hp_status hp_hmm_filter(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                        double *filt, double *loglik) {
    if (!filt) return HP_ERR_INVALID;
    const size_t n = hmm ? hmm->n_states : 0;
    double *la = hp_doubles(T * n);
    double ll = 0.0;
    hp_status st = hp_hmm_forward_log(hmm, obs, T, la, &ll);
    if (st == HP_OK) {
        for (size_t t = 0; t < T; ++t) {
            const double *row = la + t * n;
            double lse = hp_logsumexp(row, n);
            if (!isfinite(lse)) {
                st = HP_ERR_DOMAIN;
                break;
            }
            for (size_t i = 0; i < n; ++i) filt[t * n + i] = exp(row[i] - lse);
        }
    }
    if (loglik) *loglik = ll;
    hp_free(la);
    return st;
}

hp_status hp_hmm_viterbi(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                         uint32_t *path, double *logp) {
    if (!hmm || !obs || !path || T == 0) return HP_ERR_INVALID;
    HP_TRY(hp_hmm_validate(hmm));
    if (!obs_valid(obs, T, hmm->n_symbols)) return HP_ERR_RANGE;
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    const double *pi = hmm->pi, *A = hmm->a, *B = hmm->b;
    double *delta = hp_doubles(n);
    double *next = hp_doubles(n);
    uint32_t *psi = hp_u32s(T * n);

    for (size_t j = 0; j < n; ++j)
        delta[j] = safe_log(pi[j]) + safe_log(B[j * m + obs[0]]);
    for (size_t t = 1; t < T; ++t) {
        for (size_t j = 0; j < n; ++j) {
            double best = -INFINITY;
            uint32_t arg = 0;
            for (size_t i = 0; i < n; ++i) {
                double aij = A[i * n + j];
                if (aij <= 0.0 || delta[i] == -INFINITY) continue;
                double v = delta[i] + log(aij);
                if (v > best) {
                    best = v;
                    arg = (uint32_t)i;
                }
            }
            next[j] = best + safe_log(B[j * m + obs[t]]);
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
        hp_free(delta);
        hp_free(next);
        hp_free(psi);
        return HP_ERR_DOMAIN;
    }
    path[T - 1] = arg;
    for (size_t t = T - 1; t-- > 0;) path[t] = psi[(t + 1) * n + path[t + 1]];
    if (logp) *logp = best;
    hp_free(delta);
    hp_free(next);
    hp_free(psi);
    return HP_OK;
}

size_t hp_hmm_n_params(size_t n_states, size_t n_symbols) {
    if (n_states == 0 || n_symbols == 0) return 0;
    return (n_states - 1) + n_states * (n_states - 1) +
           n_states * (n_symbols - 1);
}

void hp_hmm_train_opts_default(hp_hmm_train_opts *opts) {
    if (!opts) return;
    opts->max_iter = 200;
    opts->tol = 1e-6;
    opts->n_restarts = 4;
    opts->seed = 42;
    opts->alpha = 1.0;
    opts->init = HP_HMM_INIT_QUANTILE;
    opts->verbose = 0;
}

static void hmm_init_quantile(hp_hmm *hmm, const uint32_t *obs, size_t T) {
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    size_t *cnt = hp_sizes(m);
    for (size_t t = 0; t < T; ++t) cnt[obs[t]]++;
    for (size_t i = 0; i < n; ++i) hmm->pi[i] = 1.0 / (double)n;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            hmm->a[i * n + j] = (i == j) ? 0.9 : (n > 1 ? 0.1 / (double)(n - 1) : 0.0);
        if (n == 1) hmm->a[i] = 1.0;
    }
    double acc = 0.0;
    double *mass = hp_doubles(m);
    for (size_t k = 0; k < m; ++k) mass[k] = (double)cnt[k] / (double)T;
    for (size_t k = 0; k < m; ++k) {
        double mid = acc + 0.5 * mass[k];
        acc += mass[k];
        size_t s = (size_t)(mid * (double)n);
        if (s >= n) s = n - 1;
        for (size_t j = 0; j < n; ++j) hmm->b[j * m + k] = (j == s) ? 1.0 : 0.01;
    }
    for (size_t j = 0; j < n; ++j) {
        double s = 0.0;
        for (size_t k = 0; k < m; ++k) s += hmm->b[j * m + k];
        for (size_t k = 0; k < m; ++k) hmm->b[j * m + k] /= s;
    }
    hp_free(cnt);
    hp_free(mass);
}

static int cmp_key(const double *key, size_t a, size_t b) {
    if (key[a] < key[b]) return -1;
    if (key[a] > key[b]) return 1;
    return (a > b) - (a < b);
}

/* Relabel states by increasing expected symbol. */
static void hmm_canonicalize(hp_hmm *hmm) {
    const size_t n = hmm->n_states, m = hmm->n_symbols;
    double *key = hp_doubles(n);
    size_t *order = hp_sizes(n);
    for (size_t j = 0; j < n; ++j) {
        double e = 0.0;
        for (size_t k = 0; k < m; ++k) e += (double)k * hmm->b[j * m + k];
        key[j] = e;
        order[j] = j;
    }
    for (size_t i = 1; i < n; ++i) {
        size_t v = order[i];
        size_t j = i;
        while (j > 0 && cmp_key(key, v, order[j - 1]) < 0) {
            order[j] = order[j - 1];
            --j;
        }
        order[j] = v;
    }
    double *pi = hp_doubles(n), *a = hp_doubles(n * n), *b = hp_doubles(n * m);
    for (size_t j = 0; j < n; ++j) {
        pi[j] = hmm->pi[order[j]];
        for (size_t k = 0; k < m; ++k) b[j * m + k] = hmm->b[order[j] * m + k];
    }
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            a[i * n + j] = hmm->a[order[i] * n + order[j]];
    memcpy(hmm->pi, pi, n * sizeof(double));
    memcpy(hmm->a, a, n * n * sizeof(double));
    memcpy(hmm->b, b, n * m * sizeof(double));
    hp_free(key);
    hp_free(order);
    hp_free(pi);
    hp_free(a);
    hp_free(b);
}

static hp_status hmm_em(const uint32_t *obs, size_t T, const hp_hmm *init,
                        const hp_hmm_train_opts *opts, hp_hmm *out,
                        unsigned *iterations, int *converged, double *loglik) {
    const size_t n = init->n_states, m = init->n_symbols;
    hp_hmm *cur = hp_hmm_clone(init);
    hp_hmm *best = hp_hmm_clone(init);
    double *la = hp_doubles(T * n);
    double *lb = hp_doubles(T * n);
    double *gamma = hp_doubles(T * n);
    double *xi = hp_doubles(n * n);
    double *gsum = hp_doubles(n);
    double *logA = hp_doubles(n * n);
    double *logB = hp_doubles(n * m);
    double *a_new = hp_doubles(n * n);
    double *b_new = hp_doubles(n * m);

    double ll = 0.0;
    hp_status st = hp_hmm_forward_log(cur, obs, T, la, &ll);
    if (st != HP_OK) goto done;
    double best_ll = ll;
    unsigned it = 0;
    int conv = 0;

    for (it = 1; it <= opts->max_iter; ++it) {
        st = hp_hmm_backward_log(cur, obs, T, lb);
        if (st != HP_OK) break;

        hmm_build_logs(cur, logA, logB);
        memset(xi, 0, n * n * sizeof(double));
        memset(gsum, 0, n * sizeof(double));

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
                goto em_break;
            }
            for (size_t i = 0; i < n; ++i) {
                g[i] /= s;
                gsum[i] += g[i];
            }
        }

        /* xi needs the log alpha/beta, not gamma, to stay accurate when a
         * transition is very unlikely. */
        for (size_t t = 0; t + 1 < T; ++t) {
            uint32_t o = obs[t + 1];
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    double lbj = logB[j * m + o];
                    if (lbj == -INFINITY) continue;
                    double lv = la[t * n + i] + logA[i * n + j] + lbj +
                                lb[(t + 1) * n + j] - ll;
                    if (lv > -745.0) xi[i * n + j] += exp(lv);
                }
            }
        }

        /* pi */
        for (size_t i = 0; i < n; ++i) cur->pi[i] = gamma[i];
        double pis = 0.0;
        for (size_t i = 0; i < n; ++i) pis += cur->pi[i];
        if (pis > 0.0)
            for (size_t i = 0; i < n; ++i) cur->pi[i] /= pis;
        else
            for (size_t i = 0; i < n; ++i) cur->pi[i] = 1.0 / (double)n;

        /* A: denominator excludes the last time step. */
        for (size_t i = 0; i < n; ++i) {
            double den = gsum[i] - gamma[(T - 1) * n + i];
            double row_sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                double v;
                if (den > 1e-300)
                    v = xi[i * n + j] / den;
                else
                    v = cur->a[i * n + j];
                if (!(v >= 0.0) || !isfinite(v)) v = 0.0;
                if (v < HP_HMM_FLOOR) v = HP_HMM_FLOOR;
                a_new[i * n + j] = v;
                row_sum += v;
            }
            for (size_t j = 0; j < n; ++j) a_new[i * n + j] /= row_sum;
        }

        /* B */
        memset(b_new, 0, n * m * sizeof(double));
        for (size_t t = 0; t < T; ++t) {
            uint32_t o = obs[t];
            for (size_t j = 0; j < n; ++j)
                b_new[j * m + o] += gamma[t * n + j];
        }
        for (size_t j = 0; j < n; ++j) {
            double den = gsum[j];
            double row_sum = 0.0;
            for (size_t k = 0; k < m; ++k) {
                double v;
                if (den > 1e-300)
                    v = b_new[j * m + k] / den;
                else
                    v = cur->b[j * m + k];
                if (!(v >= 0.0) || !isfinite(v)) v = 0.0;
                if (v < HP_HMM_FLOOR) v = HP_HMM_FLOOR;
                b_new[j * m + k] = v;
                row_sum += v;
            }
            for (size_t k = 0; k < m; ++k) b_new[j * m + k] /= row_sum;
        }
        memcpy(cur->a, a_new, n * n * sizeof(double));
        memcpy(cur->b, b_new, n * m * sizeof(double));

        double ll_new = 0.0;
        st = hp_hmm_forward_log(cur, obs, T, la, &ll_new);
        if (st != HP_OK) break;
        if (!isfinite(ll_new) || ll_new < best_ll - 1e-9) {
            /* Numerical regression: keep the best parameters seen. */
            hp_status s2 = hp_hmm_resize(cur, n, m);
            if (s2 != HP_OK) {
                st = s2;
                break;
            }
            memcpy(cur->pi, best->pi, n * sizeof(double));
            memcpy(cur->a, best->a, n * n * sizeof(double));
            memcpy(cur->b, best->b, n * m * sizeof(double));
            st = HP_OK;
            break;
        }
        if (ll_new > best_ll + 1e-12) {
            best_ll = ll_new;
            memcpy(best->pi, cur->pi, n * sizeof(double));
            memcpy(best->a, cur->a, n * n * sizeof(double));
            memcpy(best->b, cur->b, n * m * sizeof(double));
        }
        double improvement = ll_new - ll;
        ll = ll_new;
        if (improvement >= 0.0 &&
            improvement < opts->tol * (1.0 + fabs(ll_new))) {
            conv = 1;
            break;
        }
    }

em_break:
    if (it > opts->max_iter) it = opts->max_iter;
    hmm_canonicalize(best);
    if (st == HP_OK || st == HP_ERR_DOMAIN) {
        hp_status s2 = hp_hmm_resize(out, n, m);
        if (s2 == HP_OK) {
            memcpy(out->pi, best->pi, n * sizeof(double));
            memcpy(out->a, best->a, n * n * sizeof(double));
            memcpy(out->b, best->b, n * m * sizeof(double));
        } else {
            st = s2;
        }
    }
    *iterations = it;
    *converged = conv;
    *loglik = best_ll;
    if (st == HP_ERR_DOMAIN) st = HP_OK;
done:
    hp_hmm_free(cur);
    hp_hmm_free(best);
    hp_free(la);
    hp_free(lb);
    hp_free(gamma);
    hp_free(xi);
    hp_free(gsum);
    hp_free(logA);
    hp_free(logB);
    hp_free(a_new);
    hp_free(b_new);
    return st;
}

hp_status hp_hmm_report_fill(const hp_hmm *hmm, size_t T, double loglik,
                             unsigned iterations, int converged,
                             unsigned restarts, hp_hmm_report *report) {
    if (!hmm || !report || T == 0) return HP_ERR_INVALID;
    size_t k = hp_hmm_n_params(hmm->n_states, hmm->n_symbols);
    report->loglik = loglik;
    report->n_params = k;
    report->bic = -2.0 * loglik + (double)k * log((double)T);
    report->aic = -2.0 * loglik + 2.0 * (double)k;
    report->iterations = iterations;
    report->converged = converged;
    report->restarts = restarts;
    return HP_OK;
}

hp_status hp_hmm_train(const uint32_t *obs, size_t T, size_t n_states,
                       size_t n_symbols, const hp_hmm_train_opts *opts,
                       hp_hmm *out, hp_hmm_report *report) {
    if (!obs || !out || T < 2 || n_states == 0 || n_symbols == 0)
        return HP_ERR_INVALID;
    if (!obs_valid(obs, T, n_symbols)) return HP_ERR_RANGE;
    hp_hmm_train_opts def;
    if (!opts) {
        hp_hmm_train_opts_default(&def);
        opts = &def;
    }
    if (opts->max_iter == 0 || !(opts->tol > 0.0) || !(opts->alpha > 0.0))
        return HP_ERR_INVALID;
    unsigned restarts = opts->n_restarts ? opts->n_restarts : 1;

    hp_hmm *init = hp_hmm_new(n_states, n_symbols);
    hp_hmm *best = hp_hmm_new(n_states, n_symbols);
    hp_hmm *run = hp_hmm_new(n_states, n_symbols);
    double best_ll = -INFINITY;
    unsigned best_it = 0;
    int best_conv = 0;
    hp_status st = HP_OK;

    for (unsigned r = 0; r < restarts; ++r) {
        if (r == 0 && opts->init == HP_HMM_INIT_QUANTILE) {
            hmm_init_quantile(init, obs, T);
        } else {
            hp_rng rng;
            hp_rng_seed_stream(&rng, opts->seed, 0x1000u + r);
            hp_hmm_randomize(init, &rng, opts->alpha);
        }
        unsigned it = 0;
        int conv = 0;
        double ll = -INFINITY;
        hp_status s = hmm_em(obs, T, init, opts, run, &it, &conv, &ll);
        if (s != HP_OK) {
            st = s;
            continue;
        }
        if (ll > best_ll) {
            best_ll = ll;
            best_it = it;
            best_conv = conv;
            memcpy(best->pi, run->pi, n_states * sizeof(double));
            memcpy(best->a, run->a, n_states * n_states * sizeof(double));
            memcpy(best->b, run->b, n_states * n_symbols * sizeof(double));
        }
        if (opts->verbose)
            HP_LOGI("hmm train n=%zu m=%zu restart %u: loglik=%.10f it=%u conv=%d",
                    n_states, n_symbols, r, ll, it, conv);
    }
    if (!isfinite(best_ll)) {
        hp_hmm_free(init);
        hp_hmm_free(best);
        hp_hmm_free(run);
        return st != HP_OK ? st : HP_ERR_DOMAIN;
    }
    hp_status s2 = hp_hmm_resize(out, n_states, n_symbols);
    if (s2 != HP_OK) {
        hp_hmm_free(init);
        hp_hmm_free(best);
        hp_hmm_free(run);
        return s2;
    }
    memcpy(out->pi, best->pi, n_states * sizeof(double));
    memcpy(out->a, best->a, n_states * n_states * sizeof(double));
    memcpy(out->b, best->b, n_states * n_symbols * sizeof(double));
    if (report)
        hp_hmm_report_fill(out, T, best_ll, best_it, best_conv, restarts, report);
    hp_hmm_free(init);
    hp_hmm_free(best);
    hp_hmm_free(run);
    return HP_OK;
}

hp_status hp_hmm_select(const uint32_t *obs, size_t T, size_t n_symbols,
                        size_t n_min, size_t n_max,
                        const hp_hmm_train_opts *opts, hp_hmm *best,
                        size_t *best_n, hp_hmm_model_info *table) {
    if (!obs || !best || !table || n_min < 1 || n_max < n_min)
        return HP_ERR_INVALID;
    if (n_max > T + 8) return HP_ERR_RANGE;
    hp_hmm *cur = hp_hmm_new(n_min, n_symbols);
    double best_bic = INFINITY;
    size_t best_states = 0;
    hp_status st = HP_OK;
    for (size_t n = n_min; n <= n_max; ++n) {
        hp_hmm_report rep;
        memset(&rep, 0, sizeof(rep));
        hp_status s = hp_hmm_train(obs, T, n, n_symbols, opts, cur, &rep);
        hp_hmm_model_info *info = &table[n - n_min];
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
            hp_status s2 = hp_hmm_resize(best, n, n_symbols);
            if (s2 != HP_OK) {
                hp_hmm_free(cur);
                return s2;
            }
            best_bic = rep.bic;
            best_states = n;
            memcpy(best->pi, cur->pi, n * sizeof(double));
            memcpy(best->a, cur->a, n * n * sizeof(double));
            memcpy(best->b, cur->b, n * n_symbols * sizeof(double));
        }
    }
    hp_hmm_free(cur);
    if (best_states == 0) return st != HP_OK ? st : HP_ERR_DOMAIN;
    if (best_n) *best_n = best_states;
    return HP_OK;
}

/* ---------- persistence ---------- */

static int next_token(char **p, char **tok) {
    while (**p && isspace((unsigned char)**p)) ++(*p);
    if (!**p) return 0;
    *tok = *p;
    while (**p && !isspace((unsigned char)**p)) ++(*p);
    if (**p) *(*p)++ = '\0';
    return 1;
}

static int expect_token(char **p, const char *want) {
    char *tok = NULL;
    if (!next_token(p, &tok)) return 0;
    return strcmp(tok, want) == 0;
}

static int read_size_tok(char **p, size_t *out) {
    char *tok = NULL;
    if (!next_token(p, &tok)) return 0;
    char *end = NULL;
    unsigned long long v = strtoull(tok, &end, 10);
    if (end == tok || (end && *end)) return 0;
    *out = (size_t)v;
    return 1;
}

static int read_double_tok(char **p, double *out) {
    char *tok = NULL;
    if (!next_token(p, &tok)) return 0;
    char *end = NULL;
    double v = strtod(tok, &end);
    if (end == tok || (end && *end)) return 0;
    *out = v;
    return 1;
}

static void write_array(FILE *f, const double *v, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        fprintf(f, "%.17g%c", v[i], (i % 10 == 9 || i + 1 == n) ? '\n' : ' ');
    }
}

hp_status hp_hmm_save(const hp_hmm *hmm, const char *path) {
    if (!hmm || !path) return HP_ERR_INVALID;
    HP_TRY(hp_hmm_validate(hmm));
    FILE *f = fopen(path, "wb");
    if (!f) return HP_ERR_IO;
    fprintf(f, "hp-hmm 1\nstates %zu\nsymbols %zu\npi\n", hmm->n_states,
            hmm->n_symbols);
    write_array(f, hmm->pi, hmm->n_states);
    fprintf(f, "a\n");
    write_array(f, hmm->a, hmm->n_states * hmm->n_states);
    fprintf(f, "b\n");
    write_array(f, hmm->b, hmm->n_states * hmm->n_symbols);
    int bad = ferror(f);
    if (fclose(f) != 0 || bad) return HP_ERR_IO;
    return HP_OK;
}

hp_hmm *hp_hmm_load(const char *path) {
    if (!path) return NULL;
    char *buf = NULL;
    if (hp_read_file(path, &buf, NULL) != HP_OK) return NULL;
    char *p = buf;
    size_t n = 0, m = 0;
    hp_hmm *h = NULL;
    int ok = 1;
    ok = expect_token(&p, "hp-hmm");
    size_t version = 0;
    ok = ok && read_size_tok(&p, &version) && version == 1;
    ok = ok && expect_token(&p, "states") && read_size_tok(&p, &n);
    ok = ok && expect_token(&p, "symbols") && read_size_tok(&p, &m);
    if (!ok || n == 0 || m == 0) goto fail;
    h = hp_hmm_new(n, m);
    ok = expect_token(&p, "pi");
    for (size_t i = 0; ok && i < n; ++i) ok = read_double_tok(&p, &h->pi[i]);
    ok = ok && expect_token(&p, "a");
    for (size_t i = 0; ok && i < n * n; ++i) ok = read_double_tok(&p, &h->a[i]);
    ok = ok && expect_token(&p, "b");
    for (size_t i = 0; ok && i < n * m; ++i) ok = read_double_tok(&p, &h->b[i]);
    if (!ok || hp_hmm_validate(h) != HP_OK) goto fail;
    hp_free(buf);
    return h;
fail:
    hp_hmm_free(h);
    hp_free(buf);
    return NULL;
}
