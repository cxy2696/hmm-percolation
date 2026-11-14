/* HP · Gaussian-mixture HMM with diagonal covariance.
 *
 *   emission  b_j(x) = sum_k mix[j,k] * prod_d N(x_d; mean[j,k,d], var[j,k,d])
 *
 * Diagonal covariances are deliberate: market features are few and strongly
 * collinear (return / volatility), and a diagonal mixture stays positive
 * definite without shrinkage or Cholesky repairs.  The recursions are the same
 * as the full-covariance case; only the density changes.
 *
 * All inference is in the log domain.  x is row-major T x dim and must be
 * finite (the app layer trims warm-up rows and rejects holes).
 */
#ifndef HP_HMM_GHMM_H
#define HP_HMM_GHMM_H

#include <stddef.h>
#include <stdint.h>

#include "hp/core/rng.h"
#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t n_states;
    size_t dim;
    size_t n_mix;
    double *pi;   /* n */
    double *a;    /* n*n */
    double *mix;  /* n*k */
    double *mean; /* n*k*d */
    double *var;  /* n*k*d */
} hp_ghmm;

hp_ghmm *hp_ghmm_new(size_t n_states, size_t dim, size_t n_mix);
hp_ghmm *hp_ghmm_clone(const hp_ghmm *m);
void hp_ghmm_free(hp_ghmm *m);
hp_status hp_ghmm_resize(hp_ghmm *m, size_t n_states, size_t dim, size_t n_mix);
hp_status hp_ghmm_validate(const hp_ghmm *m);
hp_status hp_ghmm_set_uniform(hp_ghmm *m);
hp_status hp_ghmm_randomize(hp_ghmm *m, hp_rng *rng, double alpha,
                            const double *global_var);

size_t hp_ghmm_n_params(size_t n_states, size_t dim, size_t n_mix);
/* Mixture-weighted mean of component `dim` for state j. */
double hp_ghmm_state_mean(const hp_ghmm *m, size_t state, size_t dim);

/* log b_j(x) for a single observation. */
hp_status hp_ghmm_log_emission(const hp_ghmm *m, const double *x, double *logb);

hp_status hp_ghmm_forward_log(const hp_ghmm *m, const double *x, size_t T,
                              double *logalpha, double *loglik);
hp_status hp_ghmm_backward_log(const hp_ghmm *m, const double *x, size_t T,
                               double *logbeta);
hp_status hp_ghmm_posterior(const hp_ghmm *m, const double *x, size_t T,
                            double *gamma, double *loglik);
hp_status hp_ghmm_filter(const hp_ghmm *m, const double *x, size_t T,
                         double *filt, double *loglik);
/* One causal filtering step: filt is P(q_t | x_1..x_t) on input and
 * P(q_{t+1} | x_1..x_{t+1}) on output.  loglik_step = log p(x_{t+1} | past). */
hp_status hp_ghmm_filter_step(const hp_ghmm *m, const double *x, double *filt,
                              double *loglik_step);
hp_status hp_ghmm_viterbi(const hp_ghmm *m, const double *x, size_t T,
                          uint32_t *path, double *logp);

typedef enum {
    HP_GHMM_INIT_QUANTILE = 0,
    HP_GHMM_INIT_RANDOM = 1
} hp_ghmm_init;

typedef struct {
    unsigned max_iter;
    double tol;
    unsigned n_restarts;
    uint64_t seed;
    double alpha;      /* Dirichlet concentration, random restarts */
    double var_floor;  /* absolute variance floor, in feature units */
    double mix_floor;
    hp_ghmm_init init;
    int verbose;
} hp_ghmm_train_opts;

void hp_ghmm_train_opts_default(hp_ghmm_train_opts *opts);

typedef struct {
    double loglik;
    double bic;
    double aic;
    size_t n_params;
    unsigned iterations;
    int converged;
    unsigned restarts;
} hp_ghmm_report;

hp_status hp_ghmm_train(const double *x, size_t T, size_t dim,
                        size_t n_states, size_t n_mix,
                        const hp_ghmm_train_opts *opts, hp_ghmm *out,
                        hp_ghmm_report *report);

hp_status hp_ghmm_report_fill(const hp_ghmm *m, size_t T, double loglik,
                              unsigned iterations, int converged,
                              unsigned restarts, hp_ghmm_report *report);

typedef struct {
    size_t n_states;
    double loglik;
    double bic;
    double aic;
    size_t n_params;
    int converged;
} hp_ghmm_model_info;

hp_status hp_ghmm_select(const double *x, size_t T, size_t dim, size_t n_min,
                         size_t n_max, size_t n_mix,
                         const hp_ghmm_train_opts *opts, hp_ghmm *best,
                         size_t *best_n, hp_ghmm_model_info *table);

hp_status hp_ghmm_save(const hp_ghmm *m, const char *path);
hp_ghmm *hp_ghmm_load(const char *path);

#ifdef __cplusplus
}
#endif

#endif /* HP_HMM_GHMM_H */
