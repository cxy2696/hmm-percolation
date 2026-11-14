/* HP · discrete hidden Markov model.
 *
 * Conventions
 *   obs[t] in [0, n_symbols)
 *   pi[i] = P(q_1 = i)
 *   a[i*n + j] = P(q_{t+1} = j | q_t = i)
 *   b[j*m + k] = P(o_t = k | q_t = j)
 * All inference runs in the log domain (log-sum-exp), so long sequences do
 * not underflow; there is no separate scaled implementation to drift out of
 * sync.
 *
 * Ownership: constructors return heap objects, free with hp_hmm_free.  The
 * buffers passed to inference/training are caller-owned.
 */
#ifndef HP_HMM_HMM_H
#define HP_HMM_HMM_H

#include <stddef.h>
#include <stdint.h>

#include "hp/core/rng.h"
#include "hp/core/status.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t n_states;
    size_t n_symbols;
    double *pi; /* n */
    double *a;  /* n*n */
    double *b;  /* n*m */
} hp_hmm;

hp_hmm *hp_hmm_new(size_t n_states, size_t n_symbols);
hp_hmm *hp_hmm_clone(const hp_hmm *hmm);
void hp_hmm_free(hp_hmm *hmm);
/* Reallocates parameter arrays and resets to the uniform model. */
hp_status hp_hmm_resize(hp_hmm *hmm, size_t n_states, size_t n_symbols);

/* Checks sizes, finiteness, non-negativity and row sums (tolerance 1e-6). */
hp_status hp_hmm_validate(const hp_hmm *hmm);
hp_status hp_hmm_set_uniform(hp_hmm *hmm);
hp_status hp_hmm_randomize(hp_hmm *hmm, hp_rng *rng, double alpha);

/* log P(o_1..o_T, q_t=j) in logalpha[t*n+j]; returns log P(O). */
hp_status hp_hmm_forward_log(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                             double *logalpha, double *loglik);
hp_status hp_hmm_backward_log(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                              double *logbeta);
/* gamma[t*n+i] = P(q_t=i | O).  Normalised per time step. */
hp_status hp_hmm_posterior(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                           double *gamma, double *loglik);
/* Causal filtered probabilities P(q_t=i | o_1..o_t). */
hp_status hp_hmm_filter(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                        double *filt, double *loglik);
/* Viterbi path and log joint probability. */
hp_status hp_hmm_viterbi(const hp_hmm *hmm, const uint32_t *obs, size_t T,
                         uint32_t *path, double *logp);

size_t hp_hmm_n_params(size_t n_states, size_t n_symbols);

typedef enum {
    HP_HMM_INIT_QUANTILE = 0,
    HP_HMM_INIT_RANDOM = 1
} hp_hmm_init;

typedef struct {
    unsigned max_iter;
    double tol;
    unsigned n_restarts;
    uint64_t seed;
    double alpha;      /* Dirichlet concentration for random restarts */
    hp_hmm_init init;  /* restart 0 uses this; later restarts are random */
    int verbose;
} hp_hmm_train_opts;

void hp_hmm_train_opts_default(hp_hmm_train_opts *opts);

typedef struct {
    double loglik;
    double bic;
    double aic;
    size_t n_params;
    unsigned iterations;
    int converged;
    unsigned restarts;
} hp_hmm_report;

/* EM (Baum-Welch) with multiple restarts; keeps the best log-likelihood and
 * relabels states by increasing expected symbol so results are comparable
 * across runs. */
hp_status hp_hmm_train(const uint32_t *obs, size_t T, size_t n_states,
                       size_t n_symbols, const hp_hmm_train_opts *opts,
                       hp_hmm *out, hp_hmm_report *report);

hp_status hp_hmm_report_fill(const hp_hmm *hmm, size_t T, double loglik,
                             unsigned iterations, int converged,
                             unsigned restarts, hp_hmm_report *report);

typedef struct {
    size_t n_states;
    double loglik;
    double bic;
    double aic;
    size_t n_params;
    int converged;
} hp_hmm_model_info;

/* Fits n_states = n_min..n_max and returns the BIC-minimising model.
 * table must hold (n_max-n_min+1) entries. */
hp_status hp_hmm_select(const uint32_t *obs, size_t T, size_t n_symbols,
                        size_t n_min, size_t n_max,
                        const hp_hmm_train_opts *opts, hp_hmm *best,
                        size_t *best_n, hp_hmm_model_info *table);

hp_status hp_hmm_save(const hp_hmm *hmm, const char *path);
hp_hmm *hp_hmm_load(const char *path);

#ifdef __cplusplus
}
#endif

#endif /* HP_HMM_HMM_H */
