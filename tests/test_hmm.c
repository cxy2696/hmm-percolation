#include <math.h>
#include <string.h>

#include "hp/hmm/hmm.h"
#include "hp/core/memory.h"
#include "hp_test.h"

static void make_model(hp_hmm *h) {
    h->n_states = 3;
    h->n_symbols = 3;
    double pi[3] = {0.4, 0.35, 0.25};
    double a[9] = {0.70, 0.20, 0.10,
                   0.15, 0.75, 0.10,
                   0.10, 0.20, 0.70};
    double b[9] = {0.60, 0.30, 0.10,
                   0.10, 0.70, 0.20,
                   0.20, 0.20, 0.60};
    memcpy(h->pi, pi, sizeof(pi));
    memcpy(h->a, a, sizeof(a));
    memcpy(h->b, b, sizeof(b));
}

/* Enumerates all N^T state paths; fills marginals and best path. */
static void brute_force(const hp_hmm *h, const uint32_t *obs, size_t T,
                        double *gamma, double *loglik, uint32_t *best_path,
                        double *best_logp) {
    const size_t n = h->n_states;
    size_t paths = 1;
    for (size_t t = 0; t < T; ++t) paths *= n;
    double total = 0.0;
    double best = -INFINITY;
    size_t best_code = 0;
    if (gamma) memset(gamma, 0, T * n * sizeof(double));
    uint32_t *state = hp_u32s(T);
    for (size_t code = 0; code < paths; ++code) {
        size_t c = code;
        for (size_t t = 0; t < T; ++t) {
            state[t] = (uint32_t)(c % n);
            c /= n;
        }
        double p = h->pi[state[0]] * h->b[state[0] * h->n_symbols + obs[0]];
        for (size_t t = 1; t < T; ++t) {
            p *= h->a[state[t - 1] * n + state[t]];
            p *= h->b[state[t] * h->n_symbols + obs[t]];
        }
        total += p;
        if (gamma)
            for (size_t t = 0; t < T; ++t) gamma[t * n + state[t]] += p;
        if (p > best) {
            best = p;
            best_code = code;
        }
    }
    if (loglik) *loglik = log(total);
    if (best_logp) *best_logp = log(best);
    if (best_path) {
        size_t c = best_code;
        for (size_t t = 0; t < T; ++t) {
            best_path[t] = (uint32_t)(c % n);
            c /= n;
        }
    }
    hp_free(state);
}

HP_TEST(test_discrete_forward_vs_bruteforce) {
    hp_hmm *h = hp_hmm_new(3, 3);
    HP_CHECK(h != NULL);
    make_model(h);
    HP_CHECK(hp_hmm_validate(h) == HP_OK);
    const uint32_t obs[6] = {0, 2, 1, 0, 0, 2};
    const size_t T = 6;
    double logalpha[18];
    double ll = 0;
    HP_CHECK(hp_hmm_forward_log(h, obs, T, logalpha, &ll) == HP_OK);
    double brute_ll = 0;
    brute_force(h, obs, T, NULL, &brute_ll, NULL, NULL);
    HP_CHECK_NEAR(ll, brute_ll, 1e-12);

    double logbeta[18];
    HP_CHECK(hp_hmm_backward_log(h, obs, T, logbeta) == HP_OK);
    double gamma[18];
    double post_ll = 0;
    HP_CHECK(hp_hmm_posterior(h, obs, T, gamma, &post_ll) == HP_OK);
    HP_CHECK_NEAR(post_ll, brute_ll, 1e-12);
    double gamma_brute[18];
    brute_force(h, obs, T, gamma_brute, NULL, NULL, NULL);
    double total = exp(brute_ll);
    for (size_t t = 0; t < T; ++t) {
        double s = 0;
        for (size_t i = 0; i < 3; ++i) {
            HP_CHECK_NEAR(gamma[t * 3 + i], gamma_brute[t * 3 + i] / total, 1e-10);
            s += gamma[t * 3 + i];
        }
        HP_CHECK_NEAR(s, 1.0, 1e-12);
    }
    hp_hmm_free(h);
    return 0;
}

HP_TEST(test_discrete_viterbi_vs_bruteforce) {
    hp_hmm *h = hp_hmm_new(3, 3);
    make_model(h);
    const uint32_t obs[6] = {0, 2, 1, 0, 0, 2};
    uint32_t path[6], bpath[6];
    double lp = 0, blp = 0;
    HP_CHECK(hp_hmm_viterbi(h, obs, 6, path, &lp) == HP_OK);
    brute_force(h, obs, 6, NULL, NULL, bpath, &blp);
    HP_CHECK_NEAR(lp, blp, 1e-12);
    for (size_t t = 0; t < 6; ++t) HP_CHECK(path[t] == bpath[t]);
    hp_hmm_free(h);
    return 0;
}

HP_TEST(test_discrete_filter_vs_prefix) {
    hp_hmm *h = hp_hmm_new(3, 3);
    make_model(h);
    const uint32_t obs[6] = {0, 2, 1, 0, 0, 2};
    double filt[18];
    HP_CHECK(hp_hmm_filter(h, obs, 6, filt, NULL) == HP_OK);
    for (size_t t = 0; t < 6; ++t) {
        double g[3 * 6];
        brute_force(h, obs, t + 1, g, NULL, NULL, NULL);
        double sum = 0;
        for (size_t i = 0; i < 3; ++i) sum += g[t * 3 + i];
        for (size_t i = 0; i < 3; ++i)
            HP_CHECK_NEAR(filt[t * 3 + i], g[t * 3 + i] / sum, 1e-10);
    }
    hp_hmm_free(h);
    return 0;
}

static void sample_two_state(uint32_t *obs, uint32_t *truth, size_t T,
                             uint64_t seed) {
    hp_hmm *h = hp_hmm_new(2, 3);
    double pi[2] = {0.6, 0.4};
    double a[4] = {0.97, 0.03, 0.05, 0.95};
    double b[6] = {0.80, 0.15, 0.05, 0.05, 0.15, 0.80};
    memcpy(h->pi, pi, sizeof(pi));
    memcpy(h->a, a, sizeof(a));
    memcpy(h->b, b, sizeof(b));
    hp_rng rng;
    hp_rng_seed(&rng, seed);
    uint32_t state = hp_rng_uniform(&rng) < pi[0] ? 0 : 1;
    for (size_t t = 0; t < T; ++t) {
        if (t > 0) {
            double u = hp_rng_uniform(&rng);
            state = (u < a[state * 2 + state]) ? state : 1 - state;
        }
        truth[t] = state;
        double u = hp_rng_uniform(&rng);
        double c = 0;
        uint32_t sym = 0;
        for (uint32_t k = 0; k < 3; ++k) {
            c += b[state * 3 + k];
            if (u <= c) {
                sym = k;
                break;
            }
            sym = k;
        }
        obs[t] = sym;
    }
    hp_hmm_free(h);
}

HP_TEST(test_discrete_train_recovery) {
    const size_t T = 3000;
    uint32_t *obs = hp_u32s(T), *truth = hp_u32s(T);
    sample_two_state(obs, truth, T, 20240601);

    hp_hmm *uniform = hp_hmm_new(2, 3);
    double *la = hp_doubles(T * 2);
    double ull = 0;
    HP_CHECK(hp_hmm_forward_log(uniform, obs, T, la, &ull) == HP_OK);

    hp_hmm *model = hp_hmm_new(2, 3);
    hp_hmm_train_opts opts;
    hp_hmm_train_opts_default(&opts);
    opts.max_iter = 200;
    opts.n_restarts = 4;
    opts.seed = 123;
    hp_hmm_report rep;
    HP_CHECK(hp_hmm_train(obs, T, 2, 3, &opts, model, &rep) == HP_OK);
    HP_CHECK(rep.loglik > ull);
    HP_CHECK(model->b[0] > 0.6);  /* state 0 dominates symbol 0  */
    HP_CHECK(model->b[5] > 0.6);  /* state 1 dominates symbol 2  */
    HP_CHECK(model->a[0] > 0.8);
    HP_CHECK(model->a[3] > 0.8);
    HP_CHECK(rep.converged == 1);

    double e0 = model->b[0] * 0 + model->b[1] * 1 + model->b[2] * 2;
    double e1 = model->b[3] * 0 + model->b[4] * 1 + model->b[5] * 2;
    HP_CHECK(e0 < e1);

    /* EM monotonicity: more iterations cannot lower the likelihood. */
    hp_hmm_train_opts o1 = opts;
    o1.max_iter = 1;
    hp_hmm *m1 = hp_hmm_new(2, 3);
    hp_hmm_report r1;
    HP_CHECK(hp_hmm_train(obs, T, 2, 3, &o1, m1, &r1) == HP_OK);
    HP_CHECK(rep.loglik >= r1.loglik - 1e-9);

    /* Viterbi decoding accuracy against the latent truth. */
    uint32_t *path = hp_u32s(T);
    HP_CHECK(hp_hmm_viterbi(model, obs, T, path, NULL) == HP_OK);
    size_t correct = 0;
    for (size_t t = 0; t < T; ++t) correct += (path[t] == truth[t]);
    double acc = (double)correct / (double)T;
    HP_CHECK(acc > 0.85);

    /* Save / load is exact. */
    const char *path_file = "build/test_tmp/hmm_model.txt";
    HP_CHECK(hp_hmm_save(model, path_file) == HP_OK);
    hp_hmm *loaded = hp_hmm_load(path_file);
    HP_CHECK(loaded != NULL);
    HP_CHECK(loaded->n_states == 2 && loaded->n_symbols == 3);
    for (size_t i = 0; i < 2; ++i) HP_CHECK_NEAR(loaded->pi[i], model->pi[i], 1e-15);
    for (size_t i = 0; i < 4; ++i) HP_CHECK_NEAR(loaded->a[i], model->a[i], 1e-15);
    for (size_t i = 0; i < 6; ++i) HP_CHECK_NEAR(loaded->b[i], model->b[i], 1e-15);
    double ll2 = 0;
    HP_CHECK(hp_hmm_forward_log(loaded, obs, T, la, &ll2) == HP_OK);
    HP_CHECK_NEAR(ll2, rep.loglik, 1e-9);

    hp_hmm_free(loaded);
    hp_hmm_free(m1);
    hp_free(path);
    hp_free(la);
    hp_free(obs);
    hp_free(truth);
    hp_hmm_free(model);
    hp_hmm_free(uniform);
    return 0;
}

HP_TEST(test_discrete_select_bic) {
    const size_t T = 2500;
    uint32_t *obs = hp_u32s(T), *truth = hp_u32s(T);
    sample_two_state(obs, truth, T, 777);
    hp_hmm_train_opts opts;
    hp_hmm_train_opts_default(&opts);
    opts.n_restarts = 3;
    opts.max_iter = 150;
    opts.seed = 9;
    hp_hmm *best = hp_hmm_new(1, 3);
    hp_hmm_model_info table[4];
    size_t best_n = 0;
    HP_CHECK(hp_hmm_select(obs, T, 3, 1, 4, &opts, best, &best_n, table) == HP_OK);
    HP_CHECK(best_n == 2);
    HP_CHECK(table[1].bic < table[0].bic);
    HP_CHECK(table[1].bic < table[2].bic);
    HP_CHECK(best->n_states == 2);
    hp_hmm_free(best);
    hp_free(obs);
    hp_free(truth);
    return 0;
}

void hp_test_hmm_register(void) {
    HP_REG(test_discrete_forward_vs_bruteforce);
    HP_REG(test_discrete_viterbi_vs_bruteforce);
    HP_REG(test_discrete_filter_vs_prefix);
    HP_REG(test_discrete_train_recovery);
    HP_REG(test_discrete_select_bic);
}
