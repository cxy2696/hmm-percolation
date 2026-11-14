#include <math.h>
#include <string.h>

#include "hp/core/memory.h"
#include "hp/core/rng.h"
#include "hp/core/stats.h"
#include "hp/hmm/ghmm.h"
#include "hp_test.h"

static void set_two_state(hp_ghmm *m) {
    m->n_states = 2;
    m->dim = 1;
    m->n_mix = 1;
    m->pi[0] = 0.5;
    m->pi[1] = 0.5;
    m->a[0] = 0.90;
    m->a[1] = 0.10;
    m->a[2] = 0.20;
    m->a[3] = 0.80;
    m->mix[0] = 1.0;
    m->mix[1] = 1.0;
    m->mean[0] = -1.0;
    m->mean[1] = 1.0;
    m->var[0] = 0.5;
    m->var[1] = 0.5;
}

static void brute_force(const hp_ghmm *m, const double *x, size_t T,
                        double *gamma, double *loglik, uint32_t *best_path,
                        double *best_logp) {
    const size_t n = m->n_states;
    size_t paths = 1;
    for (size_t t = 0; t < T; ++t) paths *= n;
    double total = 0.0, best = -INFINITY;
    size_t best_code = 0;
    if (gamma) memset(gamma, 0, T * n * sizeof(double));
    uint32_t *state = (uint32_t *)hp_u32s(T);
    for (size_t code = 0; code < paths; ++code) {
        size_t c = code;
        for (size_t t = 0; t < T; ++t) {
            state[t] = (uint32_t)(c % n);
            c /= n;
        }
        double p = m->pi[state[0]] *
                   hp_norm_pdf(x[0], m->mean[state[0]], m->var[state[0]]);
        for (size_t t = 1; t < T; ++t) {
            p *= m->a[state[t - 1] * n + state[t]];
            p *= hp_norm_pdf(x[t], m->mean[state[t]], m->var[state[t]]);
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

HP_TEST(test_ghmm_forward_bruteforce) {
    hp_ghmm *m = hp_ghmm_new(2, 1, 1);
    set_two_state(m);
    HP_CHECK(hp_ghmm_validate(m) == HP_OK);
    double x[5] = {-1.2, 0.8, 1.1, -0.9, 0.7};
    double la[10], ll = 0;
    HP_CHECK(hp_ghmm_forward_log(m, x, 5, la, &ll) == HP_OK);
    double brute = 0;
    brute_force(m, x, 5, NULL, &brute, NULL, NULL);
    HP_CHECK_NEAR(ll, brute, 1e-12);

    double gamma[10], gbrute[10];
    double pll = 0;
    HP_CHECK(hp_ghmm_posterior(m, x, 5, gamma, &pll) == HP_OK);
    HP_CHECK_NEAR(pll, brute, 1e-12);
    brute_force(m, x, 5, gbrute, NULL, NULL, NULL);
    double total = exp(brute);
    for (size_t t = 0; t < 5; ++t) {
        for (size_t i = 0; i < 2; ++i)
            HP_CHECK_NEAR(gamma[t * 2 + i], gbrute[t * 2 + i] / total, 1e-10);
    }
    hp_ghmm_free(m);
    return 0;
}

HP_TEST(test_ghmm_viterbi_bruteforce) {
    hp_ghmm *m = hp_ghmm_new(2, 1, 1);
    set_two_state(m);
    double x[5] = {-1.2, 0.8, 1.1, -0.9, 0.7};
    uint32_t path[5], bpath[5];
    double lp = 0, blp = 0;
    HP_CHECK(hp_ghmm_viterbi(m, x, 5, path, &lp) == HP_OK);
    brute_force(m, x, 5, NULL, NULL, bpath, &blp);
    HP_CHECK_NEAR(lp, blp, 1e-12);
    for (size_t t = 0; t < 5; ++t) HP_CHECK(path[t] == bpath[t]);
    hp_ghmm_free(m);
    return 0;
}

HP_TEST(test_ghmm_filter_and_step) {
    hp_ghmm *m = hp_ghmm_new(2, 1, 1);
    set_two_state(m);
    double x[5] = {-1.2, 0.8, 1.1, -0.9, 0.7};
    double filt[10];
    HP_CHECK(hp_ghmm_filter(m, x, 5, filt, NULL) == HP_OK);
    for (size_t t = 0; t < 5; ++t) {
        double prefix_gamma[10];
        brute_force(m, x, t + 1, prefix_gamma, NULL, NULL, NULL);
        double s = 0;
        for (size_t i = 0; i < 2; ++i) s += prefix_gamma[t * 2 + i];
        for (size_t i = 0; i < 2; ++i)
            HP_CHECK_NEAR(filt[t * 2 + i], prefix_gamma[t * 2 + i] / s, 1e-10);
    }
    /* t=0 has no transition; seed from pi * emission, then step. */
    double step[2];
    double b0[2];
    HP_CHECK(hp_ghmm_log_emission(m, x, b0) == HP_OK);
    double z = 0.0;
    for (size_t i = 0; i < 2; ++i) {
        step[i] = m->pi[i] * exp(b0[i]);
        z += step[i];
    }
    for (size_t i = 0; i < 2; ++i) step[i] /= z;
    HP_CHECK_NEAR(step[0], filt[0], 1e-10);
    HP_CHECK_NEAR(step[1], filt[1], 1e-10);
    for (size_t t = 1; t < 5; ++t) {
        HP_CHECK(hp_ghmm_filter_step(m, x + t, step, NULL) == HP_OK);
        for (size_t i = 0; i < 2; ++i)
            HP_CHECK_NEAR(step[i], filt[t * 2 + i], 1e-10);
    }
    hp_ghmm_free(m);
    return 0;
}

static void sample_regimes(double *x, uint32_t *truth, size_t T, uint64_t seed) {
    hp_rng rng;
    hp_rng_seed(&rng, seed);
    uint32_t state = hp_rng_uniform(&rng) < 0.5 ? 0 : 1;
    for (size_t t = 0; t < T; ++t) {
        if (t > 0) {
            double stay = state == 0 ? 0.97 : 0.93;
            if (hp_rng_uniform(&rng) >= stay) state = 1 - state;
        }
        truth[t] = state;
        double mu = state == 0 ? -2.0 : 2.0;
        x[t] = mu + 0.7 * hp_rng_normal(&rng);
    }
}

HP_TEST(test_ghmm_train_recovery) {
    const size_t T = 2000;
    double *x = hp_doubles(T);
    uint32_t *truth = hp_u32s(T);
    sample_regimes(x, truth, T, 20240701);

    hp_ghmm *model = hp_ghmm_new(1, 1, 1);
    hp_ghmm_train_opts opts;
    hp_ghmm_train_opts_default(&opts);
    opts.max_iter = 150;
    opts.n_restarts = 3;
    opts.seed = 11;
    hp_ghmm_report rep;
    HP_CHECK(hp_ghmm_train(x, T, 1, 2, 1, &opts, model, &rep) == HP_OK);
    HP_CHECK(model->mean[0] < -1.5);
    HP_CHECK(model->mean[1] > 1.5);
    HP_CHECK(model->a[0] > 0.85);
    HP_CHECK(model->a[3] > 0.85);
    HP_CHECK(rep.converged == 1);
    HP_CHECK(rep.loglik > -4000.0);

    uint32_t *path = hp_u32s(T);
    HP_CHECK(hp_ghmm_viterbi(model, x, T, path, NULL) == HP_OK);
    size_t correct = 0;
    for (size_t t = 0; t < T; ++t) correct += (path[t] == truth[t]);
    HP_CHECK((double)correct / (double)T > 0.90);

    const char *path_file = "build/test_tmp/ghmm_model.txt";
    HP_CHECK(hp_ghmm_save(model, path_file) == HP_OK);
    hp_ghmm *loaded = hp_ghmm_load(path_file);
    HP_CHECK(loaded != NULL);
    HP_CHECK(loaded->n_states == 2 && loaded->dim == 1 && loaded->n_mix == 1);
    for (size_t i = 0; i < 2; ++i) HP_CHECK_NEAR(loaded->mean[i], model->mean[i], 1e-15);
    for (size_t i = 0; i < 4; ++i) HP_CHECK_NEAR(loaded->a[i], model->a[i], 1e-15);
    double ll1 = 0, ll2 = 0;
    double *la = hp_doubles(T * 2);
    HP_CHECK(hp_ghmm_forward_log(model, x, T, la, &ll1) == HP_OK);
    HP_CHECK(hp_ghmm_forward_log(loaded, x, T, la, &ll2) == HP_OK);
    HP_CHECK_NEAR(ll1, ll2, 1e-9);
    HP_CHECK_NEAR(ll1, rep.loglik, 1e-6);

    hp_ghmm_free(loaded);
    hp_free(path);
    hp_free(la);
    hp_free(truth);
    hp_free(x);
    hp_ghmm_free(model);
    return 0;
}

HP_TEST(test_ghmm_select_bic) {
    const size_t T = 2000;
    double *x = hp_doubles(T);
    uint32_t *truth = hp_u32s(T);
    sample_regimes(x, truth, T, 99);
    hp_ghmm_train_opts opts;
    hp_ghmm_train_opts_default(&opts);
    opts.n_restarts = 2;
    opts.max_iter = 100;
    opts.seed = 5;
    hp_ghmm *best = hp_ghmm_new(1, 1, 1);
    hp_ghmm_model_info table[3];
    size_t best_n = 0;
    HP_CHECK(hp_ghmm_select(x, T, 1, 1, 3, 1, &opts, best, &best_n, table) == HP_OK);
    HP_CHECK(best_n == 2);
    HP_CHECK(table[1].bic < table[0].bic);
    HP_CHECK(table[1].bic < table[2].bic);
    hp_ghmm_free(best);
    hp_free(x);
    hp_free(truth);
    return 0;
}

HP_TEST(test_ghmm_mixture_recovery) {
    const size_t T = 4000;
    double *x = hp_doubles(T);
    hp_rng rng;
    hp_rng_seed(&rng, 1234);
    for (size_t t = 0; t < T; ++t) {
        double mu = hp_rng_uniform(&rng) < 0.5 ? -3.0 : 3.0;
        x[t] = mu + 0.5 * hp_rng_normal(&rng);
    }
    hp_ghmm *m = hp_ghmm_new(1, 1, 2);
    hp_ghmm_train_opts opts;
    hp_ghmm_train_opts_default(&opts);
    opts.n_restarts = 4;
    opts.max_iter = 200;
    opts.seed = 21;
    hp_ghmm_report rep;
    HP_CHECK(hp_ghmm_train(x, T, 1, 1, 2, &opts, m, &rep) == HP_OK);
    double m0 = m->mean[0], m1 = m->mean[1];
    HP_CHECK(m0 < -2.0 || m1 < -2.0);
    HP_CHECK(m0 > 2.0 || m1 > 2.0);
    HP_CHECK_NEAR(m->mix[0] + m->mix[1], 1.0, 1e-12);
    HP_CHECK(m->mix[0] > 0.3 && m->mix[1] > 0.3);
    hp_ghmm_free(m);
    hp_free(x);
    return 0;
}

void hp_test_ghmm_register(void) {
    HP_REG(test_ghmm_forward_bruteforce);
    HP_REG(test_ghmm_viterbi_bruteforce);
    HP_REG(test_ghmm_filter_and_step);
    HP_REG(test_ghmm_train_recovery);
    HP_REG(test_ghmm_select_bic);
    HP_REG(test_ghmm_mixture_recovery);
}
