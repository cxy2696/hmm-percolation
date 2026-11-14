/* One-dimensional two-regime Gaussian HMM.
 *
 * This example deliberately uses raw (non-standardised) data to show the
 * variance floor; the application pipeline standardises features instead.
 *
 * Build: make examples
 * Run:   ./build/ex_gaussian_hmm
 */
#include <math.h>
#include <stdio.h>

#include "hp/hp.h"

enum { T = 1200, N = 2, D = 1, K = 1 };

int main(void) {
    double *x = hp_doubles(T);
    uint8_t *truth = hp_u8s(T);

    hp_rng rng;
    hp_rng_seed(&rng, 7);
    uint8_t state = hp_rng_uniform(&rng) < 0.5 ? 0 : 1;
    for (int t = 0; t < T; ++t) {
        if (t > 0) {
            double stay = state == 0 ? 0.97 : 0.94;
            if (hp_rng_uniform(&rng) >= stay) state = 1 - state;
        }
        truth[t] = state;
        x[t] = (state == 0 ? -1.5 : 1.5) + 0.6 * hp_rng_normal(&rng);
    }

    hp_ghmm *m = hp_ghmm_new(N, D, K);
    hp_ghmm_train_opts opt;
    hp_ghmm_train_opts_default(&opt);
    opt.max_iter = 150;
    opt.n_restarts = 4;
    opt.seed = 11;
    opt.var_floor = 0.05;

    hp_ghmm_report rep;
    hp_status st = hp_ghmm_train(x, T, D, N, K, &opt, m, &rep);
    if (st != HP_OK) {
        fprintf(stderr, "train failed: %s\n", hp_status_name(st));
        hp_ghmm_free(m);
        hp_free(x);
        hp_free(truth);
        return 1;
    }

    uint32_t *path = hp_u32s(T);
    hp_ghmm_viterbi(m, x, T, path, NULL);
    size_t same = 0, flipped = 0;
    for (int t = 0; t < T; ++t) {
        same += path[t] == truth[t];
        flipped += path[t] == (uint32_t)(1 - truth[t]);
    }
    double acc = (double)(same > flipped ? same : flipped) / (double)T;

    printf("logL = %.3f, converged = %d, iterations = %u\n", rep.loglik,
           rep.converged, rep.iterations);
    printf("means = [%.3f %.3f], vars = [%.3f %.3f]\n", m->mean[0], m->mean[1],
           m->var[0], m->var[1]);
    printf("A = [[%.3f %.3f] [%.3f %.3f]]\n", m->a[0], m->a[1], m->a[2], m->a[3]);
    printf("best-label accuracy = %.3f\n", acc);

    hp_free(path);
    hp_ghmm_free(m);
    hp_free(x);
    hp_free(truth);
    return 0;
}
