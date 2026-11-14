/* Two-state, three-symbol HMM on a synthetic categorical sequence.
 *
 * Build: make examples
 * Run:   ./build/ex_discrete_hmm
 */
#include <math.h>
#include <stdio.h>

#include "hp/hp.h"

enum { T = 2000, N = 2, M = 3 };

int main(void) {
    double pi[N] = {0.6, 0.4};
    double A[N * N] = {0.97, 0.03, 0.05, 0.95};
    double B[N * M] = {0.80, 0.15, 0.05, 0.05, 0.15, 0.80};

    uint32_t *obs = hp_u32s(T);
    uint8_t *truth = hp_u8s(T);

    hp_rng rng;
    hp_rng_seed(&rng, 20240601);
    uint32_t state = hp_rng_uniform(&rng) < pi[0] ? 0 : 1;
    for (int t = 0; t < T; ++t) {
        if (t > 0) {
            double u = hp_rng_uniform(&rng);
            if (u >= A[state * N + state]) state = 1 - state;
        }
        truth[t] = (uint8_t)state;
        double u = hp_rng_uniform(&rng);
        double acc = 0.0;
        uint32_t sym = M - 1;
        for (uint32_t k = 0; k < M; ++k) {
            acc += B[state * M + k];
            if (u <= acc) {
                sym = k;
                break;
            }
        }
        obs[t] = sym;
    }

    hp_hmm *hmm = hp_hmm_new(N, M);
    hp_hmm_train_opts opt;
    hp_hmm_train_opts_default(&opt);
    opt.max_iter = 200;
    opt.n_restarts = 4;
    opt.seed = 123;

    hp_hmm_report rep;
    hp_status st = hp_hmm_train(obs, T, N, M, &opt, hmm, &rep);
    if (st != HP_OK) {
        fprintf(stderr, "train failed: %s\n", hp_status_name(st));
        hp_hmm_free(hmm);
        hp_free(obs);
        hp_free(truth);
        return 1;
    }

    uint32_t *path = hp_u32s(T);
    st = hp_hmm_viterbi(hmm, obs, T, path, NULL);
    if (st != HP_OK) {
        fprintf(stderr, "viterbi failed: %s\n", hp_status_name(st));
        hp_free(path);
        hp_hmm_free(hmm);
        hp_free(obs);
        hp_free(truth);
        return 1;
    }

    size_t correct = 0;
    for (int t = 0; t < T; ++t) correct += path[t] == truth[t];

    printf("logL = %.3f, BIC = %.3f, converged = %d\n", rep.loglik, rep.bic,
           rep.converged);
    printf("A = [[%.3f %.3f] [%.3f %.3f]]\n", hmm->a[0], hmm->a[1], hmm->a[2],
           hmm->a[3]);
    printf("B = [[%.3f %.3f %.3f] [%.3f %.3f %.3f]]\n", hmm->b[0], hmm->b[1],
           hmm->b[2], hmm->b[3], hmm->b[4], hmm->b[5]);
    printf("Viterbi accuracy = %.3f\n", (double)correct / (double)T);

    hp_free(path);
    hp_hmm_free(hmm);
    hp_free(obs);
    hp_free(truth);
    return 0;
}
