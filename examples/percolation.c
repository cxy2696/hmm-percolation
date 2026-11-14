/* A few percolation calculations, independent of finance.
 *
 * Build: make examples
 * Run:   ./build/ex_percolation
 */
#include <math.h>
#include <stdio.h>

#include "hp/hp.h"

int main(void) {
    const int L = 64;
    const int reps = 50;

    hp_lattice *lat = hp_lattice_new(L, HP_PERC_BOND, 0);
    if (!lat) return 1;
    hp_rng rng;
    hp_rng_seed_stream(&rng, 1, 0);
    double s = 0.0, span = 0.0, chi = 0.0;
    for (int r = 0; r < reps; ++r) {
        hp_lattice_bernoulli(lat, 0.55, &rng);
        hp_lattice_stats st;
        hp_lattice_analyze(lat, &st);
        s += st.s_max;
        span += st.spanning ? 1.0 : 0.0;
        chi += st.susceptibility;
    }
    printf("2D bond p=0.55, L=%d, reps=%d: S=%.3f span=%.2f chi=%.3f\n", L, reps,
           s / reps, span / reps, chi / reps);
    hp_lattice_free(lat);

    hp_perc_exp e = {.L = L, .kind = HP_PERC_BOND, .periodic = 0,
                     .reps = reps, .seed = 2};
    double pc = 0.0, sp = 0.0;
    if (hp_perc_crossing(&e, 0.35, 0.65, 12, &pc, &sp) == HP_OK)
        printf("bond L=%d crossing = %.4f (span at cross %.2f)\n", L, pc, sp);

    printf("ER c=2 giant fraction = %.6f\n", hp_perc_er_giant_theory(2.0));

    enum { KMAX = 30 };
    double pk[KMAX + 1];
    double fact = 1.0, term = exp(-2.0), sum = 0.0;
    pk[0] = term;
    for (int k = 1; k <= KMAX; ++k) {
        fact *= (double)k;
        term = exp(-2.0) * pow(2.0, (double)k) / fact;
        pk[k] = term;
        sum += term;
    }
    double total = pk[0] + sum;
    for (int k = 0; k <= KMAX; ++k) pk[k] /= total;

    double giant = 0.0, r0 = 0.0, mean = 0.0;
    if (hp_perc_config_giant(pk, KMAX, 0.6, &giant, &r0, &mean) == HP_OK)
        printf("config model c=2, T=0.6: R0=%.3f S=%.3f\n", r0, giant);
    return 0;
}
