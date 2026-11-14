/* HP · deterministic synthetic market with latent regimes.
 *
 * The generator exists to make every number in this repository reproducible
 * and to give the HMM a ground truth to recover.  It is not a market model:
 * three sticky regimes set the drift/volatility, assets are one-factor with
 * regime-dependent idiosyncratic risk (correlations rise in the bear state),
 * and innovations are Student-t to keep fat tails.
 */
#ifndef HP_MARKET_SYNTH_H
#define HP_MARKET_SYNTH_H

#include <stdint.h>

#include "hp/market/panel.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t n_days;
    size_t n_assets;
    size_t n_regimes; /* 2 or 3 */
    uint64_t seed;
} hp_synth_opts;

void hp_synth_opts_default(hp_synth_opts *opts);

typedef struct {
    hp_panel *panel;
    uint8_t *regime;    /* T */
    double *market_ret; /* T */
    size_t n_regimes;
    double *factor_beta; /* N */
} hp_synth_data;

hp_synth_data *hp_synth_generate(const hp_synth_opts *opts);
void hp_synth_free(hp_synth_data *d);
hp_status hp_synth_write_truth(const hp_synth_data *d, const char *path);

#ifdef __cplusplus
}
#endif

#endif /* HP_MARKET_SYNTH_H */
