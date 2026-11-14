#include <math.h>
#include <string.h>

#include "hp/core/csv.h"
#include "hp/core/io.h"
#include "hp/core/memory.h"
#include "hp/core/rng.h"
#include "hp/core/stats.h"
#include "hp/cross/state_graph.h"
#include "hp/market/backtest.h"
#include "hp/market/features.h"
#include "hp/market/panel.h"
#include "hp/market/synth.h"
#include "hp_test.h"

HP_TEST(test_synth_determinism) {
    hp_synth_opts o;
    hp_synth_opts_default(&o);
    o.n_days = 300;
    o.n_assets = 4;
    o.n_regimes = 3;
    o.seed = 123;
    hp_synth_data *a = hp_synth_generate(&o);
    hp_synth_data *b = hp_synth_generate(&o);
    HP_CHECK(a && b);
    HP_CHECK(a->panel->T == 300 && a->panel->N == 4);
    HP_CHECK_STREQ(a->panel->symbols[0], "A1");
    HP_CHECK_STREQ(a->panel->dates[0], "2012-01-02");
    for (size_t i = 0; i < 300 * 4; ++i)
        HP_CHECK(a->panel->close[i] == b->panel->close[i]);
    for (size_t t = 0; t < 300; ++t) {
        HP_CHECK(a->regime[t] == b->regime[t]);
        HP_CHECK(a->regime[t] < 3);
        HP_CHECK(a->panel->close[t * 4] > 0.0);
    }
    o.seed = 124;
    hp_synth_data *c = hp_synth_generate(&o);
    int differs = 0;
    for (size_t i = 0; i < 300 * 4; ++i)
        if (a->panel->close[i] != c->panel->close[i]) differs = 1;
    HP_CHECK(differs);
    hp_synth_free(a);
    hp_synth_free(b);
    hp_synth_free(c);
    return 0;
}

HP_TEST(test_panel_csv_roundtrip) {
    hp_synth_opts o;
    hp_synth_opts_default(&o);
    o.n_days = 150;
    o.n_assets = 3;
    o.seed = 77;
    hp_synth_data *d = hp_synth_generate(&o);
    HP_CHECK(d != NULL);
    const char *path = "build/test_tmp/panel.csv";
    HP_CHECK(hp_panel_save_csv(d->panel, path) == HP_OK);
    hp_status st = HP_ERR_IO;
    hp_panel *p = hp_panel_load_csv(path, &st);
    HP_CHECK(st == HP_OK && p != NULL);
    HP_CHECK(p->T == d->panel->T && p->N == d->panel->N);
    HP_CHECK_STREQ(p->symbols[2], d->panel->symbols[2]);
    HP_CHECK_STREQ(p->dates[149], d->panel->dates[149]);
    for (size_t i = 0; i < p->T * p->N; ++i)
        HP_CHECK_NEAR(p->close[i], d->panel->close[i], 1e-9);
    for (size_t i = p->N; i < p->T * p->N; ++i)
        HP_CHECK_NEAR(p->ret[i], d->panel->ret[i], 1e-9);
    double ew[150];
    HP_CHECK(hp_panel_ew_returns(p, ew) == HP_OK);
    HP_CHECK(isfinite(ew[10]));
    HP_CHECK(hp_panel_symbol_index(p, "A2") == 1);
    HP_CHECK(hp_panel_symbol_index(p, "zz") == -1);
    hp_panel_free(p);

    const char *truth = "build/test_tmp/truth.csv";
    HP_CHECK(hp_synth_write_truth(d, truth) == HP_OK);
    hp_csv csv;
    HP_CHECK(hp_csv_read(truth, &csv) == HP_OK);
    HP_CHECK(csv.n_rows == 150);
    HP_CHECK(csv.n_cols == 6); /* date, regime, market_ret, 3 assets */
    HP_CHECK_STREQ(csv.header[1], "regime");
    hp_csv_free(&csv);
    hp_free(d);
    return 0;
}

HP_TEST(test_features) {
    const size_t T = 300, w = 20;
    double *ret = hp_doubles(T);
    hp_rng rng;
    hp_rng_seed(&rng, 5);
    for (size_t t = 0; t < T; ++t) ret[t] = 0.001 + 0.01 * hp_rng_normal(&rng);
    hp_features *f = hp_features_from_returns(ret, T, w);
    HP_CHECK(f != NULL);
    HP_CHECK(f->dim == 2);
    HP_CHECK(f->start == w);
    HP_CHECK(isnan(f->x[(w - 1) * 2]));
    HP_CHECK(isfinite(f->x[w * 2]));
    HP_CHECK_NEAR(f->x[w * 2], ret[w] * 100.0, 1e-12);
    HP_CHECK(f->vol[w] > 0.0);
    HP_CHECK_NEAR(f->x[w * 2 + 1], log(f->vol[w] * 100.0), 1e-12);
    HP_CHECK(isnan(f->x[0]));
    hp_features_free(f);
    hp_free(ret);
    return 0;
}

HP_TEST(test_backtest_timing) {
    const size_t T = 8;
    double ret[8] = {NAN, 0, 0, 0, 0, 0.01, 0, 0};
    double sig[8] = {0, 0, 0, 0, 1, 0, 0, 0};
    hp_bt_opts o;
    hp_bt_opts_default(&o);
    o.apply_costs = 0;
    hp_bt_result r;
    hp_bt_metrics m;
    HP_CHECK(hp_backtest(ret, T, sig, &o, &r, &m) == HP_OK);
    HP_CHECK_NEAR(r.equity[4], 1.0, 1e-15);
    HP_CHECK_NEAR(r.equity[5], 1.01, 1e-12);
    HP_CHECK_NEAR(r.equity[7], 1.01, 1e-12);
    HP_CHECK_NEAR(m.total_return, 0.01, 1e-12);
    HP_CHECK(m.n_days == 7);
    hp_bt_result_free(&r);

    /* Same return, but the signal only appears after the move: no capture. */
    double sig2[8] = {0, 0, 0, 0, 0, 1, 0, 0};
    HP_CHECK(hp_backtest(ret, T, sig2, &o, &r, &m) == HP_OK);
    HP_CHECK_NEAR(r.equity[7], 1.0, 1e-15);
    hp_bt_result_free(&r);

    double ret2[4] = {NAN, 0.01, 0.01, 0.01};
    double hold[4] = {1, 1, 1, 1};
    o.apply_costs = 1;
    o.cost_bps = 10.0;
    HP_CHECK(hp_backtest(ret2, 4, hold, &o, &r, &m) == HP_OK);
    HP_CHECK_NEAR(r.equity[1], 1.0 * (1.0 + 0.01 - 0.001), 1e-12);
    HP_CHECK(m.avg_turnover > 0.0);
    HP_CHECK_NEAR(m.exposure, 1.0, 1e-15);
    HP_CHECK_NEAR(m.hit_rate, 1.0, 1e-15);
    hp_bt_result_free(&r);
    return 0;
}

HP_TEST(test_transition_graph) {
    double a[16] = {0.9, 0.09, 0.005, 0.005,
                    0.09, 0.9, 0.005, 0.005,
                    0.005, 0.005, 0.8, 0.19,
                    0.005, 0.005, 0.19, 0.8};
    hp_components c;
    HP_CHECK(hp_transition_components(a, 4, 0.05, &c) == HP_OK);
    HP_CHECK(c.n_components == 2);
    HP_CHECK(c.giant == 2);
    hp_components_free(&c);
    HP_CHECK(hp_transition_components(a, 4, 0.25, &c) == HP_OK);
    HP_CHECK(c.n_components == 4);
    hp_components_free(&c);
    return 0;
}

void hp_test_market_register(void) {
    HP_REG(test_synth_determinism);
    HP_REG(test_panel_csv_roundtrip);
    HP_REG(test_features);
    HP_REG(test_backtest_timing);
    HP_REG(test_transition_graph);
}
