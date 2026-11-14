#include <math.h>

#include "hp/market/backtest.h"
#include "hp/market/strategy.h"
#include "hp_test.h"

HP_TEST(test_vol_target) {
    const size_t T = 10;
    double ret[10], pos[10];
    for (size_t t = 0; t < T; ++t) ret[t] = 0.001;

    hp_vol_target_opts o;
    hp_vol_target_opts_default(&o);
    o.lookback = 4;
    o.target_vol = 0.10;
    HP_CHECK(hp_signal_vol_target(ret, T, &o, pos) == HP_OK);
    HP_CHECK(isnan(pos[0]) && isnan(pos[2]));
    HP_CHECK_NEAR(pos[3], 1.0, 1e-15);

    for (size_t t = 0; t < T; ++t) ret[t] = (t % 2 == 0) ? 0.01 : -0.01;
    HP_CHECK(hp_signal_vol_target(ret, T, &o, pos) == HP_OK);
    HP_CHECK(pos[3] > 0.5 && pos[3] < 0.8);

    o.min_position = 0.2;
    o.max_position = 0.2;
    HP_CHECK(hp_signal_vol_target(ret, T, &o, pos) == HP_OK);
    HP_CHECK_NEAR(pos[4], 0.2, 1e-15);
    return 0;
}

HP_TEST(test_momentum) {
    const size_t T = 10;
    double ret[10], pos[10];
    for (size_t t = 0; t < 5; ++t) ret[t] = 0.01;
    for (size_t t = 5; t < T; ++t) ret[t] = -0.01;

    hp_momentum_opts o;
    hp_momentum_opts_default(&o);
    o.lookback = 5;
    HP_CHECK(hp_signal_momentum(ret, T, &o, pos) == HP_OK);
    HP_CHECK(isnan(pos[3]));
    HP_CHECK_NEAR(pos[4], 1.0, 1e-15);
    HP_CHECK_NEAR(pos[9], 0.0, 1e-15);

    double with_nan[10];
    with_nan[0] = NAN;
    for (size_t t = 1; t < T; ++t) with_nan[t] = 0.01;
    HP_CHECK(hp_signal_momentum(with_nan, T, &o, pos) == HP_OK);
    HP_CHECK_NEAR(pos[5], 1.0, 1e-15);
    return 0;
}

HP_TEST(test_state_exposure_and_multiply) {
    double filtered[6] = {1, 0, 0.5, 0.5, 0, 1};
    double exposure[2] = {0.1, 1.0};
    double pos[3];
    HP_CHECK(hp_signal_state_exposure(filtered, 3, 2, exposure, pos) == HP_OK);
    HP_CHECK_NEAR(pos[0], 0.1, 1e-15);
    HP_CHECK_NEAR(pos[1], 0.55, 1e-15);
    HP_CHECK_NEAR(pos[2], 1.0, 1e-15);

    double a[3] = {1.0, 0.5, NAN};
    double b[3] = {0.5, 0.5, 1.0};
    double out[3];
    HP_CHECK(hp_signal_multiply(a, b, 3, out) == HP_OK);
    HP_CHECK_NEAR(out[0], 0.5, 1e-15);
    HP_CHECK_NEAR(out[1], 0.25, 1e-15);
    HP_CHECK(isnan(out[2]));
    return 0;
}

HP_TEST(test_signal_no_lookahead) {
    const size_t T = 120;
    double a[120], b[120], pa[120], pb[120];
    for (size_t t = 0; t < T; ++t) {
        a[t] = 0.001 * sin((double)t);
        b[t] = a[t];
    }
    for (size_t t = 80; t < T; ++t) b[t] = -0.05;

    hp_vol_target_opts vo;
    hp_vol_target_opts_default(&vo);
    vo.lookback = 20;
    HP_CHECK(hp_signal_vol_target(a, T, &vo, pa) == HP_OK);
    HP_CHECK(hp_signal_vol_target(b, T, &vo, pb) == HP_OK);
    for (size_t t = 0; t < 80; ++t)
        HP_CHECK((isnan(pa[t]) && isnan(pb[t])) || pa[t] == pb[t]);

    hp_momentum_opts mo;
    hp_momentum_opts_default(&mo);
    mo.lookback = 20;
    HP_CHECK(hp_signal_momentum(a, T, &mo, pa) == HP_OK);
    HP_CHECK(hp_signal_momentum(b, T, &mo, pb) == HP_OK);
    for (size_t t = 0; t < 80; ++t)
        HP_CHECK((isnan(pa[t]) && isnan(pb[t])) || pa[t] == pb[t]);
    return 0;
}

HP_TEST(test_bt_opts_periods) {
    hp_bt_opts o;
    hp_bt_opts_default(&o);
    HP_CHECK_NEAR(o.periods_per_year, 252.0, 1e-15);
    o.periods_per_year = 12.0;
    o.apply_costs = 0;
    double ret[3] = {NAN, 0.01, 0.01};
    double sig[3] = {1, 1, 1};
    hp_bt_metrics m;
    hp_bt_result r;
    HP_CHECK(hp_backtest(ret, 3, sig, &o, &r, &m) == HP_OK);
    HP_CHECK_NEAR(m.final_equity, 1.0201, 1e-12);
    hp_bt_result_free(&r);
    return 0;
}

void hp_test_strategy_register(void) {
    HP_REG(test_vol_target);
    HP_REG(test_momentum);
    HP_REG(test_state_exposure_and_multiply);
    HP_REG(test_signal_no_lookahead);
    HP_REG(test_bt_opts_periods);
}
