#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hp/app/args.hpp"
#include "hp/app/connectedness.hpp"
#include "hp/app/regime.hpp"
#include "hp/app/report.hpp"
#include "hp/app/strategies.hpp"
#include "hp/app/util.hpp"
#include "hp/core/io.h"
#include "hp_test.h"

using namespace hp::app;

namespace {

PanelPtr clone_panel(const hp_panel *src) {
    hp_panel *p = hp_panel_new(src->T, src->N);
    for (size_t t = 0; t < src->T; ++t) p->dates[t] = hp_strdup(src->dates[t]);
    for (size_t j = 0; j < src->N; ++j) p->symbols[j] = hp_strdup(src->symbols[j]);
    std::memcpy(p->close, src->close, src->T * src->N * sizeof(double));
    std::memcpy(p->ret, src->ret, src->T * src->N * sizeof(double));
    return PanelPtr(p, hp_panel_free);
}

}  // namespace

HP_TEST(test_args) {
    char a0[] = "prog", a1[] = "--csv", a2[] = "x.csv", a3[] = "--flag",
         a4[] = "--n", a5[] = "3", a6[] = "--rate=0.5";
    char *argv[] = {a0, a1, a2, a3, a4, a5, a6};
    Args args(7, argv);
    HP_CHECK(args.get("csv") == "x.csv");
    HP_CHECK(args.has("flag"));
    HP_CHECK(args.get_int("n", 0) == 3);
    HP_CHECK_NEAR(args.get_double("rate", 0.0), 0.5, 1e-15);
    HP_CHECK(!args.has("missing"));
    return 0;
}

HP_TEST(test_json) {
    Json j;
    j.obj();
    j.kv("a", 1);
    j.key("b");
    j.arr();
    j.val(1.5);
    j.val("x");
    j.end();
    j.kv("bad", NAN);
    j.end();
    HP_CHECK(j.dump() == "{\"a\":1,\"b\":[1.5,\"x\"],\"bad\":null}");
    return 0;
}

HP_TEST(test_connectedness_deterministic) {
    const std::string path = "data/synthetic_panel.csv";
    if (!std::filesystem::exists(path)) return 0;
    PanelPtr p = load_panel(path);
    ConnectednessConfig cfg;
    cfg.window = 60;
    cfg.step = 5;
    cfg.theta_steps = 31;
    Connectedness conn(cfg);
    ConnectednessResult a = conn.run(p.get());
    ConnectednessResult b = conn.run(p.get());
    HP_CHECK(!a.points.empty());
    HP_CHECK(a.points.size() == b.points.size());
    for (size_t k = 0; k < a.points.size(); ++k) {
        HP_CHECK(a.points[k].theta_c == b.points[k].theta_c);
        HP_CHECK(a.points[k].giant_ref == b.points[k].giant_ref);
        HP_CHECK(a.points[k].theta_c > 0.5 && a.points[k].theta_c < 1.0);
        HP_CHECK(a.points[k].giant_ref >= 0.0 && a.points[k].giant_ref <= 1.0);
    }
    return 0;
}

HP_TEST(test_regime_no_lookahead) {
    const std::string path = "data/synthetic_panel.csv";
    if (!std::filesystem::exists(path)) return 0;
    hp_status st = HP_OK;
    hp_panel *raw = hp_panel_load_csv(path.c_str(), &st);
    HP_CHECK(raw != nullptr && st == HP_OK);
    PanelPtr p1(raw, hp_panel_free);
    PanelPtr p2 = clone_panel(p1.get());
    const size_t T = p2->T, N = p2->N;
    const size_t cut = T - 120;
    for (size_t t = cut; t < T; ++t)
        for (size_t j = 0; j < N; ++j) p2->ret[t * N + j] *= -4.0;

    RegimeConfig cfg;
    cfg.states = 2;
    cfg.max_iter = 60;
    cfg.restarts = 2;
    cfg.refit_every = 250;
    cfg.use_percolation = false;
    cfg.train_frac = 0.6;
    RegimeResult a = run_regime(p1.get(), cfg);
    RegimeResult b = run_regime(p2.get(), cfg);
    HP_CHECK(a.n_states == 2 && b.n_states == 2);
    HP_CHECK(a.position_hmm.size() == T);
    HP_CHECK(a.bt_hmm.n_days > 0);
    for (size_t t = 0; t < cut; ++t)
        HP_CHECK_NEAR(a.position_hmm[t], b.position_hmm[t], 1e-12);
    return 0;
}

HP_TEST(test_regime_deterministic) {
    const std::string path = "data/synthetic_panel.csv";
    if (!std::filesystem::exists(path)) return 0;
    PanelPtr p = load_panel(path);
    RegimeConfig cfg;
    cfg.states = 2;
    cfg.max_iter = 40;
    cfg.restarts = 2;
    cfg.refit_every = 500;
    cfg.use_percolation = false;
    RegimeResult a = run_regime(p.get(), cfg);
    RegimeResult b = run_regime(p.get(), cfg);
    HP_CHECK_NEAR(a.report.loglik, b.report.loglik, 1e-12);
    for (size_t t = 0; t < p->T; ++t)
        HP_CHECK(a.viterbi[t] == b.viterbi[t]);
    HP_CHECK_NEAR(a.bt_hmm.total_return, b.bt_hmm.total_return, 1e-12);
    return 0;
}

HP_TEST(test_strategy_comparison) {
    const std::string path = "data/synthetic_panel.csv";
    if (!std::filesystem::exists(path)) return 0;
    PanelPtr p = load_panel(path);
    RegimeConfig cfg;
    cfg.states = 2;
    cfg.max_iter = 60;
    cfg.restarts = 2;
    cfg.refit_every = 250;
    cfg.use_percolation = false;
    RegimeResult r = run_regime(p.get(), cfg);

    StrategyConfig sc;
    sc.cost_bps = cfg.cost_bps;
    StrategyComparison c = compare_strategies(r, sc);
    HP_CHECK(c.rows.size() == 9);
    HP_CHECK(c.dates.size() == r.bt_dates.size());
    HP_CHECK(c.rows[0].name == "buy_hold");
    for (const StrategyResult &row : c.rows) {
        HP_CHECK(row.metrics.n_days == r.bt_hmm.n_days);
        HP_CHECK(row.equity.size() == r.bt_ret.size());
        for (double x : row.position)
            HP_CHECK(std::isfinite(x) && x >= 0.0 && x <= 1.0);
    }
    HP_CHECK_NEAR(c.rows[0].metrics.total_return, r.bt_bh.total_return, 1e-12);
    HP_CHECK_NEAR(c.rows[6].metrics.total_return, r.bt_comb.total_return, 1e-12);

    StrategyComparison d = compare_strategies(r, sc);
    for (size_t k = 0; k < c.rows.size(); ++k)
        HP_CHECK_NEAR(c.rows[k].metrics.total_return, d.rows[k].metrics.total_return, 1e-12);
    return 0;
}

int main() {
    struct Entry {
        const char *name;
        int (*fn)();
    } tests[] = {
        {"test_args", test_args},
        {"test_json", test_json},
        {"test_connectedness_deterministic", test_connectedness_deterministic},
        {"test_regime_no_lookahead", test_regime_no_lookahead},
        {"test_regime_deterministic", test_regime_deterministic},
        {"test_strategy_comparison", test_strategy_comparison},
    };
    int failed = 0;
    for (const auto &t : tests) {
        std::printf("[ RUN  ] %s\n", t.name);
        int rc = t.fn();
        std::printf(rc == 0 ? "[  OK  ] %s\n" : "[ FAIL ] %s\n", t.name);
        failed += rc != 0;
    }
    std::printf("\n%d passed, %d failed\n", (int)(sizeof(tests) / sizeof(tests[0])) - failed,
                failed);
    return failed == 0 ? 0 : 1;
}
