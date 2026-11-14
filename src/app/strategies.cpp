#include "hp/app/strategies.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "hp/app/util.hpp"
#include "hp/market/strategy.h"

namespace hp::app {
namespace {

void check(hp_status st, const char *what) {
    if (st != HP_OK)
        throw std::runtime_error(std::string(what) + ": " + hp_status_message(st));
}

void fill_warmup(std::vector<double> &v) {
    for (double &x : v)
        if (!std::isfinite(x)) x = 0.0;
}

void clip_unit(std::vector<double> &v) {
    for (double &x : v) x = std::clamp(x, 0.0, 1.0);
}

}  // namespace

StrategyComparison compare_strategies(const RegimeResult &r,
                                      const StrategyConfig &cfg) {
    if (r.bt_ret.empty() || r.bt_dates.size() != r.bt_ret.size())
        throw std::runtime_error("regime result has no OOS slice");
    const size_t s0 = r.train_end;
    const size_t T = r.ret.size();
    const size_t len = r.bt_ret.size();
    if (s0 + len != T) throw std::runtime_error("OOS slice mismatch");

    const double *ret = r.ret.data();
    std::vector<double> vol_target(T), momentum(T);
    hp_vol_target_opts vo;
    hp_vol_target_opts_default(&vo);
    vo.lookback = cfg.vol_lookback;
    vo.target_vol = cfg.target_vol;
    vo.min_position = cfg.min_position;
    vo.max_position = cfg.max_position;
    check(hp_signal_vol_target(ret, T, &vo, vol_target.data()), "vol target");
    fill_warmup(vol_target);

    hp_momentum_opts mo;
    hp_momentum_opts_default(&mo);
    mo.lookback = cfg.momentum_lookback;
    mo.threshold = cfg.momentum_threshold;
    check(hp_signal_momentum(ret, T, &mo, momentum.data()), "momentum");
    fill_warmup(momentum);

    std::vector<double> bh(len, 1.0), constant(len, cfg.constant_exposure);
    std::vector<double> gate_only(len);
    for (size_t i = 0; i < len; ++i) gate_only[i] = r.gate[s0 + i];
    std::vector<double> hmm(r.position_hmm.begin() + s0, r.position_hmm.end());
    std::vector<double> comb(r.position_comb.begin() + s0, r.position_comb.end());
    std::vector<double> vt(vol_target.begin() + s0, vol_target.end());
    std::vector<double> mom(momentum.begin() + s0, momentum.end());
    std::vector<double> hmm_vol(len), hmm_vol_gate(len);
    for (size_t i = 0; i < len; ++i) {
        hmm_vol[i] = hmm[i] * vt[i];
        hmm_vol_gate[i] = hmm_vol[i] * r.gate[s0 + i];
    }
    clip_unit(bh);
    clip_unit(constant);
    clip_unit(gate_only);
    clip_unit(hmm);
    clip_unit(comb);
    clip_unit(vt);
    clip_unit(mom);
    clip_unit(hmm_vol);
    clip_unit(hmm_vol_gate);

    struct Spec {
        const char *name;
        const std::vector<double> *sig;
    };
    const Spec specs[] = {
        {"buy_hold", &bh},
        {"constant", &constant},
        {"vol_target", &vt},
        {"momentum", &mom},
        {"gate_only", &gate_only},
        {"hmm", &hmm},
        {"hmm_gate", &comb},
        {"hmm_vol", &hmm_vol},
        {"hmm_vol_gate", &hmm_vol_gate},
    };

    hp_bt_opts bo;
    hp_bt_opts_default(&bo);
    bo.cost_bps = cfg.cost_bps;
    bo.apply_costs = 1;

    StrategyComparison out;
    out.oos_start = s0;
    out.dates = r.bt_dates;
    out.ret = r.bt_ret;
    for (const Spec &s : specs) {
        hp_bt_result br{};
        StrategyResult row;
        row.name = s.name;
        check(hp_backtest(r.bt_ret.data(), len, s.sig->data(), &bo,
                          &br, &row.metrics),
              "hp_backtest");
        row.equity.assign(br.equity, br.equity + len);
        row.position.assign(br.position, br.position + len);
        hp_bt_result_free(&br);
        out.rows.push_back(std::move(row));
    }
    return out;
}

}  // namespace hp::app
