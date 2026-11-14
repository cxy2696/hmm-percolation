#include "hp/app/regime.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "hp/core/memory.h"
#include "hp/core/stats.h"
#include "hp/cross/state_graph.h"
#include "hp/market/features.h"

namespace hp::app {
namespace {

void check(hp_status st, const char *what) {
    if (st != HP_OK)
        throw std::runtime_error(std::string(what) + " failed: " + hp_status_name(st));
}

struct FeaturesDeleter {
    void operator()(hp_features *f) const { hp_features_free(f); }
};
using FeaturesPtr = std::unique_ptr<hp_features, FeaturesDeleter>;

GhmmPtr fit_model(const double *z, size_t n_rows, size_t n_states, size_t n_mix,
                  const hp_ghmm_train_opts &opts, hp_ghmm_report *rep) {
    GhmmPtr m = make_ghmm(n_states, 2, n_mix);
    hp_status st = hp_ghmm_train(z, n_rows, 2, n_states, n_mix, &opts, m.get(), rep);
    check(st, "hp_ghmm_train");
    return m;
}

}  // namespace

std::string regime_label(size_t state, size_t n_states) {
    if (n_states == 2) return state == 0 ? "弱势" : "强势";
    if (n_states == 3) {
        static const char *labels[3] = {"弱势", "震荡", "强势"};
        return labels[state < 3 ? state : 2];
    }
    return "状态" + std::to_string(state + 1);
}

RegimeResult run_regime(const hp_panel *panel, const RegimeConfig &cfg) {
    if (!panel) throw std::runtime_error("null panel");
    if (cfg.states != 0 && cfg.states < 2) throw std::runtime_error("states must be 0 or >= 2");
    if (cfg.n_mix < 1) throw std::runtime_error("n_mix must be >= 1");
    const size_t T = panel->T;
    std::vector<double> ret(T, std::numeric_limits<double>::quiet_NaN());
    if (cfg.asset == "EW") {
        check(hp_panel_ew_returns(panel, ret.data()), "hp_panel_ew_returns");
    } else {
        int idx = hp_panel_symbol_index(panel, cfg.asset.c_str());
        if (idx < 0) throw std::runtime_error("unknown asset: " + cfg.asset);
        check(hp_panel_return_series(panel, (size_t)idx, ret.data()), "hp_panel_return_series");
    }

    FeaturesPtr feats(hp_features_from_returns(ret.data(), T, cfg.vol_window));
    if (!feats) throw std::runtime_error("cannot build features");
    const size_t start = feats->start;
    if (start >= T || T - start < 300)
        throw std::runtime_error("not enough valid feature rows (need >= 300)");
    for (size_t t = start; t < T; ++t) {
        if (!std::isfinite(feats->x[t * 2]) || !std::isfinite(feats->x[t * 2 + 1]))
            throw std::runtime_error("feature row " + std::to_string(t) + " is not finite");
    }
    const size_t n_days = T - start;

    size_t train_len = (size_t)((double)n_days * cfg.train_frac);
    if (train_len < 150) train_len = 150;
    if (train_len > n_days - 60) train_len = n_days > 60 ? n_days - 60 : n_days;
    const size_t train_end = start + train_len;

    std::vector<double> mean(2, 0.0), sd(2, 1.0);
    for (size_t d = 0; d < 2; ++d) {
        double s = 0.0, s2 = 0.0;
        for (size_t t = start; t < train_end; ++t) {
            double v = feats->x[t * 2 + d];
            s += v;
            s2 += v * v;
        }
        mean[d] = s / (double)train_len;
        double var = (s2 - (double)train_len * mean[d] * mean[d]) / (double)(train_len - 1);
        sd[d] = std::sqrt(std::max(var, 1e-12));
    }
    std::vector<double> z(n_days * 2);
    for (size_t t = start; t < T; ++t)
        for (size_t d = 0; d < 2; ++d)
            z[(t - start) * 2 + d] = (feats->x[t * 2 + d] - mean[d]) / sd[d];

    hp_ghmm_train_opts opts;
    hp_ghmm_train_opts_default(&opts);
    opts.max_iter = cfg.max_iter;
    opts.n_restarts = cfg.restarts;
    opts.seed = cfg.seed;
    opts.init = HP_GHMM_INIT_QUANTILE;
    opts.verbose = 0;

    RegimeResult out;
    out.start = start;
    out.train_end = train_end;
    out.oos_start = train_end + 1; /* first day with an out-of-sample return */
    out.feature_mean = mean;
    out.feature_sd = sd;
    out.ret = ret;
    out.vol.assign(T, std::numeric_limits<double>::quiet_NaN());
    for (size_t t = 0; t < T; ++t) out.vol[t] = feats->vol[t];

    size_t n = cfg.states;
    if (n == 0) {
        size_t nmin = 2, nmax = 5;
        if (nmax > train_len / 40) nmax = std::max<size_t>(nmin, train_len / 40);
        std::vector<hp_ghmm_model_info> table(nmax - nmin + 1);
        size_t best_n = 0;
        GhmmPtr best = make_ghmm(nmin, 2, cfg.n_mix);
        hp_status st = hp_ghmm_select(z.data(), train_len, 2, nmin, nmax, cfg.n_mix,
                                      &opts, best.get(), &best_n, table.data());
        check(st, "hp_ghmm_select");
        n = best_n;
        out.model_train = std::move(best);
        out.bic_table = std::move(table);
        const hp_ghmm_model_info &info = out.bic_table[n - nmin];
        out.report.loglik = info.loglik;
        out.report.bic = info.bic;
        out.report.aic = info.aic;
        out.report.n_params = info.n_params;
        out.report.iterations = 0;
        out.report.converged = info.converged;
        out.report.restarts = cfg.restarts;
    } else {
        hp_ghmm_report rep{};
        out.model_train = fit_model(z.data(), train_len, n, cfg.n_mix, opts, &rep);
        out.report = rep;
        hp_ghmm_model_info info{};
        info.n_states = n;
        info.loglik = rep.loglik;
        info.bic = rep.bic;
        info.aic = rep.aic;
        info.n_params = rep.n_params;
        info.converged = rep.converged;
        out.bic_table.push_back(info);
    }
    if (n < 2 || n > 8) throw std::runtime_error("unreasonable selected state count");
    out.n_states = n;
    /* Structural diagnostic: threshold A at its median off-diagonal
     * transition, i.e. keep the stronger half of the inter-state links. */
    std::vector<double> offdiag;
    offdiag.reserve(n * (n - 1));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (i != j) offdiag.push_back(out.model_train->a[i * n + j]);
    std::sort(offdiag.begin(), offdiag.end());
    if (!offdiag.empty()) {
        size_t mid = offdiag.size() / 2;
        out.state_graph_theta = offdiag.size() % 2
                                    ? offdiag[mid]
                                    : 0.5 * (offdiag[mid - 1] + offdiag[mid]);
    }
    hp_components sgc;
    if (hp_transition_components(out.model_train->a, n, out.state_graph_theta, &sgc) == HP_OK) {
        out.state_n_components = sgc.n_components;
        out.state_giant = sgc.giant;
        out.state_giant_fraction = sgc.giant_fraction;
        hp_components_free(&sgc);
    }
    out.model_final = GhmmPtr(hp_ghmm_clone(out.model_train.get()), hp_ghmm_free);

    out.posterior.assign(T * n, std::numeric_limits<double>::quiet_NaN());
    out.filtered.assign(T * n, std::numeric_limits<double>::quiet_NaN());
    out.viterbi.assign(T, UINT32_MAX);
    std::vector<double> post(n_days * n), filt(n_days * n);
    std::vector<uint32_t> vpath(n_days);
    double ll = 0.0;
    check(hp_ghmm_posterior(out.model_train.get(), z.data(), n_days, post.data(), &ll),
          "hp_ghmm_posterior");
    check(hp_ghmm_filter(out.model_train.get(), z.data(), n_days, filt.data(), nullptr),
          "hp_ghmm_filter");
    check(hp_ghmm_viterbi(out.model_train.get(), z.data(), n_days, vpath.data(), nullptr),
          "hp_ghmm_viterbi");
    for (size_t t = start; t < T; ++t) {
        for (size_t j = 0; j < n; ++j) {
            out.posterior[t * n + j] = post[(t - start) * n + j];
            out.filtered[t * n + j] = filt[(t - start) * n + j];
        }
        out.viterbi[t] = vpath[t - start];
    }

    std::vector<double> expo(n);
    for (size_t j = 0; j < n; ++j)
        expo[j] = n == 1 ? 1.0 : 0.1 + 0.9 * (double)j / (double)(n - 1);

    out.theta_c.assign(T, std::numeric_limits<double>::quiet_NaN());
    out.giant_ref.assign(T, std::numeric_limits<double>::quiet_NaN());
    out.gate.assign(T, 1.0);
    out.position_hmm.assign(T, 0.0);
    out.position_comb.assign(T, 0.0);

    if (cfg.use_percolation) {
        Connectedness conn(cfg.conn);
        ConnectednessResult cr = conn.run(panel);
        out.conn_points = cr.points.size();
        for (size_t t = 0; t < T; ++t) {
            out.theta_c[t] = cr.theta_c[t];
            out.giant_ref[t] = cr.giant_ref[t];
        }
        std::vector<double> hist;
        hist.reserve(cfg.risk_lookback);
        for (size_t t = 0; t < T; ++t) {
            if (std::isfinite(out.theta_c[t])) {
                hist.clear();
                for (size_t s = t; s-- > 0 && hist.size() < cfg.risk_lookback;)
                    if (std::isfinite(out.theta_c[s])) hist.push_back(out.theta_c[s]);
                if (hist.size() >= cfg.risk_min_obs) {
                    std::sort(hist.begin(), hist.end());
                    double q = hp_quantile_sorted(hist.data(), hist.size(), cfg.risk_quantile);
                    if (out.theta_c[t] > q) out.gate[t] = cfg.risk_gate;
                }
            }
        }
    }

    for (size_t t = start; t < T; ++t) {
        double p = 0.0;
        for (size_t j = 0; j < n; ++j) p += out.filtered[t * n + j] * expo[j];
        out.position_hmm[t] = p;
        out.position_comb[t] = p * out.gate[t];
    }

    /* Walk-forward refits and causal filtered probabilities. */
    size_t block = train_end;
    std::vector<double> fbuf(n_days * n);
    while (block < T) {
        size_t end = std::min(T, block + std::max<size_t>(1, cfg.refit_every));
        hp_ghmm_report rep2{};
        GhmmPtr m = fit_model(z.data(), block - start, n, cfg.n_mix, opts, &rep2);
        for (size_t t = block; t < end; ++t) {
            size_t len = t - start + 1;
            check(hp_ghmm_filter(m.get(), z.data(), len, fbuf.data(), nullptr),
                  "hp_ghmm_filter(walk-forward)");
            for (size_t j = 0; j < n; ++j) out.filtered[t * n + j] = fbuf[(len - 1) * n + j];
            double p = 0.0;
            for (size_t j = 0; j < n; ++j) p += out.filtered[t * n + j] * expo[j];
            out.position_hmm[t] = p;
            out.position_comb[t] = p * out.gate[t];
        }
        out.model_final = std::move(m);
        block = end;
    }

    /* Descriptive regime statistics (full sample, smoothed). */
    out.stats.resize(n);
    std::vector<double> sumw(n, 0.0), sumr(n, 0.0), sumr2(n, 0.0), sumv(n, 0.0);
    for (size_t t = start; t < T; ++t) {
        for (size_t j = 0; j < n; ++j) {
            double g = out.posterior[t * n + j];
            sumw[j] += g;
            sumr[j] += g * ret[t];
            sumr2[j] += g * ret[t] * ret[t];
            sumv[j] += g * feats->vol[t];
        }
    }
    std::vector<double> run_sum(n, 0.0);
    std::vector<size_t> run_cnt(n, 0);
    size_t run_len = 0;
    uint32_t run_state = UINT32_MAX;
    for (size_t t = start; t < T; ++t) {
        uint32_t s = out.viterbi[t];
        if (s == run_state) {
            ++run_len;
        } else {
            if (run_state != UINT32_MAX && run_state < n) {
                run_sum[run_state] += (double)run_len;
                ++run_cnt[run_state];
            }
            run_state = s;
            run_len = 1;
        }
    }
    if (run_state != UINT32_MAX && run_state < n) {
        run_sum[run_state] += (double)run_len;
        ++run_cnt[run_state];
    }
    for (size_t j = 0; j < n; ++j) {
        RegimeStat &s = out.stats[j];
        s.state = j;
        s.self_transition = out.model_train->a[j * n + j];
        s.expected_duration = s.self_transition < 1.0 ? 1.0 / (1.0 - s.self_transition) : std::numeric_limits<double>::infinity();
        s.observed_duration = run_cnt[j] ? run_sum[j] / (double)run_cnt[j] : 0.0;
        if (sumw[j] > 1e-12) {
            s.mean_daily = sumr[j] / sumw[j];
            double var = sumr2[j] / sumw[j] - s.mean_daily * s.mean_daily;
            s.ann_return = s.mean_daily * 252.0;
            s.ann_vol = std::sqrt(std::max(var, 0.0) * 252.0);
            s.sharpe = s.ann_vol > 0.0 ? s.ann_return / s.ann_vol : 0.0;
            s.frac_days = sumw[j] / (double)n_days;
            s.mean_vol = sumv[j] / sumw[j];
        }
    }

    /* OOS backtest: the first out-of-sample decision is made at the close of
     * train_end, so the first P&L day is train_end+1.  The slice therefore
     * starts one day early and carries one initial flat equity point. */
    if (train_end + 1 >= T) throw std::runtime_error("not enough out-of-sample days");
    size_t s0 = train_end;
    size_t len = T - s0;
    std::vector<double> br(len), bh(len, 1.0), sh(len), sc(len);
    for (size_t i = 0; i < len; ++i) {
        br[i] = ret[s0 + i];
        sh[i] = out.position_hmm[s0 + i];
        sc[i] = out.position_comb[s0 + i];
        out.bt_dates.push_back(panel->dates[s0 + i]);
    }
    out.bt_ret = br;
    hp_bt_opts bo;
    hp_bt_opts_default(&bo);
    bo.cost_bps = cfg.cost_bps;
    bo.apply_costs = 1;

    auto run_bt = [&](const std::vector<double> &sig, hp_bt_metrics *metrics,
                      std::vector<double> *eq) {
        hp_bt_result r{};
        check(hp_backtest(br.data(), len, sig.data(), &bo, &r, metrics), "hp_backtest");
        eq->assign(r.equity, r.equity + len);
        hp_bt_result_free(&r);
    };
    run_bt(bh, &out.bt_bh, &out.eq_bh);
    run_bt(sh, &out.bt_hmm, &out.eq_hmm);
    run_bt(sc, &out.bt_comb, &out.eq_comb);

    return out;
}

}  // namespace hp::app
