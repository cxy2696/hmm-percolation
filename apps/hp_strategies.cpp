/* hp_strategies · compare a few causal position rules on one OOS window.
 *
 * The point is not to find "the best" rule but to keep the comparison honest:
 * same dates, same costs, same execution lag, no per-strategy date fiddling.
 */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hp/app/args.hpp"
#include "hp/app/regime.hpp"
#include "hp/app/report.hpp"
#include "hp/app/strategies.hpp"
#include "hp/app/util.hpp"
#include "hp/core/csv.h"

using namespace hp::app;

namespace {

std::string pct(double v) {
    if (!std::isfinite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f%%", 100.0 * v);
    return std::string(buf);
}

std::string num(double v, int p = 3) {
    if (!std::isfinite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", p, v);
    return std::string(buf);
}

std::string row(const std::string &name, const hp_bt_metrics &m) {
    std::ostringstream os;
    os << "| " << name << " | " << pct(m.total_return) << " | " << pct(m.cagr)
       << " | " << pct(m.ann_vol) << " | " << num(m.sharpe, 2) << " | "
       << pct(m.max_drawdown) << " | " << num(m.avg_turnover, 3) << " | "
       << pct(m.exposure) << " |\n";
    return os.str();
}

void json_metrics(Json &j, const hp_bt_metrics &m) {
    j.obj();
    j.kv("total_return", m.total_return);
    j.kv("cagr", m.cagr);
    j.kv("ann_vol", m.ann_vol);
    j.kv("sharpe", m.sharpe);
    j.kv("max_drawdown", m.max_drawdown);
    j.kv("calmar", m.calmar);
    j.kv("avg_turnover", m.avg_turnover);
    j.kv("exposure", m.exposure);
    j.kv("hit_rate", m.hit_rate);
    j.kv("final_equity", m.final_equity);
    j.kv("n_days", m.n_days);
    j.end();
}

}  // namespace

int main(int argc, char **argv) {
    try {
        Args args(argc, argv);
        std::string csv = args.require("csv");
        RegimeConfig rc;
        rc.asset = args.get("asset", rc.asset);
        rc.states = (size_t)args.get_int("states", (long)rc.states);
        rc.n_mix = (size_t)args.get_int("mix", (long)rc.n_mix);
        rc.vol_window = (size_t)args.get_int("vol-window", (long)rc.vol_window);
        rc.train_frac = args.get_double("train-frac", rc.train_frac);
        rc.refit_every = (size_t)args.get_int("refit", (long)rc.refit_every);
        rc.max_iter = (unsigned)args.get_int("max-iter", (long)rc.max_iter);
        rc.restarts = (unsigned)args.get_int("restarts", (long)rc.restarts);
        rc.seed = (uint64_t)args.get_int("seed", (long)rc.seed);
        rc.cost_bps = args.get_double("cost-bps", rc.cost_bps);
        rc.use_percolation = !args.has("no-perc");
        rc.conn.window = (size_t)args.get_int("conn-window", (long)rc.conn.window);
        rc.conn.theta_min = args.get_double("conn-theta-min", rc.conn.theta_min);
        rc.conn.theta_max = args.get_double("conn-theta-max", rc.conn.theta_max);
        rc.conn.theta_ref = args.get_double("conn-theta-ref", rc.conn.theta_ref);
        rc.conn.theta_steps = (size_t)args.get_int("conn-theta-steps", (long)rc.conn.theta_steps);
        rc.risk_quantile = args.get_double("risk-quantile", rc.risk_quantile);
        rc.risk_gate = args.get_double("risk-gate", rc.risk_gate);

        StrategyConfig sc;
        sc.constant_exposure = args.get_double("constant", sc.constant_exposure);
        sc.vol_lookback = (size_t)args.get_int("vol-lookback", (long)sc.vol_lookback);
        sc.target_vol = args.get_double("vol-target", sc.target_vol);
        sc.min_position = args.get_double("min-position", sc.min_position);
        sc.max_position = args.get_double("max-position", sc.max_position);
        sc.momentum_lookback = (size_t)args.get_int("momentum-lookback", (long)sc.momentum_lookback);
        sc.momentum_threshold = args.get_double("momentum-threshold", sc.momentum_threshold);
        sc.cost_bps = rc.cost_bps;

        std::string outdir = args.get("outdir", "out/strategies");
        ensure_dir(outdir);
        PanelPtr panel = load_panel(csv);
        RegimeResult r = run_regime(panel.get(), rc);
        StrategyComparison c = compare_strategies(r, sc);

        /* metrics.csv */
        std::string mpath = path_join(outdir, "strategies.csv");
        const char *mhead[] = {"strategy", "total_return", "cagr", "ann_vol",
                               "sharpe", "max_drawdown", "calmar", "turnover",
                               "exposure", "hit_rate", "n_days"};
        hp_csv_writer *mw = hp_csv_writer_open(mpath.c_str(), mhead, 11, 1);
        if (!mw) throw std::runtime_error("cannot write " + mpath);
        for (const StrategyResult &s : c.rows) {
            double vals[10] = {s.metrics.total_return, s.metrics.cagr,
                               s.metrics.ann_vol, s.metrics.sharpe,
                               s.metrics.max_drawdown, s.metrics.calmar,
                               s.metrics.avg_turnover, s.metrics.exposure,
                               s.metrics.hit_rate, (double)s.metrics.n_days};
            if (hp_csv_writer_row_label(mw, s.name.c_str(), vals, 10) != HP_OK)
                throw std::runtime_error("cannot write strategy row");
        }
        hp_csv_writer_close(mw);

        /* equity.csv and positions.csv */
        std::vector<const char *> eh = {"date"};
        std::vector<const char *> ph = {"date"};
        for (const StrategyResult &s : c.rows) {
            eh.push_back(s.name.c_str());
            ph.push_back(s.name.c_str());
        }
        hp_csv_writer *ew = hp_csv_writer_open(path_join(outdir, "equity.csv").c_str(),
                                               eh.data(), eh.size(), 1);
        hp_csv_writer *pw = hp_csv_writer_open(path_join(outdir, "positions.csv").c_str(),
                                               ph.data(), ph.size(), 1);
        if (!ew || !pw) throw std::runtime_error("cannot write equity/positions");
        std::vector<double> ev(c.rows.size()), pv(c.rows.size());
        for (size_t i = 0; i < c.dates.size(); ++i) {
            for (size_t j = 0; j < c.rows.size(); ++j) {
                ev[j] = c.rows[j].equity[i];
                pv[j] = c.rows[j].position[i];
            }
            hp_csv_writer_row_label(ew, c.dates[i].c_str(), ev.data(), ev.size());
            hp_csv_writer_row_label(pw, c.dates[i].c_str(), pv.data(), pv.size());
        }
        hp_csv_writer_close(ew);
        hp_csv_writer_close(pw);

        /* summary.json */
        Json j;
        j.obj();
        j.kv("strategy_count", c.rows.size());
        j.kv("oos_start", c.oos_start);
        j.kv("oos_days", c.ret.size());
        j.key("config"); j.obj();
        j.kv("cost_bps", sc.cost_bps);
        j.kv("constant", sc.constant_exposure);
        j.kv("vol_lookback", sc.vol_lookback);
        j.kv("target_vol", sc.target_vol);
        j.kv("min_position", sc.min_position);
        j.kv("max_position", sc.max_position);
        j.kv("momentum_lookback", sc.momentum_lookback);
        j.kv("momentum_threshold", sc.momentum_threshold);
        j.kv("hmm_states", r.n_states);
        j.kv("use_percolation", rc.use_percolation);
        j.end();
        j.key("strategies"); j.arr();
        for (const StrategyResult &s : c.rows) {
            j.obj();
            j.kv("name", s.name);
            j.key("metrics"); json_metrics(j, s.metrics); j.end();
        }
        j.end();
        j.end();
        write_text(path_join(outdir, "summary.json"), j.dump());

        /* report.md */
        std::ostringstream md;
        md << "# Strategy comparison\n\n";
        const size_t pnl_days = c.rows.empty() ? 0 : c.rows[0].metrics.n_days;
        md << "First decision: " << (c.dates.empty() ? "n/a" : c.dates.front())
           << ". Signals: " << c.ret.size() << " rows; P&L days: " << pnl_days
           << ". Cost: "
           << num(sc.cost_bps, 1) << " bps per unit turnover. "
           << "All rows use signal[t] on return[t+1].\n\n";
        md << "| strategy | total | CAGR | vol | Sharpe | max DD | turnover | exposure |\n";
        md << "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n";
        for (const StrategyResult &s : c.rows) md << row(s.name, s.metrics);
        md << "\nThe HMM rows are the only ones that use a fitted state model; "
              "`buy_hold`, `constant`, `vol_target` and `momentum` are there to "
              "make the comparison falsifiable. A higher Sharpe on one synthetic "
              "sample is not evidence of a durable edge.\n";
        write_text(path_join(outdir, "report.md"), md.str());

        const size_t pnl_n = c.rows.empty() ? 0 : c.rows[0].metrics.n_days;
        std::printf("first decision %s, signals %zu, P&L days %zu\n",
                    c.dates.empty() ? "n/a" : c.dates.front().c_str(),
                    c.ret.size(), pnl_n);
        std::printf("%-12s %8s %8s %8s %8s %8s\n", "strategy", "total", "CAGR",
                    "vol", "Sharpe", "maxDD");
        for (const StrategyResult &s : c.rows) {
            std::printf("%-12s %7.2f%% %7.2f%% %7.2f%% %8.2f %7.2f%%\n",
                        s.name.c_str(), 100 * s.metrics.total_return,
                        100 * s.metrics.cagr, 100 * s.metrics.ann_vol,
                        s.metrics.sharpe, 100 * s.metrics.max_drawdown);
        }
        std::printf("wrote %s\n", outdir.c_str());
    } catch (const std::exception &e) {
        std::fprintf(stderr, "hp_strategies: %s\n", e.what());
        return 1;
    }
    return 0;
}
