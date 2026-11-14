#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "hp/app/args.hpp"
#include "hp/app/report.hpp"
#include "hp/app/regime.hpp"
#include "hp/app/util.hpp"
#include "hp/core/csv.h"
#include "hp/core/stats.h"

using namespace hp::app;

namespace {

const char *nfn(double v, int prec = 4) {
    static char bufs[8][64];
    static int slot = 0;
    slot = (slot + 1) % 8;
    if (!std::isfinite(v)) {
        std::snprintf(bufs[slot], sizeof(bufs[slot]), "n/a");
    } else {
        std::snprintf(bufs[slot], sizeof(bufs[slot]), "%.*f", prec, v);
    }
    return bufs[slot];
}

std::string pct(double v) {
    if (!std::isfinite(v)) return "n/a";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f%%", 100.0 * v);
    return std::string(buf);
}

void json_doubles(Json &j, const std::vector<double> &v) {
    j.arr();
    for (double x : v) j.val(x);
    j.end();
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

std::string metrics_row(const std::string &name, const hp_bt_metrics &m) {
    std::ostringstream os;
    os << "| " << name << " | " << pct(m.total_return) << " | " << pct(m.cagr) << " | "
       << pct(m.ann_vol) << " | " << nfn(m.sharpe, 2) << " | " << pct(m.max_drawdown)
       << " | " << nfn(m.calmar, 2) << " | " << nfn(m.avg_turnover, 3) << " | "
       << pct(m.exposure) << " | " << pct(m.hit_rate) << " |\n";
    return os.str();
}

}  // namespace

int main(int argc, char **argv) {
    try {
        Args args(argc, argv);
        std::string csv = args.require("csv");
        RegimeConfig cfg;
        cfg.asset = args.get("asset", cfg.asset);
        cfg.states = (size_t)args.get_int("states", (long)cfg.states);
        cfg.n_mix = (size_t)args.get_int("mix", (long)cfg.n_mix);
        cfg.vol_window = (size_t)args.get_int("vol-window", (long)cfg.vol_window);
        cfg.train_frac = args.get_double("train-frac", cfg.train_frac);
        cfg.refit_every = (size_t)args.get_int("refit", (long)cfg.refit_every);
        cfg.max_iter = (unsigned)args.get_int("max-iter", (long)cfg.max_iter);
        cfg.restarts = (unsigned)args.get_int("restarts", (long)cfg.restarts);
        cfg.seed = (uint64_t)args.get_int("seed", (long)cfg.seed);
        cfg.cost_bps = args.get_double("cost-bps", cfg.cost_bps);
        cfg.use_percolation = !args.has("no-perc");
        cfg.conn.window = (size_t)args.get_int("conn-window", (long)cfg.conn.window);
        cfg.conn.theta_min = args.get_double("conn-theta-min", cfg.conn.theta_min);
        cfg.conn.theta_max = args.get_double("conn-theta-max", cfg.conn.theta_max);
        cfg.conn.theta_ref = args.get_double("conn-theta-ref", cfg.conn.theta_ref);
        cfg.conn.theta_steps = (size_t)args.get_int("conn-theta-steps", (long)cfg.conn.theta_steps);
        cfg.risk_quantile = args.get_double("risk-quantile", cfg.risk_quantile);
        cfg.risk_gate = args.get_double("risk-gate", cfg.risk_gate);
        std::string outdir = args.get("outdir", "out/regime");
        ensure_dir(outdir);

        PanelPtr panel = load_panel(csv);
        RegimeResult r = run_regime(panel.get(), cfg);
        const size_t T = panel->T, n = r.n_states, dim = 2;

        /* states.csv */
        std::vector<std::string> names = {"date", "ret", "vol"};
        for (size_t j = 0; j < n; ++j) names.push_back("p" + std::to_string(j));
        names.push_back("viterbi");
        names.push_back("position_hmm");
        names.push_back("gate");
        names.push_back("position_comb");
        names.push_back("theta_c");
        names.push_back("giant_ref");
        std::vector<const char *> cnames;
        for (const auto &s : names) cnames.push_back(s.c_str());
        std::string spath = path_join(outdir, "states.csv");
        hp_csv_writer *w = hp_csv_writer_open(spath.c_str(), cnames.data(), cnames.size(), 1);
        if (!w) throw std::runtime_error("cannot write " + spath);
        std::vector<double> vals(cnames.size() - 1, NAN);
        for (size_t t = 0; t < T; ++t) {
            std::fill(vals.begin(), vals.end(), NAN);
            vals[0] = r.ret[t];
            if (t >= r.start) {
                vals[1] = r.vol[t];
                for (size_t j = 0; j < n; ++j) vals[2 + j] = r.posterior[t * n + j];
                if (r.viterbi[t] != std::numeric_limits<uint32_t>::max())
                    vals[2 + n] = (double)r.viterbi[t];
                vals[3 + n] = r.position_hmm[t];
                vals[4 + n] = r.gate[t];
                vals[5 + n] = r.position_comb[t];
                vals[6 + n] = r.theta_c[t];
                vals[7 + n] = r.giant_ref[t];
            }
            hp_status st = hp_csv_writer_row_label(w, panel->dates[t], vals.data(), vals.size());
            if (st != HP_OK) {
                hp_csv_writer_abort(w);
                throw std::runtime_error("write states.csv");
            }
        }
        if (hp_csv_writer_close(w) != HP_OK) throw std::runtime_error("close states.csv");

        /* equity.csv */
        std::string epath = path_join(outdir, "equity.csv");
        const char *enames[8] = {"date", "ret", "buy_hold", "hmm", "hmm_perc",
                                 "pos_hmm", "pos_comb", "gate"};
        hp_csv_writer *ew = hp_csv_writer_open(epath.c_str(), enames, 8, 1);
        if (!ew) throw std::runtime_error("cannot write " + epath);
        for (size_t i = 0; i < r.eq_bh.size(); ++i) {
            size_t t = r.oos_start - 1 + i;
            double ev[7] = {r.bt_ret[i], r.eq_bh[i], r.eq_hmm[i], r.eq_comb[i],
                            r.position_hmm[t], r.position_comb[t], r.gate[t]};
            hp_status st = hp_csv_writer_row_label(ew, r.bt_dates[i].c_str(), ev, 7);
            if (st != HP_OK) {
                hp_csv_writer_abort(ew);
                throw std::runtime_error("write equity.csv");
            }
        }
        if (hp_csv_writer_close(ew) != HP_OK) throw std::runtime_error("close equity.csv");

        /* bic.csv */
        std::string bpath = path_join(outdir, "bic.csv");
        const char *bnames[6] = {"n_states", "params", "loglik", "bic", "aic", "converged"};
        hp_csv_writer *bw = hp_csv_writer_open(bpath.c_str(), bnames, 6, 0);
        if (!bw) throw std::runtime_error("cannot write " + bpath);
        for (const auto &info : r.bic_table) {
            double bv[6] = {(double)info.n_states, (double)info.n_params, info.loglik,
                            info.bic, info.aic, (double)info.converged};
            hp_csv_writer_row(bw, bv, 6);
        }
        hp_csv_writer_close(bw);

        /* models */
        std::string m1 = path_join(outdir, "model_train.txt");
        std::string m2 = path_join(outdir, "model_final.txt");
        if (hp_ghmm_save(r.model_train.get(), m1.c_str()) != HP_OK)
            throw std::runtime_error("save model_train.txt");
        if (hp_ghmm_save(r.model_final.get(), m2.c_str()) != HP_OK)
            throw std::runtime_error("save model_final.txt");

        /* percolation summary */
        std::vector<double> theta_finite;
        for (size_t t = r.oos_start; t < T; ++t)
            if (std::isfinite(r.theta_c[t])) theta_finite.push_back(r.theta_c[t]);
        std::sort(theta_finite.begin(), theta_finite.end());
        double theta_med = theta_finite.empty() ? NAN
                          : hp_quantile_sorted(theta_finite.data(), theta_finite.size(), 0.5);
        double theta_mean = theta_finite.empty() ? NAN : hp_mean(theta_finite.data(), theta_finite.size());
        size_t gate_days = 0;
        for (size_t t = r.oos_start; t < T; ++t)
            if (r.gate[t] < 1.0) ++gate_days;
        double gate_frac = (double)gate_days / (double)(T - r.oos_start);

        /* metrics.json */
        Json j;
        j.obj();
        j.kv("tool", "hp_regime");
        j.kv("version", "1.0.0");
        j.kv("generated_at", now_utc());
        j.kv("csv", csv);
        j.key("data");
        j.obj();
        j.kv("asset", cfg.asset);
        j.kv("T", T);
        j.kv("N", panel->N);
        j.kv("start", r.start);
        j.kv("train_end", r.train_end);
        j.kv("oos_start", r.oos_start);
        j.kv("start_date", std::string(panel->dates[r.start]));
        j.kv("train_end_date", std::string(panel->dates[r.train_end]));
        j.kv("end_date", std::string(panel->dates[T - 1]));
        j.end();
        j.key("config");
        j.obj();
        j.kv("states", n);
        j.kv("mix", cfg.n_mix);
        j.kv("vol_window", cfg.vol_window);
        j.kv("train_frac", cfg.train_frac);
        j.kv("refit_every", cfg.refit_every);
        j.kv("seed", cfg.seed);
        j.kv("cost_bps", cfg.cost_bps);
        j.kv("use_percolation", cfg.use_percolation);
        j.kv("risk_quantile", cfg.risk_quantile);
        j.kv("risk_gate", cfg.risk_gate);
        j.end();
        j.key("features");
        j.obj();
        j.key("mean");
        json_doubles(j, r.feature_mean);
        j.key("sd");
        json_doubles(j, r.feature_sd);
        j.end();
        j.key("model");
        j.obj();
        j.kv("n_states", n);
        j.kv("n_mix", cfg.n_mix);
        j.kv("dim", dim);
        j.kv("loglik", r.report.loglik);
        j.kv("bic", r.report.bic);
        j.kv("aic", r.report.aic);
        j.kv("n_params", r.report.n_params);
        j.kv("converged", r.report.converged != 0);
        j.kv("iterations", r.report.iterations);
        j.key("pi");
        j.arr();
        for (size_t i = 0; i < n; ++i) j.val(r.model_train->pi[i]);
        j.end();
        j.key("A");
        j.arr();
        for (size_t i = 0; i < n; ++i) {
            j.arr();
            for (size_t k = 0; k < n; ++k) j.val(r.model_train->a[i * n + k]);
            j.end();
        }
        j.end();
        j.key("mix");
        j.arr();
        for (size_t i = 0; i < n; ++i) {
            j.arr();
            for (size_t c = 0; c < cfg.n_mix; ++c) j.val(r.model_train->mix[i * cfg.n_mix + c]);
            j.end();
        }
        j.end();
        j.key("mean");
        j.arr();
        for (size_t i = 0; i < n; ++i) {
            j.arr();
            for (size_t c = 0; c < cfg.n_mix; ++c) {
                j.arr();
                for (size_t d = 0; d < dim; ++d)
                    j.val(r.model_train->mean[(i * cfg.n_mix + c) * dim + d]);
                j.end();
            }
            j.end();
        }
        j.end();
        j.key("var");
        j.arr();
        for (size_t i = 0; i < n; ++i) {
            j.arr();
            for (size_t c = 0; c < cfg.n_mix; ++c) {
                j.arr();
                for (size_t d = 0; d < dim; ++d)
                    j.val(r.model_train->var[(i * cfg.n_mix + c) * dim + d]);
                j.end();
            }
            j.end();
        }
        j.end();
        j.end();
        j.key("state_graph");
        j.obj();
        j.kv("theta", r.state_graph_theta);
        j.kv("n_components", r.state_n_components);
        j.kv("giant", r.state_giant);
        j.kv("giant_fraction", r.state_giant_fraction);
        j.end();
        j.key("bic_table");
        j.arr();
        for (const auto &info : r.bic_table) {
            j.obj();
            j.kv("n_states", info.n_states);
            j.kv("n_params", info.n_params);
            j.kv("loglik", info.loglik);
            j.kv("bic", info.bic);
            j.kv("aic", info.aic);
            j.kv("converged", info.converged != 0);
            j.end();
        }
        j.end();
        j.key("regimes");
        j.arr();
        for (const auto &s : r.stats) {
            j.obj();
            j.kv("state", s.state);
            j.kv("label", regime_label(s.state, n));
            j.kv("mean_daily", s.mean_daily);
            j.kv("ann_return", s.ann_return);
            j.kv("ann_vol", s.ann_vol);
            j.kv("sharpe", s.sharpe);
            j.kv("frac_days", s.frac_days);
            j.kv("self_transition", s.self_transition);
            j.kv("expected_duration", s.expected_duration);
            j.kv("observed_duration", s.observed_duration);
            j.kv("mean_vol", s.mean_vol);
            j.end();
        }
        j.end();
        j.key("percolation");
        j.obj();
        j.kv("enabled", cfg.use_percolation);
        j.kv("window", cfg.conn.window);
        j.kv("points", r.conn_points);
        j.kv("theta_c_mean", theta_mean);
        j.kv("theta_c_median", theta_med);
        j.kv("theta_c_last", std::isfinite(r.theta_c[T - 1]) ? r.theta_c[T - 1] : NAN);
        j.kv("gate_days", gate_days);
        j.kv("gate_fraction", gate_frac);
        j.end();
        j.key("backtest");
        j.obj();
        j.kv("oos_days", r.bt_hmm.n_days);
        j.key("buy_hold");
        json_metrics(j, r.bt_bh);
        j.key("hmm");
        json_metrics(j, r.bt_hmm);
        j.key("hmm_perc");
        json_metrics(j, r.bt_comb);
        j.end();
        j.key("manifest");
        j.obj();
        j.kv("seed", cfg.seed);
        j.kv("asset", cfg.asset);
        j.kv("states", n);
        j.kv("cost_bps", cfg.cost_bps);
        j.end();
        j.end();
        write_text(path_join(outdir, "metrics.json"), j.dump());

        /* report.md */
        std::ostringstream md;
        md << "# HMM × 渗流 · 状态与回测报告\n\n";
        md << "- 数据：`" << csv << "`，资产 `" << cfg.asset << "`，共 " << T << " 个交易日\n";
        md << "- 训练区间：" << panel->dates[r.start] << " → " << panel->dates[r.train_end]
           << "（" << r.train_end - r.start << " 日）\n";
        md << "- 样本外区间：" << panel->dates[r.oos_start] << " → " << panel->dates[T - 1]
           << "（" << T - r.oos_start << " 日），滚动重估每 " << cfg.refit_every << " 日\n";
        md << "- 模型：" << n << " 状态、混合数 " << cfg.n_mix
           << "、特征 = 日收益(%) + log 年化滚动波动(" << cfg.vol_window << "日)\n";
        md << "- 成本：" << cfg.cost_bps << " bps / 换手；状态按特征均值升序排列（0 = 最弱）\n\n";
        md << "## 模型选择\n\n";
        if (r.bic_table.size() > 1) {
            md << "| 状态数 | 参数数 | logL | BIC | 收敛 |\n| --- | --- | --- | --- | --- |\n";
            for (const auto &info : r.bic_table)
                md << "| " << info.n_states << " | " << info.n_params << " | "
                   << nfn(info.loglik, 2) << " | " << nfn(info.bic, 2) << " | "
                   << (info.converged ? "是" : "否") << " |\n";
            md << "\n按 BIC 选择 " << n << " 个状态。\n\n";
        } else {
            md << "固定 " << n << " 个状态；logL = " << nfn(r.report.loglik, 2)
               << "，BIC = " << nfn(r.report.bic, 2) << "。\n\n";
        }
        md << "## 状态画像\n\n";
        md << "| 状态 | 标签 | 日均收益 | 年化收益 | 年化波动 | Sharpe | 天数占比 | 平均持续(日) | 自持概率 |\n";
        md << "| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n";
        for (const auto &s : r.stats)
            md << "| " << s.state << " | " << regime_label(s.state, n) << " | "
               << pct(s.mean_daily) << " | " << pct(s.ann_return) << " | " << pct(s.ann_vol)
               << " | " << nfn(s.sharpe, 2) << " | " << pct(s.frac_days) << " | "
               << nfn(s.observed_duration, 1) << " | " << nfn(s.self_transition, 3)
               << " |\n";
        md << "\n转移矩阵：\n\n| from \\ to |";
        for (size_t k = 0; k < n; ++k) md << " " << regime_label(k, n) << " |";
        md << "\n| --- |";
        for (size_t k = 0; k < n; ++k) md << " --- |";
        md << "\n";
        for (size_t i = 0; i < n; ++i) {
            md << "| " << regime_label(i, n) << " |";
            for (size_t k = 0; k < n; ++k)
                md << " " << nfn(r.model_train->a[i * n + k], 3) << " |";
            md << "\n";
        }
        md << "\n状态图诊断（转移矩阵按其中位非对角转移概率 theta = " << nfn(r.state_graph_theta, 4)
           << " 阈值化）：" << r.state_n_components << " 个连通分量，巨分量 "
           << r.state_giant << "/" << n << "（" << pct(r.state_giant_fraction) << "）。\n";
        md << "\n## 最新状态\n\n";
        size_t last = T - 1;
        size_t best = 0;
        for (size_t k = 1; k < n; ++k)
            if (r.filtered[last * n + k] > r.filtered[last * n + best]) best = k;
        md << "- 日期：" << panel->dates[last] << "\n";
        md << "- 滤波状态：" << regime_label(best, n) << "（p = "
           << nfn(r.filtered[last * n + best], 3) << "）\n";
        md << "- 目标仓位：HMM " << pct(r.position_hmm[last]) << "；叠加渗流闸门后 "
           << pct(r.position_comb[last]) << "\n";
        if (std::isfinite(r.theta_c[last]))
            md << "- 渗流临界阈值 theta_c：" << nfn(r.theta_c[last], 3)
               << "（越高表示横截面越难被打散）\n";
        md << "\n## 样本外回测\n\n";
        md << "| 策略 | 总收益 | CAGR | 年化波动 | Sharpe | 最大回撤 | Calmar | 平均换手 | 平均仓位 | 胜率 |\n";
        md << "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n";
        md << metrics_row("买入持有", r.bt_bh);
        md << metrics_row("HMM 状态", r.bt_hmm);
        md << metrics_row("HMM + 渗流闸门", r.bt_comb);
        md << "\n信号在 t 日收盘生成、作用于 t+1 日收益；所有指标只统计样本外。\n\n";
        md << "## 渗流结构\n\n";
        md << "- 窗口 " << cfg.conn.window << " 日，相关阈值扫描 [" << cfg.conn.theta_min
           << ", " << cfg.conn.theta_max << "]，共 " << cfg.conn.theta_steps << " 个点\n";
        md << "- theta_c 中位数 " << nfn(theta_med, 3) << "，样本外均值 " << nfn(theta_mean, 3)
           << "，闸门触发天数占比 " << pct(gate_frac) << "\n";
        md << "- 闸门规则：theta_c 高于过去 " << (int)(100 * cfg.risk_quantile)
           << "% 分位时仓位乘以 " << nfn(cfg.risk_gate, 2) << "（至少 "
           << cfg.risk_min_obs << " 个历史点后才生效）\n\n";
        md << "## 读法与边界\n\n";
        md << "- 状态是统计标签，不是基本面判断；HMM 的切换通常滞后于价格。\n";
        md << "- 渗流阈值是横截面相关结构的描述量，在资产数很少时粒度过粗。\n";
        md << "- 合成数据上的结果用来验证算法可复现，不代表真实市场收益。\n";
        write_text(path_join(outdir, "report.md"), md.str());

        std::printf("regime: %zu states, train %zu days, OOS %zu days, seed %llu\n", n,
                    r.train_end - r.start, T - r.oos_start,
                    (unsigned long long)cfg.seed);
        for (const auto &s : r.stats)
            std::printf("  state %zu %-6s  ann_ret %7s  ann_vol %7s  days %6s  dur %5s\n",
                        s.state, regime_label(s.state, n).c_str(), pct(s.ann_return).c_str(),
                        pct(s.ann_vol).c_str(), pct(s.frac_days).c_str(),
                        nfn(s.observed_duration, 1));
        std::printf("  OOS  buy&hold: ret %7s  Sharpe %5s  MDD %7s\n",
                    pct(r.bt_bh.total_return).c_str(), nfn(r.bt_bh.sharpe, 2),
                    pct(r.bt_bh.max_drawdown).c_str());
        std::printf("       HMM      : ret %7s  Sharpe %5s  MDD %7s\n",
                    pct(r.bt_hmm.total_return).c_str(), nfn(r.bt_hmm.sharpe, 2),
                    pct(r.bt_hmm.max_drawdown).c_str());
        std::printf("       HMM+perc : ret %7s  Sharpe %5s  MDD %7s\n",
                    pct(r.bt_comb.total_return).c_str(), nfn(r.bt_comb.sharpe, 2),
                    pct(r.bt_comb.max_drawdown).c_str());
        std::printf("  wrote %s/{states,equity,bic}.csv, metrics.json, report.md, model_*.txt\n",
                    outdir.c_str());
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "hp_regime: %s\n", e.what());
        return 1;
    }
}
