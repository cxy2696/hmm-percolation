#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "hp/app/args.hpp"
#include "hp/app/connectedness.hpp"
#include "hp/app/report.hpp"
#include "hp/app/util.hpp"
#include "hp/core/csv.h"
#include "hp/core/stats.h"

using namespace hp::app;

int main(int argc, char **argv) {
    try {
        Args args(argc, argv);
        std::string csv = args.require("csv");
        ConnectednessConfig cfg;
        cfg.window = (size_t)args.get_int("window", (long)cfg.window);
        cfg.step = (size_t)args.get_int("step", (long)cfg.step);
        cfg.theta_min = args.get_double("theta-min", cfg.theta_min);
        cfg.theta_max = args.get_double("theta-max", cfg.theta_max);
        cfg.theta_ref = args.get_double("theta-ref", cfg.theta_ref);
        cfg.theta_steps = (size_t)args.get_int("theta-steps", (long)cfg.theta_steps);
        std::string outdir = args.get("outdir", "out/connectedness");
        ensure_dir(outdir);

        PanelPtr panel = load_panel(csv);
        Connectedness conn(cfg);
        ConnectednessResult r = conn.run(panel.get());

        std::vector<double> finite;
        for (double v : r.theta_c)
            if (std::isfinite(v)) finite.push_back(v);
        if (finite.empty()) throw std::runtime_error("no connectedness estimate produced");

        std::string path = path_join(outdir, "connectedness.csv");
        const char *names[7] = {"date", "theta_c", "theta_chi", "chi_max",
                                "giant_ref", "mean_degree_ref", "n_edges_ref"};
        hp_csv_writer *w = hp_csv_writer_open(path.c_str(), names, 7, 1);
        if (!w) throw std::runtime_error("cannot write " + path);
        for (const auto &pt : r.points) {
            double vals[6] = {pt.theta_c, pt.theta_chi, pt.chi_max, pt.giant_ref,
                              pt.mean_degree_ref, (double)pt.n_edges_ref};
            hp_status st = hp_csv_writer_row_label(w, panel->dates[pt.t], vals, 6);
            if (st != HP_OK) {
                hp_csv_writer_abort(w);
                throw std::runtime_error("write connectedness.csv");
            }
        }
        if (hp_csv_writer_close(w) != HP_OK) throw std::runtime_error("close connectedness.csv");

        /* Correlation snapshot at the last sampled window. */
        if (!r.last_corr.empty()) {
            std::string cpath = path_join(outdir, "corr_last.csv");
            hp_csv_writer *cw = hp_csv_writer_open(cpath.c_str(), (const char *const *)panel->symbols, panel->N, 0);
            if (cw) {
                for (size_t i = 0; i < panel->N; ++i)
                    hp_csv_writer_row(cw, r.last_corr.data() + i * panel->N, panel->N);
                hp_csv_writer_close(cw);
            }
        }

        std::sort(finite.begin(), finite.end());
        double mean = hp_mean(finite.data(), finite.size());
        double median = hp_quantile_sorted(finite.data(), finite.size(), 0.5);
        double q80 = hp_quantile_sorted(finite.data(), finite.size(), 0.8);
        size_t max_t = 0;
        double max_v = -1.0;
        for (size_t t = 0; t < r.T; ++t) {
            if (std::isfinite(r.theta_c[t]) && r.theta_c[t] > max_v) {
                max_v = r.theta_c[t];
                max_t = t;
            }
        }
        size_t above = 0;
        for (double v : finite)
            if (v > q80) ++above;

        Json j;
        j.obj();
        j.kv("tool", "hp_connectedness");
        j.kv("csv", csv);
        j.kv("T", r.T);
        j.kv("N", r.N);
        j.kv("window", cfg.window);
        j.kv("step", cfg.step);
        j.kv("theta_min", cfg.theta_min);
        j.kv("theta_max", cfg.theta_max);
        j.kv("theta_steps", cfg.theta_steps);
        j.kv("theta_ref", cfg.theta_ref);
        j.kv("n_points", r.points.size());
        j.kv("theta_c_mean", mean);
        j.kv("theta_c_median", median);
        j.kv("theta_c_q80", q80);
        j.kv("theta_c_min", finite.front());
        j.kv("theta_c_max", finite.back());
        j.kv("theta_c_max_date", panel->dates[max_t]);
        j.kv("theta_c_last", r.theta_c[r.T - 1]);
        j.kv("frac_above_q80", (double)above / (double)finite.size());
        j.end();
        write_text(path_join(outdir, "summary.json"), j.dump());

        std::printf("connectedness: %zu points, theta_c mean %.4f median %.4f max %.4f (%s)\n",
                    r.points.size(), mean, median, finite.back(), panel->dates[max_t]);
        std::printf("  last theta_c %.4f, giant@ref %.3f\n", r.theta_c[r.T - 1],
                    r.giant_ref[r.T - 1]);
        std::printf("  wrote %s/summary.json, connectedness.csv, corr_last.csv\n", outdir.c_str());
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "hp_connectedness: %s\n", e.what());
        return 1;
    }
}
