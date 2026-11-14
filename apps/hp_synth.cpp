#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "hp/app/args.hpp"
#include "hp/app/util.hpp"
#include "hp/market/synth.h"

using namespace hp::app;

int main(int argc, char **argv) {
    try {
        Args args(argc, argv);
        hp_synth_opts o;
        hp_synth_opts_default(&o);
        o.n_days = (size_t)args.get_int("days", (long)o.n_days);
        o.n_assets = (size_t)args.get_int("assets", (long)o.n_assets);
        o.n_regimes = (size_t)args.get_int("regimes", (long)o.n_regimes);
        o.seed = (uint64_t)args.get_int("seed", (long)o.seed);
        std::string out = args.get("out", "data/synthetic_panel.csv");
        std::string truth = args.get("truth", "data/synthetic_truth.csv");
        if (o.n_days < 100 || o.n_assets < 2 || (o.n_regimes != 2 && o.n_regimes != 3))
            throw std::runtime_error("require --days >= 100, --assets >= 2, --regimes 2 or 3");
        hp_synth_data *d = hp_synth_generate(&o);
        if (!d) throw std::runtime_error("hp_synth_generate failed");
        ensure_parent(out);
        ensure_parent(truth);
        hp_status st = hp_panel_save_csv(d->panel, out.c_str());
        if (st != HP_OK) throw std::runtime_error(std::string("save panel: ") + hp_status_message(st));
        st = hp_synth_write_truth(d, truth.c_str());
        if (st != HP_OK) throw std::runtime_error(std::string("save truth: ") + hp_status_message(st));

        std::printf("synthetic market: %zu days x %zu assets, %zu regimes, seed %llu\n",
                    d->panel->T, d->panel->N, d->n_regimes,
                    (unsigned long long)o.seed);
        std::vector<size_t> cnt(d->n_regimes, 0);
        std::vector<double> sum(d->n_regimes, 0.0);
        for (size_t t = 0; t < d->panel->T; ++t) {
            cnt[d->regime[t]]++;
            sum[d->regime[t]] += d->market_ret[t];
        }
        for (size_t r = 0; r < d->n_regimes; ++r)
            std::printf("  regime %zu: %5.1f%% of days, mean daily market %.5f\n", r,
                        100.0 * (double)cnt[r] / (double)d->panel->T,
                        cnt[r] ? sum[r] / (double)cnt[r] : 0.0);
        std::printf("  panel: %s\n  truth: %s\n", out.c_str(), truth.c_str());
        hp_synth_free(d);
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "hp_synth: %s\n", e.what());
        return 1;
    }
}
