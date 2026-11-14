#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "hp/app/args.hpp"
#include "hp/app/report.hpp"
#include "hp/app/util.hpp"
#include "hp/core/csv.h"
#include "hp/core/stats.h"
#include "hp/perc/lattice.h"
#include "hp/perc/random_graph.h"

using namespace hp::app;

namespace {

hp_perc_kind parse_kind(const std::string &s) {
    if (s == "bond") return HP_PERC_BOND;
    if (s == "site") return HP_PERC_SITE;
    throw std::runtime_error("--kind must be bond or site");
}

struct Summary {
    double s = 0.0;
    double span = 0.0;
    double chi = 0.0;
    double clusters = 0.0;
    double occupancy = 0.0;
};

Summary run_summary(hp_perc_kind kind, int L, double p, size_t reps, uint64_t seed,
                    int periodic) {
    hp_lattice *lat = hp_lattice_new(L, kind, periodic);
    if (!lat) throw std::runtime_error("bad lattice size");
    hp_rng rng;
    hp_rng_seed_stream(&rng, seed, (uint64_t)L * 131u + (uint64_t)kind);
    Summary out;
    for (size_t r = 0; r < reps; ++r) {
        if (hp_lattice_bernoulli(lat, p, &rng) != HP_OK) throw std::runtime_error("bernoulli");
        hp_lattice_stats st;
        if (hp_lattice_analyze(lat, &st) != HP_OK) throw std::runtime_error("analyze");
        out.s += st.s_max;
        out.span += st.spanning ? 1.0 : 0.0;
        out.chi += st.susceptibility;
        out.clusters += (double)st.n_clusters;
        out.occupancy += st.occupancy;
    }
    out.s /= (double)reps;
    out.span /= (double)reps;
    out.chi /= (double)reps;
    out.clusters /= (double)reps;
    out.occupancy /= (double)reps;
    hp_lattice_free(lat);
    return out;
}

std::vector<int> parse_sizes(const std::string &s) {
    std::vector<int> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        int v = std::stoi(item);
        if (v < 8 || v > 2048) throw std::runtime_error("lattice size out of range [8,2048]");
        out.push_back(v);
    }
    if (out.empty()) throw std::runtime_error("empty --sizes");
    return out;
}

double known_pc(hp_perc_kind kind) {
    return kind == HP_PERC_BOND ? 0.5 : 0.592746;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        Args args(argc, argv);
        std::string mode = args.get("mode", "bond");
        std::string kind_s = args.get("kind", "bond");
        hp_perc_kind kind = parse_kind(kind_s);
        int L = (int)args.get_int("L", 64);
        int periodic = args.has("periodic") ? 1 : 0;
        size_t reps = (size_t)args.get_int("reps", 200);
        uint64_t seed = (uint64_t)args.get_int("seed", 1);
        if (reps == 0) throw std::runtime_error("--reps must be positive");

        if (mode == "bond" || mode == "site") {
            kind = parse_kind(mode);
            double p = args.get_double("p", 0.5);
            Summary s = run_summary(kind, L, p, reps, seed, periodic);
            Json j;
            j.obj();
            j.kv("mode", "single");
            j.kv("kind", hp_perc_kind_name(kind));
            j.kv("L", L);
            j.kv("p", p);
            j.kv("reps", reps);
            j.kv("periodic", periodic != 0);
            j.kv("s_max", s.s);
            j.kv("spanning", s.span);
            j.kv("susceptibility", s.chi);
            j.kv("clusters", s.clusters);
            j.kv("occupancy", s.occupancy);
            j.end();
            std::string jpath = args.get("json", "out/perc_single.json");
            ensure_parent(jpath);
            write_text(jpath, j.dump());
            std::printf("%s L=%d p=%.3f reps=%zu: S=%.4f span=%.3f chi=%.4f "
                        "clusters=%.2f\n",
                        hp_perc_kind_name(kind), L, p, reps, s.s, s.span, s.chi, s.clusters);
            return 0;
        }

        if (mode == "scan") {
            double pmin = args.get_double("pmin", 0.3);
            double pmax = args.get_double("pmax", 0.7);
            size_t steps = (size_t)args.get_int("steps", 21);
            if (!(pmin < pmax) || steps < 2) throw std::runtime_error("bad scan range");
            std::string out = args.get("out", "out/perc_scan.csv");
            ensure_parent(out);
            const char *names[6] = {"p", "s_max", "spanning", "susceptibility",
                                    "clusters", "occupancy"};
            hp_csv_writer *w = hp_csv_writer_open(out.c_str(), names, 6, 0);
            if (!w) throw std::runtime_error("cannot write " + out);
            for (size_t k = 0; k < steps; ++k) {
                double p = pmin + (pmax - pmin) * (double)k / (double)(steps - 1);
                Summary s = run_summary(kind, L, p, reps, seed + k, periodic);
                double vals[6] = {p, s.s, s.span, s.chi, s.clusters, s.occupancy};
                if (hp_csv_writer_row(w, vals, 6) != HP_OK) {
                    hp_csv_writer_abort(w);
                    throw std::runtime_error("write scan");
                }
            }
            if (hp_csv_writer_close(w) != HP_OK) throw std::runtime_error("close scan");
            std::printf("scan: %s L=%d %s p=[%.2f,%.2f] steps=%zu reps=%zu -> %s\n",
                        hp_perc_kind_name(kind), L, periodic ? "periodic" : "open",
                        pmin, pmax, steps, reps, out.c_str());
            return 0;
        }

        if (mode == "scaling") {
            std::vector<int> sizes = parse_sizes(args.get("sizes", "32,64,128,256"));
            unsigned iters = (unsigned)args.get_int("iters", 16);
            std::string out = args.get("out", "out/perc_scaling.csv");
            ensure_parent(out);
            std::vector<double> pcs(sizes.size()), sps(sizes.size());
            std::vector<double> x(sizes.size()), y(sizes.size());
            for (size_t k = 0; k < sizes.size(); ++k) {
                hp_perc_exp e;
                e.L = sizes[k];
                e.kind = kind;
                e.periodic = periodic;
                e.reps = reps;
                e.seed = seed + 17 * k;
                double pc = 0.0, sp = 0.0;
                double lo = kind == HP_PERC_BOND ? 0.35 : 0.45;
                double hi = kind == HP_PERC_BOND ? 0.65 : 0.75;
                if (hp_perc_crossing(&e, lo, hi, iters, &pc, &sp) != HP_OK)
                    throw std::runtime_error("crossing failed");
                pcs[k] = pc;
                sps[k] = sp;
                x[k] = std::pow((double)sizes[k], -0.75); /* 1/nu, nu=4/3 */
                y[k] = pc;
            }
            double slope = 0.0, intercept = 0.0;
            double r2 = hp_linreg(x.data(), y.data(), sizes.size(), &slope, &intercept);
            const char *names[3] = {"L", "p_cross", "span_at_cross"};
            hp_csv_writer *w = hp_csv_writer_open(out.c_str(), names, 3, 0);
            if (!w) throw std::runtime_error("cannot write " + out);
            for (size_t k = 0; k < sizes.size(); ++k) {
                double vals[3] = {(double)sizes[k], pcs[k], sps[k]};
                hp_csv_writer_row(w, vals, 3);
            }
            hp_csv_writer_close(w);
            std::string jpath = args.get("json", "out/perc_scaling.json");
            ensure_parent(jpath);
            Json j;
            j.obj();
            j.kv("mode", "scaling");
            j.kv("kind", hp_perc_kind_name(kind));
            j.kv("periodic", periodic != 0);
            j.kv("reps", reps);
            j.kv("iters", iters);
            j.kv("nu_assumed", 4.0 / 3.0);
            j.kv("pc_infinite", intercept);
            j.kv("fit_slope", slope);
            j.kv("fit_r2", r2);
            j.kv("pc_known", known_pc(kind));
            j.kv("pc_abs_error", std::fabs(intercept - known_pc(kind)));
            j.key("points");
            j.arr();
            for (size_t k = 0; k < sizes.size(); ++k) {
                j.obj();
                j.kv("L", (size_t)sizes[k]);
                j.kv("pc", pcs[k]);
                j.kv("span", sps[k]);
                j.end();
            }
            j.end();
            j.end();
            write_text(jpath, j.dump());
            std::printf("scaling %s: p_c(inf) = %.4f (known %.4f, err %.4f), R2=%.4f\n",
                        hp_perc_kind_name(kind), intercept, known_pc(kind),
                        std::fabs(intercept - known_pc(kind)), r2);
            for (size_t k = 0; k < sizes.size(); ++k)
                std::printf("  L=%4d  p_cross=%.4f\n", sizes[k], pcs[k]);
            std::printf("  wrote %s and %s\n", out.c_str(), jpath.c_str());
            return 0;
        }

        if (mode == "meanfield") {
            double cmin = args.get_double("cmin", 0.5);
            double cmax = args.get_double("cmax", 3.0);
            size_t steps = (size_t)args.get_int("steps", 26);
            if (!(cmin < cmax) || steps < 2) throw std::runtime_error("bad meanfield range");
            std::string out = args.get("out", "out/perc_meanfield.csv");
            ensure_parent(out);
            const char *names[3] = {"mean_degree", "giant_fraction", "bond_pc"};
            hp_csv_writer *w = hp_csv_writer_open(out.c_str(), names, 3, 0);
            if (!w) throw std::runtime_error("cannot write " + out);
            for (size_t k = 0; k < steps; ++k) {
                double c = cmin + (cmax - cmin) * (double)k / (double)(steps - 1);
                double vals[3] = {c, hp_perc_er_giant_theory(c), hp_perc_bond_threshold(c)};
                hp_csv_writer_row(w, vals, 3);
            }
            hp_csv_writer_close(w);
            std::printf("meanfield: S(c) for c in [%.2f,%.2f], wrote %s\n", cmin, cmax,
                        out.c_str());
            return 0;
        }

        throw std::runtime_error("unknown --mode (bond|site|scan|scaling|meanfield)");
    } catch (const std::exception &e) {
        std::fprintf(stderr, "hp_perc_lab: %s\n", e.what());
        return 1;
    }
}
