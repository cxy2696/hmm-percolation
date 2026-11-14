#include "hp/app/connectedness.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "hp/core/stats.h"
#include "hp/perc/corr.h"
#include "hp/perc/graph.h"

namespace hp::app {

Connectedness::Connectedness(ConnectednessConfig cfg) : cfg_(cfg) {
    if (cfg_.window < 10) throw std::runtime_error("connectedness window must be >= 10");
    if (cfg_.step < 1) cfg_.step = 1;
    if (!(cfg_.theta_min >= 0.0 && cfg_.theta_min < cfg_.theta_max && cfg_.theta_max <= 1.0))
        throw std::runtime_error("invalid theta range");
    if (cfg_.theta_steps < 3) throw std::runtime_error("theta_steps must be >= 3");
}

ConnectednessResult Connectedness::run(const hp_panel *panel) const {
    if (!panel || panel->N < 2) throw std::runtime_error("connectedness needs at least 2 assets");
    const size_t T = panel->T, N = panel->N;
    ConnectednessResult out;
    out.T = T;
    out.N = N;
    out.window = cfg_.window;
    out.theta_c.assign(T, std::numeric_limits<double>::quiet_NaN());
    out.theta_chi.assign(T, std::numeric_limits<double>::quiet_NaN());
    out.giant_ref.assign(T, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> corr(N * N);
    for (size_t t = cfg_.window; t < T; t += cfg_.step) {
        const double *rows = panel->ret + (t - cfg_.window + 1) * N;
        hp_corr_matrix(rows, cfg_.window, N, corr.data());
        hp_corr_critical_result crit;
        hp_status st = hp_corr_critical(corr.data(), N, cfg_.theta_min, cfg_.theta_max,
                                        cfg_.theta_steps, &crit);
        if (st != HP_OK) continue;
        hp_graph *g = hp_graph_from_corr(corr.data(), N, cfg_.theta_ref);
        if (!g) continue;
        hp_components comp;
        if (hp_graph_components(g, -INFINITY, &comp) != HP_OK) {
            hp_graph_free(g);
            continue;
        }
        ConnectednessPoint pt;
        pt.t = t;
        pt.theta_c = crit.theta_half;
        pt.theta_chi = crit.theta_c;
        pt.chi_max = crit.chi_max;
        pt.giant_ref = comp.giant_fraction;
        pt.n_edges_ref = hp_graph_edge_count(g);
        pt.mean_degree_ref = 2.0 * (double)pt.n_edges_ref / (double)N;
        out.points.push_back(pt);
        out.theta_c[t] = pt.theta_c;
        out.theta_chi[t] = pt.theta_chi;
        out.giant_ref[t] = pt.giant_ref;
        out.last_corr = corr;
        hp_components_free(&comp);
        hp_graph_free(g);
    }
    return out;
}

}  // namespace hp::app
