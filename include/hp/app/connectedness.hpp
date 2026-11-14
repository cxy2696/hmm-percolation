/* HP app · rolling cross-sectional connectedness via correlation percolation. */
#ifndef HP_APP_CONNECTEDNESS_HPP
#define HP_APP_CONNECTEDNESS_HPP

#include <cstddef>
#include <vector>

#include "hp/market/panel.h"

namespace hp::app {

struct ConnectednessConfig {
    size_t window = 60;
    size_t step = 1;
    double theta_min = 0.10;
    double theta_max = 0.99;
    double theta_ref = 0.90;
    size_t theta_steps = 45;
};

struct ConnectednessPoint {
    size_t t = 0;
    double theta_c = 0.0;   /* giant-fraction 0.5 crossing */
    double theta_chi = 0.0; /* susceptibility peak */
    double chi_max = 0.0;
    double giant_ref = 0.0;
    double mean_degree_ref = 0.0;
    size_t n_edges_ref = 0;
};

struct ConnectednessResult {
    size_t T = 0;
    size_t N = 0;
    size_t window = 0;
    std::vector<double> theta_c;   /* T, NAN before the first window */
    std::vector<double> theta_chi; /* T, susceptibility peak */
    std::vector<double> giant_ref; /* T */
    std::vector<ConnectednessPoint> points;
    std::vector<double> last_corr;  /* N*N at the last sampled date */
};

class Connectedness {
public:
    explicit Connectedness(ConnectednessConfig cfg);
    ConnectednessResult run(const hp_panel *panel) const;

private:
    ConnectednessConfig cfg_;
};

}  // namespace hp::app

#endif
