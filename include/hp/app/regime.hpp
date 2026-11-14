/* HP app · regime pipeline: features -> Gaussian HMM -> walk-forward
 * filtered probabilities -> exposure -> optional percolation risk gate ->
 * out-of-sample backtest.
 *
 * Data flow is strictly causal after `train_end`: the model used on day t was
 * fitted on data before its refit block, and the filtered state probability
 * uses observations up to t only.
 */
#ifndef HP_APP_REGIME_HPP
#define HP_APP_REGIME_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "hp/app/connectedness.hpp"
#include "hp/app/util.hpp"
#include "hp/hmm/ghmm.h"
#include "hp/market/backtest.h"
#include "hp/market/panel.h"

namespace hp::app {

struct RegimeConfig {
    std::string asset = "EW";
    size_t states = 3; /* 0 = BIC select 2..5 */
    size_t n_mix = 1;
    size_t vol_window = 20;
    double train_frac = 0.6;
    size_t refit_every = 250;
    unsigned max_iter = 150;
    unsigned restarts = 4;
    uint64_t seed = 42;
    double cost_bps = 5.0;

    bool use_percolation = true;
    ConnectednessConfig conn;

    double risk_quantile = 0.80;
    size_t risk_lookback = 500;
    size_t risk_min_obs = 120;
    double risk_gate = 0.5;
};

struct RegimeStat {
    size_t state = 0;
    double mean_daily = 0.0;
    double ann_return = 0.0;
    double ann_vol = 0.0;
    double sharpe = 0.0;
    double frac_days = 0.0;
    double self_transition = 0.0;
    double expected_duration = 0.0;
    double observed_duration = 0.0;
    double mean_vol = 0.0;
};

struct RegimeResult {
    size_t n_states = 0;
    size_t start = 0;
    size_t train_end = 0;
    size_t oos_start = 0;
    std::vector<double> ret;       /* T, asset/EW log return */
    std::vector<double> vol;       /* T, annualised rolling vol */
    std::vector<double> posterior; /* T*n, smoothed */
    std::vector<double> filtered;  /* T*n, causal; OOS rows from refits */
    std::vector<uint32_t> viterbi; /* T, UINT32_MAX before start */
    std::vector<double> position_hmm;
    std::vector<double> gate;
    std::vector<double> position_comb;
    std::vector<double> theta_c;
    std::vector<double> giant_ref;
    std::vector<double> feature_mean;
    std::vector<double> feature_sd;

    GhmmPtr model_train{nullptr, hp_ghmm_free};
    GhmmPtr model_final{nullptr, hp_ghmm_free};
    hp_ghmm_report report{};
    std::vector<hp_ghmm_model_info> bic_table;
    std::vector<RegimeStat> stats;
    size_t conn_points = 0;
    /* Transition matrix as a threshold graph (doc 4.3 diagnostic). */
    double state_graph_theta = 0.0;
    size_t state_n_components = 0;
    size_t state_giant = 0;
    double state_giant_fraction = 0.0;

    hp_bt_metrics bt_bh{};
    hp_bt_metrics bt_hmm{};
    hp_bt_metrics bt_comb{};
    std::vector<double> eq_bh;
    std::vector<double> eq_hmm;
    std::vector<double> eq_comb;
    std::vector<double> bt_ret;
    std::vector<std::string> bt_dates;
};

RegimeResult run_regime(const hp_panel *panel, const RegimeConfig &cfg);

std::string regime_label(size_t state, size_t n_states);

}  // namespace hp::app

#endif
