/* HP app · compare several causal position rules on the same OOS slice.
 *
 * This is deliberately a flat table, not a strategy framework: each row is a
 * pure function of past data, and all rows share dates, costs and the one-day
 * execution lag.  Add a rule here only after adding a test for its timing.
 */
#ifndef HP_APP_STRATEGIES_HPP
#define HP_APP_STRATEGIES_HPP

#include <string>
#include <vector>

#include "hp/app/regime.hpp"
#include "hp/market/backtest.h"

namespace hp::app {

struct StrategyConfig {
    double constant_exposure = 0.5;
    size_t vol_lookback = 20;
    double target_vol = 0.10;
    double min_position = 0.0;
    double max_position = 1.0;
    size_t momentum_lookback = 60;
    double momentum_threshold = 0.0;
    double cost_bps = 5.0;
};

struct StrategyResult {
    std::string name;
    hp_bt_metrics metrics{};
    std::vector<double> equity;
    std::vector<double> position;
};

struct StrategyComparison {
    size_t oos_start = 0;
    std::vector<std::string> dates;
    std::vector<double> ret;
    std::vector<StrategyResult> rows;
};

StrategyComparison compare_strategies(const RegimeResult &r,
                                      const StrategyConfig &cfg);

}  // namespace hp::app

#endif /* HP_APP_STRATEGIES_HPP */
