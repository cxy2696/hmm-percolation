# Strategy comparison

First decision: 2018-12-06. Signals: 1192 rows; P&L days: 1191. Cost: 5.0 bps per unit turnover. All rows use signal[t] on return[t+1].

| strategy | total | CAGR | vol | Sharpe | max DD | turnover | exposure |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| buy_hold | 92.77% | 14.90% | 18.34% | 0.85 | -33.83% | 0.001 | 100.00% |
| constant | 41.62% | 7.64% | 9.17% | 0.85 | -18.27% | 0.000 | 50.00% |
| vol_target | 56.89% | 10.00% | 11.51% | 0.89 | -17.96% | 0.027 | 63.49% |
| momentum | 37.01% | 6.89% | 14.34% | 0.54 | -24.60% | 0.056 | 70.28% |
| gate_only | 102.59% | 16.11% | 16.92% | 0.97 | -25.64% | 0.009 | 94.00% |
| hmm | 53.40% | 9.48% | 12.81% | 0.77 | -22.40% | 0.030 | 66.34% |
| hmm_gate | 60.65% | 10.55% | 12.19% | 0.88 | -16.90% | 0.032 | 63.28% |
| hmm_vol | 33.35% | 6.28% | 7.94% | 0.81 | -12.89% | 0.024 | 42.39% |
| hmm_vol_gate | 36.27% | 6.77% | 7.66% | 0.89 | -9.81% | 0.025 | 40.79% |

The HMM rows are the only ones that use a fitted state model; `buy_hold`, `constant`, `vol_target` and `momentum` are there to make the comparison falsifiable. A higher Sharpe on one synthetic sample is not evidence of a durable edge.
