# 策略对照

原先的结果页只有 buy & hold、HMM、HMM + 渗流闸门三条线。这样很容易把“回撤比 buy & hold 小”误读成 HMM 的功劳，因为仓位本身变了。这里加三个不依赖 HMM 的基线，用同一段样本外、同一笔成本、同一个执行滞后做对照。

## 规则

所有信号都满足 `signal[t]` 作用于 `return[t+1]`，只用 $t$ 及以前的数据。

| 名称 | 仓位 | 说明 |
| --- | --- | --- |
| `buy_hold` | $1$ | 满仓，基准 |
| `constant` | $0.5$ | 固定半仓，用来区分“降风险”和“择时” |
| `vol_target` | $\min(1,0.10/\hat\sigma_t)$ | 20 日年化波动，目标 10%，不加杠杆 |
| `momentum` | $1$ 或 $0$ | 过去 60 日累计对数收益为正则持有，否则空仓 |
| `gate_only` | `gate` | 只用渗流闸门：正常 1，触发时 0.5 |
| `hmm` | $0.1,0.55,1.0$ | 三状态 GMM-HMM 的滤波后验映射，按均值排序 |
| `hmm_gate` | `hmm × gate` | 相关图 `theta_c` 的 80% 分位闸门 |
| `hmm_vol` | $\min(1,hmm\times 0.10/\hat\sigma_t)$ | HMM × 波动目标 |
| `hmm_vol_gate` | $\min(1,hmm\times 0.10/\hat\sigma_t\times gate)$ | 三者相乘 |

`vol_target` 的波动用截至 $t$ 的 20 个收益计算；不足 20 天时信号为 0（空仓），不是 1。`momentum` 的窗口含 $t$，同样只用历史。

## 结果

合成面板，样本外第一次决策在 2018-12-06 收盘；有效 P&L 1191 天；单边换手 5 bps。

| 策略 | 总收益 | CAGR | 年化波动 | Sharpe | 最大回撤 |
| --- | ---: | ---: | ---: | ---: | ---: |
| buy_hold | 92.77% | 14.90% | 18.34% | 0.85 | -33.83% |
| constant | 41.62% | 7.64% | 9.17% | 0.85 | -18.27% |
| vol_target | 56.89% | 10.00% | 11.51% | 0.89 | -17.96% |
| momentum | 37.01% | 6.89% | 14.34% | 0.54 | -24.60% |
| gate_only | 102.59% | 16.11% | 16.92% | 0.97 | -25.64% |
| hmm | 53.40% | 9.48% | 12.81% | 0.77 | -22.40% |
| hmm_gate | 60.65% | 10.55% | 12.19% | 0.88 | -16.90% |
| hmm_vol | 33.35% | 6.28% | 7.94% | 0.81 | -12.89% |
| hmm_vol_gate | 36.27% | 6.77% | 7.66% | 0.89 | -9.81% |

## 怎么读

- `vol_target` 的 Sharpe 略高于 buy & hold，但这是**风险缩放**，不是方向 alpha。它的波动没有正好落在 10%，因为 20 日窗口、1 倍上限和离散调仓都有误差。
- `momentum` 在这份合成数据上只有 0.54 的 Sharpe，胜率 38.5%。它没有提供择时价值；这比“只展示表现好的策略”更有用。
- **单独闸门 `gate_only` 的 Sharpe 0.97，高于 `hmm_gate` 的 0.88。** 这是最重要的负结果：在这份合成数据上，收益主要来自“相关结构变差时减仓”，HMM 的状态判断反而拖累了它。不能把 `hmm_gate` 的改善归功于 HMM。
- `hmm_gate` 仍比单独的 `hmm` 好，但没有超过 `vol_target` 的 Sharpe。闸门是风险 overlay，不是独立的收益来源。
- `hmm_vol_gate` 的回撤最小，但绝对收益也最低。它可以作为“低风险版本”的候选，不能因此说它“最优”。
- 所有数字来自一条合成路径，不能当作样本外显著性或实盘预期。策略之间高度相关，差异很可能来自少数几个时段。

## 验证

- `tests/test_strategy.c` 检查滚动波动、动量、状态映射、NaN 传播和“篡改未来收益不改变此前信号”。
- `tests/test_app.cpp::test_strategy_comparison` 检查九条线共享同一 OOS 切片，并核对 buy & hold / HMM+gate 与 `run_regime` 的旧指标一致。
- Python、Julia 绑定有同一组信号的数值测试，避免两套绑定各写一套规则。
- 复现：

```bash
make -j
./build/hp_strategies --csv data/synthetic_panel.csv --outdir out/strategies --seed 42 --states 3
```

输出 `strategies.csv`、`equity.csv`、`positions.csv`、`summary.json`、`report.md`。

## 加自己的策略

先在 `src/market/strategy.c` 写一个纯函数，或在 `src/app/strategies.cpp` 组合已有信号；然后加一个“只用过去”的测试。具体落点见 [extending.md](extending.md)。

不要在这份样本上反复调 lookback，再把最好的结果写进 README。那会把“可复现的对照”变成“事后选参数”。
