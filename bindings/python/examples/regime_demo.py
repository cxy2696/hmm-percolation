"""End-to-end Python example: EW index -> Gaussian HMM -> percolation gate.

Run from the repository root after `make shared`:

    PYTHONPATH=bindings/python python3 bindings/python/examples/regime_demo.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "bindings" / "python"))

import hpmarket as hp


def build_features(ret, window=20):
    x, start = [], None
    for t in range(len(ret)):
        if t < window:
            x.append([math.nan, math.nan])
            continue
        win = ret[t - window + 1:t + 1]
        mean = sum(win) / window
        var = sum((v - mean) ** 2 for v in win) / (window - 1)
        vol = max(math.sqrt(var * 252.0), 1e-4)
        x.append([ret[t] * 100.0, math.log(vol * 100.0)])
        if start is None:
            start = t
    return x, start


def main():
    panel = hp.Panel.load(ROOT / "data" / "synthetic_panel.csv")
    ret = panel.ew_returns()
    x, start = build_features(ret, 20)
    train_end = start + int((len(ret) - start) * 0.6)
    mean = [sum(x[t][d] for t in range(start, train_end)) / (train_end - start) for d in range(2)]
    sd = []
    for d in range(2):
        var = sum((x[t][d] - mean[d]) ** 2 for t in range(start, train_end)) / (train_end - start - 1)
        sd.append(math.sqrt(var))
    z = [[(x[t][d] - mean[d]) / sd[d] for d in range(2)] for t in range(start, len(ret))]

    model = hp.GaussianHMM(3, 2)
    report = model.train(z[: train_end - start], max_iter=150, restarts=4, seed=42)
    print("training:", report)
    post, ll = model.posterior(z)
    path, _ = model.viterbi(z)
    filt, _ = model.filtered(z)
    print(f"full-sample loglik {ll:.3f}; states 0..2 means "
          f"{[round(model.state_mean(j), 3) for j in range(3)]}")

    exposure = [0.1, 0.55, 1.0]
    pos = [sum(filt[t - start][j] * exposure[j] for j in range(3))
           for t in range(start, len(ret))]
    oos = train_end - start
    metrics = hp.backtest(ret[train_end - 1:], pos[oos - 1:], cost_bps=5.0)
    bh = hp.backtest(ret[train_end - 1:], [1.0] * (len(ret) - train_end + 1), cost_bps=5.0)
    print(f"OOS  buy&hold: return {bh['total_return']:.2%}  Sharpe {bh['sharpe']:.2f}  "
          f"MDD {bh['max_drawdown']:.2%}")
    print(f"     HMM     : return {metrics['total_return']:.2%}  Sharpe {metrics['sharpe']:.2f}  "
          f"MDD {metrics['max_drawdown']:.2%}")

    w = 60
    window = [[panel.return_series(j)[t] for j in range(panel.N)]
              for t in range(panel.T - w, panel.T)]
    corr = hp.corr_matrix(window)
    crit = hp.corr_critical(corr)
    giant = hp.graph_giant(corr, 0.90)
    print(f"last window: theta_c {crit['theta_half']:.3f}, "
          f"giant@0.90 {giant['giant_fraction']:.2f}, edges {giant['n_edges']}")


if __name__ == "__main__":
    main()
