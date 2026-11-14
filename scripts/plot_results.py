#!/usr/bin/env python3
"""Generate the figures used by the documentation.

Run from the project root after `scripts/reproduce.sh`:

    MPLCONFIGDIR=build/mplcache build/venv/bin/python scripts/plot_results.py
"""
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def col(rows, key, default=float("nan")):
    out = []
    for r in rows:
        v = r.get(key, "")
        try:
            out.append(float(v) if v not in ("", "nan", "NA") else np.nan)
        except ValueError:
            out.append(default)
    return np.array(out)


def dates(rows):
    return [datetime.strptime(r["date"], "%Y-%m-%d") for r in rows]


def shade_bear(ax, dts, regime):
    start = None
    for i, r in enumerate(regime):
        if r == 2 and start is None:
            start = dts[i]
        elif r != 2 and start is not None:
            ax.axvspan(start, dts[i], color="0.85", alpha=0.55, lw=0)
            start = None
    if start is not None:
        ax.axvspan(start, dts[-1], color="0.85", alpha=0.55, lw=0)


def fig_regimes():
    states = read_csv(ROOT / "out/regime/states.csv")
    truth = read_csv(ROOT / "data/synthetic_truth.csv")
    meta = json.loads((ROOT / "out/regime/metrics.json").read_text())
    dts = dates(states)
    ret = col(states, "ret")
    p0 = col(states, "p0")
    theta = col(states, "theta_c")
    pos_hmm = col(states, "position_hmm")
    pos_comb = col(states, "position_comb")
    gate = col(states, "gate")
    regime = [int(float(r["regime"])) for r in truth]

    fig, axes = plt.subplots(3, 1, figsize=(11, 8.2), sharex=True,
                             gridspec_kw={"height_ratios": [2.2, 1.6, 1.4]})
    ax = axes[0]
    ax.plot(dts, np.nancumsum(np.nan_to_num(ret)) * 100, color="#1f4e79", lw=1.4,
            label="equal-weight index (cumulative log return, %)")
    shade_bear(ax, dts, regime)
    ax.set_ylabel("cum. return (%)")
    ax.legend(loc="upper left", frameon=False)
    ax.set_title("HMM regimes and percolation gate on the synthetic market")

    ax = axes[1]
    ax.plot(dts, p0, color="#b03a2e", lw=1.0, label="filtered P(weak state)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("P(weak)")
    ax2 = ax.twinx()
    ax2.plot(dts, theta, color="#117a65", lw=0.9, alpha=0.85, label="theta_c")
    ax2.set_ylabel("theta_c", color="#117a65")
    ax2.tick_params(axis="y", colors="#117a65")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], loc="upper left", frameon=False, fontsize=8)

    ax = axes[2]
    ax.plot(dts, pos_hmm, color="#2c3e50", lw=1.0, label="HMM exposure")
    ax.plot(dts, pos_comb, color="#e67e22", lw=1.5, label="HMM x gate")
    ax.fill_between(dts, 0, 1, where=(gate < 1), color="#e67e22", alpha=0.12, step="mid")
    ax.set_ylim(-0.03, 1.05)
    ax.set_ylabel("exposure")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    train_end = meta["data"]["train_end"]
    if train_end < len(dts):
        for a in axes:
            a.axvline(dts[train_end], color="k", ls="--", lw=0.8, alpha=0.6)
        axes[0].text(dts[train_end], axes[0].get_ylim()[1] * 0.92, "  OOS start",
                     fontsize=8, va="top")
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.tight_layout()
    fig.savefig(FIG / "regimes.png", dpi=160)
    plt.close(fig)


def fig_equity():
    eq = read_csv(ROOT / "out/regime/equity.csv")
    meta = json.loads((ROOT / "out/regime/metrics.json").read_text())
    dts = dates(eq)
    fig, ax = plt.subplots(figsize=(10, 4.6))
    styles = [("buy_hold", "buy & hold", "#7f8c8d"), ("hmm", "HMM", "#2c3e50"),
              ("hmm_perc", "HMM + percolation gate", "#e67e22")]
    for key, label, color in styles:
        b = meta["backtest"][key]
        ax.plot(dts, col(eq, key), label=f"{label}  (Sharpe {b['sharpe']:.2f}, MDD {b['max_drawdown']*100:.1f}%)",
                color=color, lw=1.5 if key != "hmm_perc" else 2.0)
    ax.set_ylabel("equity (start = 1)")
    ax.set_title("Out-of-sample equity, 5 bps cost per unit turnover")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "equity.png", dpi=160)
    plt.close(fig)


def fig_percolation():
    scan = read_csv(ROOT / "out/perc_scan_bond.csv")
    p = col(scan, "p")
    s = col(scan, "s_max")
    span = col(scan, "spanning")
    chi = col(scan, "susceptibility")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.plot(p, s, "o-", ms=3, color="#1f4e79", label="S = largest / L^2")
    ax.plot(p, span, "s-", ms=3, color="#b03a2e", label="spanning probability")
    ax.plot(p, chi / max(np.nanmax(chi), 1e-9), "^-", ms=3, color="#117a65",
            label="chi / max(chi)")
    ax.axvline(0.5, color="k", ls="--", lw=0.8, label="p_c = 1/2")
    ax.set_xlabel("bond occupation p")
    ax.set_ylabel("value")
    ax.set_title("2D bond percolation, open L=128, 200 reps")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for kind, color in [("bond", "#1f4e79"), ("site", "#b03a2e")]:
        rows = read_csv(ROOT / f"out/perc_scaling_{kind}.csv")
        L = col(rows, "L")
        pc = col(rows, "p_cross")
        ax.plot(L, pc, "o", color=color, label=f"{kind} p_c(L)")
        meta = json.loads((ROOT / f"out/perc_scaling_{kind}.json").read_text())
        xs = np.linspace(min(L), max(L), 100)
        ax.plot(xs, meta["pc_infinite"] + meta["fit_slope"] * xs ** (-0.75), "-",
                color=color, lw=1.2, alpha=0.7,
                label=f"{kind} fit p_c(inf)={meta['pc_infinite']:.4f}")
        ax.axhline(0.5 if kind == "bond" else 0.592746, color=color, ls=":", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("L")
    ax.set_ylabel("crossing p")
    ax.set_title("finite-size scaling, nu = 4/3")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "percolation.png", dpi=160)
    plt.close(fig)


def fig_corr():
    rows = read_csv(ROOT / "out/connectedness/corr_last.csv")
    symbols = list(rows[0].keys())
    mat = np.array([[float(r[s]) for s in symbols] for r in rows])
    meta = json.loads((ROOT / "out/connectedness/summary.json").read_text())
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(symbols)))
    ax.set_yticks(range(len(symbols)))
    ax.set_xticklabels(symbols, rotation=45, fontsize=8)
    ax.set_yticklabels(symbols, fontsize=8)
    ax.set_title(f"last-window correlation, theta_c={meta['theta_c_last']:.3f}")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(FIG / "correlation.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    fig_regimes()
    fig_equity()
    fig_percolation()
    fig_corr()
    print("wrote figures to", FIG)
