# Changelog

## 0.1.0 — 2026-11-14

First public snapshot.

- C11 core: RNG, statistics, CSV, discrete HMM, diagonal GMM-HMM.
- Percolation: union-find, 2D site/bond lattice, finite-size crossing,
  ER/configuration model, correlation-threshold graphs.
- C++17 application: walk-forward regime model, percolation gate, backtest,
  JSON/Markdown reports.
- CLIs: `hp_synth`, `hp_regime`, `hp_connectedness`, `hp_perc_lab`.
- Python (ctypes) and Julia (ccall) bindings.
- C/C++/Python/Julia tests; ASan-checked C and C++ paths.
- Mathematical notes separating proofs, citations, and numerical observations.
- Configuration-model solver computes the neighbour extinction probability
  before the giant fraction, fixing the non-Poisson case.
- Out-of-sample accounting starts at the first return after `train_end`.
