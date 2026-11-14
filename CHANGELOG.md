# Changelog

## 0.1.1 — 2026-11-14

- docs: switch to GitHub-compatible `$...$` / `$$...$$` math delimiters and add
  the blank lines display math needs.
- docs: use `\lbrace`, `\lt`, `\gt` where Markdown escaping would corrupt the
  TeX, and clarify the EM positivity assumption.
- docs: spell out that `theta_half` is a finite-sample convention and that the
  backtest signal, not the realized return, is the causal object.

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
