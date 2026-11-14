# Changelog

## 0.1.3 — 2026-11-14

- core: add `hp_version()` / `HP_VERSION_*` so bindings do not duplicate the
  library version.
- bindings: Python `hpmarket.__version__` and Julia `HP.lib_version()` now read
  the runtime version from C.
- docs: add `docs/extending.md` with the extension seams and compatibility
  rules.

## 0.1.2 — 2026-11-14

- docs: add DOI / publisher / repository links for every cited work, plus a
  machine-readable `docs/references.bib`.
- docs: fix the Dempster–Laird–Rubin page range (1–22) and the
  Erdős–Rényi page range (17–61).
- docs: render the formulas in README, design.md and results.md with the same
  GitHub-compatible delimiters as mathematics.md.

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
