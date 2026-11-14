#!/usr/bin/env bash
# Full deterministic experiment set.  Usage:
#   ./scripts/reproduce.sh
#   JOBS=8 ./scripts/reproduce.sh
set -euo pipefail
cd "$(dirname "$0")/.."

JOBS="${JOBS:-4}"
PYTHON="${PYTHON:-python3}"
JULIA="${JULIA:-julia}"

echo "== build =="
make -j"$JOBS"
make shared
make check

echo "== synthetic market =="
mkdir -p data out
./build/hp_synth --days 3000 --assets 10 --regimes 3 --seed 20240101 \
  --out data/synthetic_panel.csv --truth data/synthetic_truth.csv

echo "== connectedness =="
./build/hp_connectedness --csv data/synthetic_panel.csv --window 60 \
  --outdir out/connectedness

echo "== regime variants =="
./build/hp_regime --csv data/synthetic_panel.csv --states 3 --outdir out/regime
./build/hp_regime --csv data/synthetic_panel.csv --states 2 --outdir out/regime_s2
./build/hp_regime --csv data/synthetic_panel.csv --states 0 --outdir out/regime_s0
./build/hp_regime --csv data/synthetic_panel.csv --states 3 --no-perc \
  --outdir out/regime_noperc
./build/hp_regime --csv data/synthetic_panel.csv --states 3 --conn-window 30 \
  --outdir out/regime_w30

echo "== percolation lab =="
./build/hp_perc_lab --mode scan --kind bond --L 128 --pmin 0.35 --pmax 0.65 \
  --steps 31 --reps 200 --out out/perc_scan_bond.csv
./build/hp_perc_lab --mode scaling --kind bond --sizes 32,64,128,256 \
  --reps 200 --iters 16 --out out/perc_scaling_bond.csv \
  --json out/perc_scaling_bond.json
./build/hp_perc_lab --mode scaling --kind site --sizes 32,64,128,256 \
  --reps 200 --iters 16 --out out/perc_scaling_site.csv \
  --json out/perc_scaling_site.json
./build/hp_perc_lab --mode meanfield --cmin 0.5 --cmax 3 --steps 26 \
  --out out/perc_meanfield.csv

echo "== Python bindings =="
PYTHONPATH=bindings/python "$PYTHON" -m unittest discover -s bindings/python/tests -v

echo "== Julia bindings =="
if "$JULIA" --version >/dev/null 2>&1; then
  "$JULIA" --startup-file=no bindings/julia/test/runtests.jl
else
  echo "julia not runnable; skipping (set JULIA=... to override)"
fi

echo "== figures =="
if [ -x build/venv/bin/python ]; then
  MPLCONFIGDIR=build/mplcache build/venv/bin/python scripts/plot_results.py
else
  echo "build/venv missing; skipping figures (see docs/05 for setup)"
fi

echo
echo "== summary =="
"$PYTHON" - <<'PY'
import json
from pathlib import Path
for d in ["out/regime", "out/regime_s2", "out/regime_s0"]:
    p = Path(d) / "metrics.json"
    if not p.exists():
        continue
    m = json.loads(p.read_text())
    print(f"{d}: states={m['model']['n_states']} "
          f"oos={m['backtest']['oos_days']} "
          f"hmm Sharpe={m['backtest']['hmm']['sharpe']:.2f} "
          f"gate Sharpe={m['backtest']['hmm_perc']['sharpe']:.2f} "
          f"MDD={m['backtest']['hmm_perc']['max_drawdown']*100:.1f}%")
for k in ["bond", "site"]:
    p = Path(f"out/perc_scaling_{k}.json")
    if p.exists():
        m = json.loads(p.read_text())
        print(f"{k} p_c(inf)={m['pc_infinite']:.4f} known={m['pc_known']:.4f} R2={m['fit_r2']:.3f}")
PY
