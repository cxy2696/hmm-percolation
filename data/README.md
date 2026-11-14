# Synthetic data

`hp_synth` generates this panel. It is not market data.

Default command:

```bash
./build/hp_synth --days 3000 --assets 10 --regimes 3 --seed 20240101 \
  --out data/synthetic_panel.csv --truth data/synthetic_truth.csv
```

Files:

- `synthetic_panel.csv`: `date,A1,...,A10`, close prices. 3000 weekdays,
  2012-01-02 through 2023-06-30. No holidays.
- `synthetic_truth.csv`: `date,regime,market_ret,A1,...,A10`; regime 0/1/2
  and log returns.

The generator has three sticky regimes with regime-dependent drift,
volatility, and idiosyncratic risk. In the bear state idiosyncratic risk falls,
so cross-sectional correlations rise. See `src/market/synth.c`.

Regenerating the files with another seed invalidates `out/`; rerun
`scripts/reproduce.sh` afterwards.
