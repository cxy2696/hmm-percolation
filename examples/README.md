# Examples

These are deliberately small and do not use financial data.

```bash
make examples
./build/ex_discrete_hmm
./build/ex_gaussian_hmm
./build/ex_percolation
```

- `discrete_hmm.c` — two-state, three-symbol HMM; train, Viterbi, accuracy.
- `gaussian_hmm.c` — one-dimensional two-regime GMM-HMM with a variance floor.
- `percolation.c` — 2D bond lattice, finite-size crossing, ER and configuration model.

The finance-specific pipeline lives in `apps/hp_regime.cpp`; the library itself
does not know about returns or assets.
