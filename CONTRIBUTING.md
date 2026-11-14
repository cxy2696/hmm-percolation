# Contributing

This is a small research library, not a framework. Changes are welcome if they
keep the core dependency-free and the claims checkable.

## Build and test

```bash
make clean && make -j
make check
make shared
PYTHONPATH=bindings/python python3 -m unittest discover -s bindings/python/tests -v
julia --startup-file=no bindings/julia/test/runtests.jl
```

If a change touches numerical paths, also run the ASan build from
`docs/reproduction.md` and `./scripts/reproduce.sh`.

## Style

- C11, no compiler extensions. Keep `-Wall -Wextra -Wpedantic` clean.
- No external numerical libraries in `src/core`, `src/hmm`, `src/perc`.
- Public C structs are append-only. Changing a field order or a function
  signature is an ABI break: bump the model/format version and document it.
- Comments should explain a non-obvious decision, not restate the code.
- Chinese or English are both fine in docs; keep one language per file.

## Claims and proofs

Every nontrivial claim in `docs/` must be in one of three buckets:

1. theorem with proof;
2. citation to a specific result;
3. explicit numerical observation.

Do not move an item from bucket 3 to bucket 1 without a proof. In particular,
finite-size estimates and backtest numbers are not theorems.

## Data and reproducibility

- All random experiments take an explicit seed.
- Do not commit generated caches; `data/` and `out/` are small, fixed-seed
  artifacts and may be regenerated.
- If you change a default that affects `out/`, regenerate the artifacts and
  update `docs/results.md` and the checkpoints in `docs/reproduction.md`.
- Never commit personal paths, tokens, or machine-specific output.

