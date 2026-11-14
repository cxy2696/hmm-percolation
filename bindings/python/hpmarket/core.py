"""Thin, RAII-style Python objects over the C library.

All numerical work happens in C; these wrappers only marshal arrays.  Lists
are returned instead of numpy arrays so the package has no third-party
dependency (numpy inputs are accepted and flattened).
"""
from __future__ import annotations

import ctypes as C
import math
import os
from typing import Iterable, Sequence

from . import _capi as capi
from ._capi import HPError, check, lib

INIT_QUANTILE = 0
INIT_RANDOM = 1
SITE = 0
BOND = 1

__all__ = [
    "HPError", "RNG", "DiscreteHMM", "GaussianHMM", "Panel", "Lattice",
    "INIT_QUANTILE", "INIT_RANDOM", "SITE", "BOND",
    "corr_matrix", "corr_critical", "corr_sweep", "graph_giant",
    "perc_crossing", "perc_spanning_prob", "er_giant", "config_giant",
    "bond_threshold", "backtest", "select_ghmm",
    "vol_target_signal", "momentum_signal", "state_exposure_signal",
    "multiply_signals", "constant_signal",
]


def _flat(x) -> list[float]:
    if hasattr(x, "ravel"):
        return [float(v) for v in x.ravel()]
    out: list[float] = []
    for row in x:
        if isinstance(row, (list, tuple)):
            out.extend(float(v) for v in row)
        else:
            out.append(float(row))
    return out


def _matrix(x) -> tuple[list[float], int, int]:
    if hasattr(x, "shape") and len(x.shape) == 2:
        T, d = int(x.shape[0]), int(x.shape[1])
        return _flat(x), T, d
    rows = list(x)
    T = len(rows)
    if T == 0:
        return [], 0, 0
    d = len(rows[0])
    out: list[float] = []
    for row in rows:
        if len(row) != d:
            raise ValueError("ragged input matrix")
        out.extend(float(v) for v in row)
    return out, T, d


def _path(p) -> bytes:
    return os.fspath(p).encode()


def _dvec(values: Iterable[float]):
    vals = [float(v) for v in values]
    return (C.c_double * len(vals))(*vals)


def _u32vec(values: Iterable[int]):
    vals = [int(v) for v in values]
    return (C.c_uint32 * len(vals))(*vals)


def _read(ptr, n: int) -> list[float]:
    return [float(ptr[i]) for i in range(n)]


def _take(ptr, base: int, n: int) -> list[float]:
    return [float(ptr[base + i]) for i in range(n)]


def _cstr(p) -> str:
    return C.cast(p, C.c_char_p).value.decode("utf-8")


class RNG:
    """xoshiro256++ stream; same seed/stream => same numbers as C."""

    def __init__(self, seed: int = 42, stream: int = 0):
        self._rng = capi.RNG()
        lib.hp_rng_seed_stream(C.byref(self._rng), C.c_uint64(seed), C.c_uint64(stream))

    @property
    def raw(self):
        return self._rng

    def uniform(self) -> float:
        return float(lib.hp_rng_uniform(C.byref(self._rng)))

    def normal(self) -> float:
        return float(lib.hp_rng_normal(C.byref(self._rng)))


class DiscreteHMM:
    def __init__(self, n_states: int, n_symbols: int):
        self._ptr = lib.hp_hmm_new(n_states, n_symbols)
        if not self._ptr:
            raise ValueError("bad HMM dimensions")

    @classmethod
    def _from_ptr(cls, ptr):
        obj = cls.__new__(cls)
        obj._ptr = ptr
        return obj

    @classmethod
    def load(cls, path: str) -> "DiscreteHMM":
        ptr = lib.hp_hmm_load(_path(path))
        if not ptr:
            raise HPError(f"cannot load discrete HMM from {path}")
        return cls._from_ptr(ptr)

    def close(self) -> None:
        if getattr(self, "_ptr", None):
            lib.hp_hmm_free(self._ptr)
            self._ptr = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def n_states(self) -> int:
        return int(self._ptr.contents.n_states)

    @property
    def n_symbols(self) -> int:
        return int(self._ptr.contents.n_symbols)

    def params(self) -> dict:
        n, m = self.n_states, self.n_symbols
        c = self._ptr.contents
        return {
            "pi": _read(c.pi, n),
            "a": [[float(c.a[i * n + j]) for j in range(n)] for i in range(n)],
            "b": [[float(c.b[i * m + k]) for k in range(m)] for i in range(n)],
        }

    def validate(self) -> None:
        check(lib.hp_hmm_validate(self._ptr))

    def train(self, obs: Sequence[int], *, max_iter: int = 200, tol: float = 1e-6,
              restarts: int = 4, seed: int = 42, alpha: float = 1.0,
              init: int = INIT_QUANTILE, verbose: int = 0) -> dict:
        T = len(obs)
        o = _u32vec(obs)
        opts = capi.HMMTrainOpts(max_iter, tol, restarts, seed, alpha, init, verbose)
        rep = capi.HMMReport()
        check(lib.hp_hmm_train(o, T, self.n_states, self.n_symbols, C.byref(opts),
                               self._ptr, C.byref(rep)))
        return {
            "loglik": float(rep.loglik), "bic": float(rep.bic), "aic": float(rep.aic),
            "n_params": int(rep.n_params), "iterations": int(rep.iterations),
            "converged": bool(rep.converged), "restarts": int(rep.restarts),
        }

    def forward(self, obs: Sequence[int]) -> tuple[list[float], float]:
        T, n = len(obs), self.n_states
        o = _u32vec(obs)
        la = (C.c_double * (T * n))()
        ll = C.c_double()
        check(lib.hp_hmm_forward_log(self._ptr, o, T, la, C.byref(ll)))
        return _read(la, T * n), float(ll.value)

    def viterbi(self, obs: Sequence[int]) -> tuple[list[int], float]:
        T, n = len(obs), self.n_states
        o = _u32vec(obs)
        path = (C.c_uint32 * T)()
        lp = C.c_double()
        check(lib.hp_hmm_viterbi(self._ptr, o, T, path, C.byref(lp)))
        return [int(path[t]) for t in range(T)], float(lp.value)

    def posterior(self, obs: Sequence[int]) -> tuple[list[list[float]], float]:
        T, n = len(obs), self.n_states
        o = _u32vec(obs)
        g = (C.c_double * (T * n))()
        ll = C.c_double()
        check(lib.hp_hmm_posterior(self._ptr, o, T, g, C.byref(ll)))
        return [_read(g, T * n)[t * n:(t + 1) * n] for t in range(T)], float(ll.value)

    def filtered(self, obs: Sequence[int]) -> tuple[list[list[float]], float]:
        T, n = len(obs), self.n_states
        o = _u32vec(obs)
        f = (C.c_double * (T * n))()
        ll = C.c_double()
        check(lib.hp_hmm_filter(self._ptr, o, T, f, C.byref(ll)))
        return [_read(f, T * n)[t * n:(t + 1) * n] for t in range(T)], float(ll.value)

    def save(self, path: str) -> None:
        check(lib.hp_hmm_save(self._ptr, _path(path)))


class GaussianHMM:
    def __init__(self, n_states: int, dim: int, n_mix: int = 1):
        self._ptr = lib.hp_ghmm_new(n_states, dim, n_mix)
        if not self._ptr:
            raise ValueError("bad GHMM dimensions")

    @classmethod
    def _from_ptr(cls, ptr):
        obj = cls.__new__(cls)
        obj._ptr = ptr
        return obj

    @classmethod
    def load(cls, path: str) -> "GaussianHMM":
        ptr = lib.hp_ghmm_load(_path(path))
        if not ptr:
            raise HPError(f"cannot load GHMM from {path}")
        return cls._from_ptr(ptr)

    def close(self) -> None:
        if getattr(self, "_ptr", None):
            lib.hp_ghmm_free(self._ptr)
            self._ptr = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def n_states(self) -> int:
        return int(self._ptr.contents.n_states)

    @property
    def dim(self) -> int:
        return int(self._ptr.contents.dim)

    @property
    def n_mix(self) -> int:
        return int(self._ptr.contents.n_mix)

    def params(self) -> dict:
        n, d, k = self.n_states, self.dim, self.n_mix
        c = self._ptr.contents
        return {
            "pi": _read(c.pi, n),
            "a": [[float(c.a[i * n + j]) for j in range(n)] for i in range(n)],
            "mix": [[float(c.mix[i * k + m]) for m in range(k)] for i in range(n)],
            "mean": [[_take(c.mean, (i * k + m) * d, d) for m in range(k)] for i in range(n)],
            "var": [[_take(c.var, (i * k + m) * d, d) for m in range(k)] for i in range(n)],
        }

    def state_mean(self, state: int, dim: int = 0) -> float:
        p = self.params()
        return sum(p["mix"][state][m] * p["mean"][state][m][dim] for m in range(self.n_mix))

    def validate(self) -> None:
        check(lib.hp_ghmm_validate(self._ptr))

    def train(self, x, *, max_iter: int = 200, tol: float = 1e-6, restarts: int = 4,
              seed: int = 42, alpha: float = 1.0, var_floor: float = 1e-4,
              mix_floor: float = 1e-8, init: int = INIT_QUANTILE,
              verbose: int = 0) -> dict:
        flat, T, d = _matrix(x)
        if d != self.dim:
            raise ValueError(f"expected dim {self.dim}, got {d}")
        xa = _dvec(flat)
        opts = capi.GHMMTrainOpts(max_iter, tol, restarts, seed, alpha, var_floor,
                                  mix_floor, init, verbose)
        rep = capi.GHMMReport()
        check(lib.hp_ghmm_train(xa, T, d, self.n_states, self.n_mix, C.byref(opts),
                                self._ptr, C.byref(rep)))
        return {
            "loglik": float(rep.loglik), "bic": float(rep.bic), "aic": float(rep.aic),
            "n_params": int(rep.n_params), "iterations": int(rep.iterations),
            "converged": bool(rep.converged), "restarts": int(rep.restarts),
        }

    def forward(self, x) -> tuple[list[float], float]:
        flat, T, d = _matrix(x)
        if d != self.dim:
            raise ValueError(f"expected dim {self.dim}, got {d}")
        la = (C.c_double * (T * self.n_states))()
        ll = C.c_double()
        check(lib.hp_ghmm_forward_log(self._ptr, _dvec(flat), T, la, C.byref(ll)))
        return _read(la, T * self.n_states), float(ll.value)

    def posterior(self, x) -> tuple[list[list[float]], float]:
        flat, T, d = _matrix(x)
        n = self.n_states
        g = (C.c_double * (T * n))()
        ll = C.c_double()
        check(lib.hp_ghmm_posterior(self._ptr, _dvec(flat), T, g, C.byref(ll)))
        vals = _read(g, T * n)
        return [vals[t * n:(t + 1) * n] for t in range(T)], float(ll.value)

    def filtered(self, x) -> tuple[list[list[float]], float]:
        flat, T, d = _matrix(x)
        n = self.n_states
        f = (C.c_double * (T * n))()
        ll = C.c_double()
        check(lib.hp_ghmm_filter(self._ptr, _dvec(flat), T, f, C.byref(ll)))
        vals = _read(f, T * n)
        return [vals[t * n:(t + 1) * n] for t in range(T)], float(ll.value)

    def viterbi(self, x) -> tuple[list[int], float]:
        flat, T, d = _matrix(x)
        path = (C.c_uint32 * T)()
        lp = C.c_double()
        check(lib.hp_ghmm_viterbi(self._ptr, _dvec(flat), T, path, C.byref(lp)))
        return [int(path[t]) for t in range(T)], float(lp.value)

    def save(self, path: str) -> None:
        check(lib.hp_ghmm_save(self._ptr, _path(path)))


def select_ghmm(x, n_min: int, n_max: int, *, n_mix: int = 1, **train_kwargs):
    """BIC model selection; returns (best_model, best_n, table)."""
    flat, T, d = _matrix(x)
    model = GaussianHMM(n_min, d, n_mix)
    opts = capi.GHMMTrainOpts(
        train_kwargs.get("max_iter", 150), train_kwargs.get("tol", 1e-6),
        train_kwargs.get("restarts", 4), train_kwargs.get("seed", 42),
        train_kwargs.get("alpha", 1.0), train_kwargs.get("var_floor", 1e-4),
        train_kwargs.get("mix_floor", 1e-8), train_kwargs.get("init", INIT_QUANTILE), 0)
    table = (capi.GHMMModelInfo * (n_max - n_min + 1))()
    best_n = C.c_size_t()
    check(lib.hp_ghmm_select(_dvec(flat), T, d, n_min, n_max, n_mix, C.byref(opts),
                             model._ptr, C.byref(best_n), table))
    rows = [{
        "n_states": int(table[i].n_states), "loglik": float(table[i].loglik),
        "bic": float(table[i].bic), "aic": float(table[i].aic),
        "n_params": int(table[i].n_params), "converged": bool(table[i].converged),
    } for i in range(n_max - n_min + 1)]
    return model, int(best_n.value), rows


class Panel:
    def __init__(self, ptr=None):
        self._ptr = ptr

    @classmethod
    def load(cls, path: str) -> "Panel":
        st = C.c_int()
        ptr = lib.hp_panel_load_csv(_path(path), C.byref(st))
        if not ptr:
            raise HPError(f"cannot load panel {path}: status {st.value}")
        return cls(ptr)

    def close(self):
        if getattr(self, "_ptr", None):
            lib.hp_panel_free(self._ptr)
            self._ptr = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def T(self) -> int:
        return int(self._ptr.contents.T)

    @property
    def N(self) -> int:
        return int(self._ptr.contents.N)

    @property
    def dates(self) -> list[str]:
        return [_cstr(self._ptr.contents.dates[t]) for t in range(self.T)]

    @property
    def symbols(self) -> list[str]:
        return [_cstr(self._ptr.contents.symbols[j]) for j in range(self.N)]

    @property
    def close_prices(self) -> list[float]:
        return _read(self._ptr.contents.close, self.T * self.N)

    @property
    def returns(self) -> list[float]:
        return _read(self._ptr.contents.ret, self.T * self.N)

    def ew_returns(self) -> list[float]:
        out = (C.c_double * self.T)()
        check(lib.hp_panel_ew_returns(self._ptr, out))
        return _read(out, self.T)

    def return_series(self, j: int) -> list[float]:
        out = (C.c_double * self.T)()
        check(lib.hp_panel_return_series(self._ptr, j, out))
        return _read(out, self.T)

    def symbol_index(self, name: str) -> int:
        low = [s.lower() for s in self.symbols]
        try:
            return low.index(name.lower())
        except ValueError:
            return -1

    def save(self, path: str) -> None:
        check(lib.hp_panel_save_csv(self._ptr, _path(path)))


class Lattice:
    def __init__(self, L: int, kind: str = "bond", periodic: bool = False):
        k = BOND if kind == "bond" else SITE
        self._ptr = lib.hp_lattice_new(L, k, 1 if periodic else 0)
        if not self._ptr:
            raise ValueError("bad lattice")
        self.L = L
        self.kind = kind
        self.periodic = periodic

    def close(self):
        if getattr(self, "_ptr", None):
            lib.hp_lattice_free(self._ptr)
            self._ptr = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def bernoulli(self, p: float, rng: RNG) -> None:
        check(lib.hp_lattice_bernoulli(self._ptr, p, C.byref(rng.raw)))

    def analyze(self) -> dict:
        st = capi.LatticeStats()
        check(lib.hp_lattice_analyze(self._ptr, C.byref(st)))
        return {
            "n_sites": int(st.n_sites), "n_occupied": int(st.n_occupied),
            "n_clusters": int(st.n_clusters), "max_cluster": int(st.max_cluster),
            "s_max": float(st.s_max), "susceptibility": float(st.susceptibility),
            "occupancy": float(st.occupancy), "spanning": bool(st.spanning),
            "spanning_x": bool(st.spanning_x), "spanning_y": bool(st.spanning_y),
        }

    def cluster_sizes(self) -> list[int]:
        n = C.c_size_t()
        ptr = lib.hp_lattice_cluster_sizes(self._ptr, C.byref(n))
        if not ptr:
            return []
        try:
            return [int(ptr[i]) for i in range(n.value)]
        finally:
            lib.hp_free(C.cast(ptr, C.c_void_p))


def perc_crossing(L: int, *, kind: str = "bond", periodic: bool = False,
                  reps: int = 100, seed: int = 1, p_lo: float = 0.35,
                  p_hi: float = 0.65, iters: int = 14) -> tuple[float, float]:
    exp = capi.PercExp(L, BOND if kind == "bond" else SITE, 1 if periodic else 0,
                       reps, seed)
    pc, sp = C.c_double(), C.c_double()
    check(lib.hp_perc_crossing(C.byref(exp), p_lo, p_hi, iters, C.byref(pc), C.byref(sp)))
    return float(pc.value), float(sp.value)


def perc_spanning_prob(L: int, p: float, *, kind: str = "bond",
                       periodic: bool = False, reps: int = 100,
                       seed: int = 1) -> float:
    exp = capi.PercExp(L, BOND if kind == "bond" else SITE, 1 if periodic else 0,
                       reps, seed)
    return float(lib.hp_perc_spanning_prob(C.byref(exp), p))


def er_giant(mean_degree: float) -> float:
    return float(lib.hp_perc_er_giant_theory(mean_degree))


def bond_threshold(mean_degree: float) -> float:
    return float(lib.hp_perc_bond_threshold(mean_degree))


def config_giant(pk: Sequence[float], transmissibility: float = 1.0) -> dict:
    arr = _dvec(pk)
    g, r0, mean = C.c_double(), C.c_double(), C.c_double()
    check(lib.hp_perc_config_giant(arr, len(pk) - 1, transmissibility,
                                   C.byref(g), C.byref(r0), C.byref(mean)))
    return {"giant": float(g.value), "r0": float(r0.value), "mean_degree": float(mean.value)}


def corr_matrix(x) -> list[list[float]]:
    flat, T, N = _matrix(x)
    out = (C.c_double * (N * N))()
    lib.hp_corr_matrix(_dvec(flat), T, N, out)
    vals = _read(out, N * N)
    return [vals[i * N:(i + 1) * N] for i in range(N)]


def corr_critical(corr, theta_min: float = 0.1, theta_max: float = 0.99,
                  steps: int = 45) -> dict:
    rows = list(corr)
    n = len(rows)
    flat = [float(v) for row in rows for v in row]
    crit = capi.CorrCritical()
    check(lib.hp_corr_critical(_dvec(flat), n, theta_min, theta_max, steps, C.byref(crit)))
    return {
        "theta_peak": float(crit.theta_c),
        "theta_half": float(crit.theta_half),
        "chi_max": float(crit.chi_max),
        "giant_at_peak": float(crit.giant_at_c),
        "censored_high": bool(crit.censored_high),
        "censored_low": bool(crit.censored_low),
    }


def corr_sweep(corr, thetas: Sequence[float]) -> list[dict]:
    rows = list(corr)
    n = len(rows)
    flat = [float(v) for row in rows for v in row]
    th = _dvec(thetas)
    out = (capi.PercPoint * len(thetas))()
    check(lib.hp_corr_sweep(_dvec(flat), n, th, len(thetas), out))
    return [{
        "theta": float(out[i].theta), "n_edges": int(out[i].n_edges),
        "n_components": int(out[i].n_components),
        "giant_fraction": float(out[i].giant_fraction),
        "susceptibility": float(out[i].susceptibility),
        "mean_degree": float(out[i].mean_degree),
    } for i in range(len(thetas))]


def graph_giant(corr, theta: float) -> dict:
    rows = list(corr)
    n = len(rows)
    flat = [float(v) for row in rows for v in row]
    g = lib.hp_graph_from_corr(_dvec(flat), n, theta)
    if not g:
        raise HPError("graph_from_corr failed")
    comp = capi.Components()
    try:
        check(lib.hp_graph_components(g, -math.inf, C.byref(comp)))
        return {
            "n_edges": int(lib.hp_graph_edge_count(g)),
            "n_components": int(comp.n_components),
            "giant": int(comp.giant),
            "giant_fraction": float(comp.giant_fraction),
            "susceptibility": float(comp.susceptibility),
            "sizes": [int(comp.sizes[i]) for i in range(comp.n_components)],
        }
    finally:
        lib.hp_components_free(C.byref(comp))
        lib.hp_graph_free(g)


def vol_target_signal(asset_ret: Sequence[float], *, lookback: int = 20,
                      target_vol: float = 0.10, min_position: float = 0.0,
                      max_position: float = 1.0,
                      periods_per_year: float = 252.0) -> list[float]:
    """Causal inverse-vol position.  NaN before enough history."""
    ret = _dvec(asset_ret)
    T = len(asset_ret)
    out = (C.c_double * T)()
    opts = capi.VolTargetOpts()
    lib.hp_vol_target_opts_default(C.byref(opts))
    opts.lookback = int(lookback)
    opts.target_vol = float(target_vol)
    opts.min_position = float(min_position)
    opts.max_position = float(max_position)
    opts.periods_per_year = float(periods_per_year)
    check(lib.hp_signal_vol_target(ret, T, C.byref(opts), out))
    return _read(out, T)


def momentum_signal(asset_ret: Sequence[float], *, lookback: int = 60,
                    threshold: float = 0.0, long_position: float = 1.0,
                    flat_position: float = 0.0) -> list[float]:
    """Long if trailing cumulative log return > threshold, else flat."""
    ret = _dvec(asset_ret)
    T = len(asset_ret)
    out = (C.c_double * T)()
    opts = capi.MomentumOpts()
    lib.hp_momentum_opts_default(C.byref(opts))
    opts.lookback = int(lookback)
    opts.threshold = float(threshold)
    opts.long_position = float(long_position)
    opts.flat_position = float(flat_position)
    check(lib.hp_signal_momentum(ret, T, C.byref(opts), out))
    return _read(out, T)


def state_exposure_signal(filtered: Sequence[Sequence[float]], n_states: int,
                          exposure: Sequence[float]) -> list[float]:
    flat, T, d = _matrix(filtered)
    if d != int(n_states) or len(exposure) != int(n_states):
        raise ValueError("state dimension mismatch")
    filt = _dvec(flat)
    expo = _dvec(exposure)
    out = (C.c_double * T)()
    check(lib.hp_signal_state_exposure(filt, T, int(n_states), expo, out))
    return _read(out, T)


def multiply_signals(a: Sequence[float], b: Sequence[float]) -> list[float]:
    if len(a) != len(b):
        raise ValueError("signal length mismatch")
    T = len(a)
    av, bv = _dvec(a), _dvec(b)
    out = (C.c_double * T)()
    check(lib.hp_signal_multiply(av, bv, T, out))
    return _read(out, T)


def constant_signal(n: int, position: float = 1.0) -> list[float]:
    if n <= 0:
        raise ValueError("n must be positive")
    out = (C.c_double * n)()
    check(lib.hp_signal_constant(n, float(position), out))
    return _read(out, n)


def backtest(asset_ret: Sequence[float], signal: Sequence[float], *,
             cost_bps: float = 5.0, apply_costs: bool = True,
             periods_per_year: float = 252.0) -> dict:
    T = len(asset_ret)
    if len(signal) != T:
        raise ValueError("signal length mismatch")
    ret = _dvec(asset_ret)
    sig = _dvec(signal)
    opts = capi.BTOpts(cost_bps, 1 if apply_costs else 0, periods_per_year)
    res, m = capi.BTResult(), capi.BTMetrics()
    check(lib.hp_backtest(ret, T, sig, C.byref(opts), C.byref(res), C.byref(m)))
    try:
        return {
            "total_return": float(m.total_return), "cagr": float(m.cagr),
            "ann_vol": float(m.ann_vol), "sharpe": float(m.sharpe),
            "max_drawdown": float(m.max_drawdown), "calmar": float(m.calmar),
            "avg_turnover": float(m.avg_turnover), "exposure": float(m.exposure),
            "hit_rate": float(m.hit_rate), "final_equity": float(m.final_equity),
            "n_days": int(m.n_days),
            "equity": _read(res.equity, T),
            "strategy_ret": _read(res.strategy_ret, T),
            "position": _read(res.position, T),
        }
    finally:
        lib.hp_bt_result_free(C.byref(res))
