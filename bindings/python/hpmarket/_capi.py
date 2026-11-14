"""ctypes declarations for the C core ``libhp``.

The C library is the single source of truth: this module only mirrors the
public structs and function signatures.  Build the shared library with
``make shared`` first, or point ``HPMARKET_LIB`` at it.
"""
from __future__ import annotations

import ctypes as C
import os
import sys
from pathlib import Path

__all__ = [
    "HPError", "lib", "check",
    "RNG", "HMM", "GHMM", "Panel",
    "HMMTrainOpts", "HMMReport", "GHMMTrainOpts", "GHMMReport", "GHMMModelInfo",
    "BTOpts", "BTMetrics", "BTResult",
    "Lattice", "LatticeStats", "PercExp", "CorrCritical", "PercPoint",
    "Components",
]

c_double_p = C.POINTER(C.c_double)
c_u32_p = C.POINTER(C.c_uint32)
c_size_p = C.POINTER(C.c_size_t)


class HPError(RuntimeError):
    """Non-zero hp_status returned by the C library."""


def _candidate_paths():
    root = Path(__file__).resolve().parents[3]
    env = os.environ.get("HPMARKET_LIB")
    if env:
        yield Path(env)
    names = ["libhp.dylib", "libhp.so", "hp.dll"]
    for name in names:
        yield root / "build" / name
    for name in names:
        yield Path(name)


def _load_library():
    last = None
    for path in _candidate_paths():
        try:
            if path.exists():
                return C.CDLL(str(path))
        except OSError as exc:  # pragma: no cover - platform detail
            last = exc
    hint = "build it with `make shared`; set HPMARKET_LIB to the dylib/so path"
    raise ImportError(f"cannot load libhp ({hint})" + (f": {last}" if last else ""))


lib = _load_library()


def check(status: int) -> None:
    if status != 0:
        msg = lib.hp_status_message(status)
        text = msg.decode("utf-8") if msg else "unknown hp_status"
        raise HPError(f"{text} (status {status})")


# ---------------------------------------------------------------- structs
class RNG(C.Structure):
    _fields_ = [("s", C.c_uint64 * 4)]


class HMM(C.Structure):
    _fields_ = [
        ("n_states", C.c_size_t),
        ("n_symbols", C.c_size_t),
        ("pi", c_double_p),
        ("a", c_double_p),
        ("b", c_double_p),
    ]


class GHMM(C.Structure):
    _fields_ = [
        ("n_states", C.c_size_t),
        ("dim", C.c_size_t),
        ("n_mix", C.c_size_t),
        ("pi", c_double_p),
        ("a", c_double_p),
        ("mix", c_double_p),
        ("mean", c_double_p),
        ("var", c_double_p),
    ]


class Panel(C.Structure):
    _fields_ = [
        ("T", C.c_size_t),
        ("N", C.c_size_t),
        ("dates", C.POINTER(C.c_char_p)),
        ("symbols", C.POINTER(C.c_char_p)),
        ("close", c_double_p),
        ("ret", c_double_p),
    ]


class HMMTrainOpts(C.Structure):
    _fields_ = [
        ("max_iter", C.c_uint),
        ("tol", C.c_double),
        ("n_restarts", C.c_uint),
        ("seed", C.c_uint64),
        ("alpha", C.c_double),
        ("init", C.c_int),
        ("verbose", C.c_int),
    ]


class GHMMTrainOpts(C.Structure):
    _fields_ = [
        ("max_iter", C.c_uint),
        ("tol", C.c_double),
        ("n_restarts", C.c_uint),
        ("seed", C.c_uint64),
        ("alpha", C.c_double),
        ("var_floor", C.c_double),
        ("mix_floor", C.c_double),
        ("init", C.c_int),
        ("verbose", C.c_int),
    ]


class HMMReport(C.Structure):
    _fields_ = [
        ("loglik", C.c_double),
        ("bic", C.c_double),
        ("aic", C.c_double),
        ("n_params", C.c_size_t),
        ("iterations", C.c_uint),
        ("converged", C.c_int),
        ("restarts", C.c_uint),
    ]


class GHMMReport(C.Structure):
    _fields_ = [
        ("loglik", C.c_double),
        ("bic", C.c_double),
        ("aic", C.c_double),
        ("n_params", C.c_size_t),
        ("iterations", C.c_uint),
        ("converged", C.c_int),
        ("restarts", C.c_uint),
    ]


class GHMMModelInfo(C.Structure):
    _fields_ = [
        ("n_states", C.c_size_t),
        ("loglik", C.c_double),
        ("bic", C.c_double),
        ("aic", C.c_double),
        ("n_params", C.c_size_t),
        ("converged", C.c_int),
    ]


class BTOpts(C.Structure):
    _fields_ = [("cost_bps", C.c_double), ("apply_costs", C.c_int)]


class BTMetrics(C.Structure):
    _fields_ = [
        ("total_return", C.c_double),
        ("cagr", C.c_double),
        ("ann_vol", C.c_double),
        ("sharpe", C.c_double),
        ("max_drawdown", C.c_double),
        ("calmar", C.c_double),
        ("avg_turnover", C.c_double),
        ("exposure", C.c_double),
        ("hit_rate", C.c_double),
        ("final_equity", C.c_double),
        ("n_days", C.c_size_t),
    ]


class BTResult(C.Structure):
    _fields_ = [("equity", c_double_p), ("strategy_ret", c_double_p), ("position", c_double_p)]


class Lattice(C.Structure):
    _fields_ = [
        ("L", C.c_int),
        ("periodic", C.c_int),
        ("kind", C.c_int),
        ("site", C.POINTER(C.c_uint8)),
        ("h", C.POINTER(C.c_uint8)),
        ("v", C.POINTER(C.c_uint8)),
    ]


class LatticeStats(C.Structure):
    _fields_ = [
        ("n_sites", C.c_size_t),
        ("n_occupied", C.c_size_t),
        ("n_clusters", C.c_size_t),
        ("max_cluster", C.c_size_t),
        ("s_max", C.c_double),
        ("susceptibility", C.c_double),
        ("occupancy", C.c_double),
        ("spanning_x", C.c_int),
        ("spanning_y", C.c_int),
        ("spanning", C.c_int),
    ]


class PercExp(C.Structure):
    _fields_ = [
        ("L", C.c_int),
        ("kind", C.c_int),
        ("periodic", C.c_int),
        ("reps", C.c_size_t),
        ("seed", C.c_uint64),
    ]


class CorrCritical(C.Structure):
    _fields_ = [
        ("theta_c", C.c_double),
        ("theta_half", C.c_double),
        ("chi_max", C.c_double),
        ("giant_at_c", C.c_double),
        ("censored_high", C.c_int),
        ("censored_low", C.c_int),
    ]


class PercPoint(C.Structure):
    _fields_ = [
        ("theta", C.c_double),
        ("n_edges", C.c_size_t),
        ("n_components", C.c_size_t),
        ("giant_fraction", C.c_double),
        ("susceptibility", C.c_double),
        ("mean_degree", C.c_double),
    ]


class Components(C.Structure):
    _fields_ = [
        ("n", C.c_size_t),
        ("n_components", C.c_size_t),
        ("giant", C.c_size_t),
        ("giant_fraction", C.c_double),
        ("susceptibility", C.c_double),
        ("labels", c_size_p),
        ("sizes", c_size_p),
    ]


def _proto(name, restype, argtypes):
    fn = getattr(lib, name)
    fn.restype = restype
    fn.argtypes = argtypes


_proto("hp_status_message", C.c_char_p, [C.c_int])
_proto("hp_version", C.c_char_p, [])
_proto("hp_free", None, [C.c_void_p])

# RNG
_proto("hp_rng_seed", None, [C.POINTER(RNG), C.c_uint64])
_proto("hp_rng_seed_stream", None, [C.POINTER(RNG), C.c_uint64, C.c_uint64])
_proto("hp_rng_uniform", C.c_double, [C.POINTER(RNG)])
_proto("hp_rng_normal", C.c_double, [C.POINTER(RNG)])

# Discrete HMM
_proto("hp_hmm_new", C.POINTER(HMM), [C.c_size_t, C.c_size_t])
_proto("hp_hmm_free", None, [C.POINTER(HMM)])
_proto("hp_hmm_validate", C.c_int, [C.POINTER(HMM)])
_proto("hp_hmm_forward_log", C.c_int, [C.POINTER(HMM), c_u32_p, C.c_size_t, c_double_p, c_double_p])
_proto("hp_hmm_backward_log", C.c_int, [C.POINTER(HMM), c_u32_p, C.c_size_t, c_double_p])
_proto("hp_hmm_posterior", C.c_int, [C.POINTER(HMM), c_u32_p, C.c_size_t, c_double_p, c_double_p])
_proto("hp_hmm_filter", C.c_int, [C.POINTER(HMM), c_u32_p, C.c_size_t, c_double_p, c_double_p])
_proto("hp_hmm_viterbi", C.c_int, [C.POINTER(HMM), c_u32_p, C.c_size_t, c_u32_p, c_double_p])
_proto("hp_hmm_train", C.c_int,
       [c_u32_p, C.c_size_t, C.c_size_t, C.c_size_t, C.POINTER(HMMTrainOpts),
        C.POINTER(HMM), C.POINTER(HMMReport)])
_proto("hp_hmm_save", C.c_int, [C.POINTER(HMM), C.c_char_p])
_proto("hp_hmm_load", C.POINTER(HMM), [C.c_char_p])

# Gaussian-mixture HMM
_proto("hp_ghmm_new", C.POINTER(GHMM), [C.c_size_t, C.c_size_t, C.c_size_t])
_proto("hp_ghmm_free", None, [C.POINTER(GHMM)])
_proto("hp_ghmm_validate", C.c_int, [C.POINTER(GHMM)])
_proto("hp_ghmm_forward_log", C.c_int, [C.POINTER(GHMM), c_double_p, C.c_size_t, c_double_p, c_double_p])
_proto("hp_ghmm_backward_log", C.c_int, [C.POINTER(GHMM), c_double_p, C.c_size_t, c_double_p])
_proto("hp_ghmm_posterior", C.c_int, [C.POINTER(GHMM), c_double_p, C.c_size_t, c_double_p, c_double_p])
_proto("hp_ghmm_filter", C.c_int, [C.POINTER(GHMM), c_double_p, C.c_size_t, c_double_p, c_double_p])
_proto("hp_ghmm_filter_step", C.c_int, [C.POINTER(GHMM), c_double_p, c_double_p, c_double_p])
_proto("hp_ghmm_viterbi", C.c_int, [C.POINTER(GHMM), c_double_p, C.c_size_t, c_u32_p, c_double_p])
_proto("hp_ghmm_train", C.c_int,
       [c_double_p, C.c_size_t, C.c_size_t, C.c_size_t, C.c_size_t,
        C.POINTER(GHMMTrainOpts), C.POINTER(GHMM), C.POINTER(GHMMReport)])
_proto("hp_ghmm_select", C.c_int,
       [c_double_p, C.c_size_t, C.c_size_t, C.c_size_t, C.c_size_t, C.c_size_t,
        C.POINTER(GHMMTrainOpts), C.POINTER(GHMM), c_size_p, C.POINTER(GHMMModelInfo)])
_proto("hp_ghmm_save", C.c_int, [C.POINTER(GHMM), C.c_char_p])
_proto("hp_ghmm_load", C.POINTER(GHMM), [C.c_char_p])

# Panel and backtest
_proto("hp_panel_load_csv", C.POINTER(Panel), [C.c_char_p, C.POINTER(C.c_int)])
_proto("hp_panel_free", None, [C.POINTER(Panel)])
_proto("hp_panel_compute_returns", C.c_int, [C.POINTER(Panel)])
_proto("hp_panel_save_csv", C.c_int, [C.POINTER(Panel), C.c_char_p])
_proto("hp_panel_ew_returns", C.c_int, [C.POINTER(Panel), c_double_p])
_proto("hp_panel_return_series", C.c_int, [C.POINTER(Panel), C.c_size_t, c_double_p])
_proto("hp_backtest", C.c_int,
       [c_double_p, C.c_size_t, c_double_p, C.POINTER(BTOpts),
        C.POINTER(BTResult), C.POINTER(BTMetrics)])
_proto("hp_bt_result_free", None, [C.POINTER(BTResult)])

# Percolation
_proto("hp_lattice_new", C.POINTER(Lattice), [C.c_int, C.c_int, C.c_int])
_proto("hp_lattice_free", None, [C.POINTER(Lattice)])
_proto("hp_lattice_bernoulli", C.c_int, [C.POINTER(Lattice), C.c_double, C.POINTER(RNG)])
_proto("hp_lattice_analyze", C.c_int, [C.POINTER(Lattice), C.POINTER(LatticeStats)])
_proto("hp_lattice_cluster_sizes", c_size_p, [C.POINTER(Lattice), c_size_p])
_proto("hp_perc_spanning_prob", C.c_double, [C.POINTER(PercExp), C.c_double])
_proto("hp_perc_crossing", C.c_int,
       [C.POINTER(PercExp), C.c_double, C.c_double, C.c_uint, c_double_p, c_double_p])
_proto("hp_perc_er_giant_theory", C.c_double, [C.c_double])
_proto("hp_perc_bond_threshold", C.c_double, [C.c_double])
_proto("hp_perc_config_giant", C.c_int,
       [c_double_p, C.c_size_t, C.c_double, c_double_p, c_double_p, c_double_p])
_proto("hp_corr_matrix", None, [c_double_p, C.c_size_t, C.c_size_t, c_double_p])
_proto("hp_corr_critical", C.c_int,
       [c_double_p, C.c_size_t, C.c_double, C.c_double, C.c_size_t, C.POINTER(CorrCritical)])
_proto("hp_corr_sweep", C.c_int,
       [c_double_p, C.c_size_t, c_double_p, C.c_size_t, C.POINTER(PercPoint)])
_proto("hp_graph_from_corr", C.c_void_p, [c_double_p, C.c_size_t, C.c_double])
_proto("hp_graph_components", C.c_int, [C.c_void_p, C.c_double, C.POINTER(Components)])
_proto("hp_graph_free", None, [C.c_void_p])
_proto("hp_components_free", None, [C.POINTER(Components)])
_proto("hp_graph_edge_count", C.c_size_t, [C.c_void_p])
