"""hpmarket · Python bindings for the C HMM / percolation trading library."""
from ._capi import HPError
from .core import (
    BOND, INIT_QUANTILE, INIT_RANDOM, SITE,
    DiscreteHMM, GaussianHMM, Lattice, Panel, RNG,
    backtest, bond_threshold, config_giant, corr_critical, corr_matrix,
    corr_sweep, er_giant, graph_giant, perc_crossing, perc_spanning_prob,
    select_ghmm, vol_target_signal, momentum_signal, state_exposure_signal,
    multiply_signals, constant_signal,
)

__all__ = [
    "HPError", "RNG", "DiscreteHMM", "GaussianHMM", "Panel", "Lattice",
    "select_ghmm", "corr_matrix", "corr_critical", "corr_sweep", "graph_giant",
    "perc_crossing", "perc_spanning_prob", "er_giant", "config_giant",
    "bond_threshold", "backtest", "SITE", "BOND", "INIT_QUANTILE", "INIT_RANDOM",
    "vol_target_signal", "momentum_signal", "state_exposure_signal",
    "multiply_signals", "constant_signal",
]
from ._capi import lib as _lib

__version__ = _lib.hp_version().decode("utf-8")
