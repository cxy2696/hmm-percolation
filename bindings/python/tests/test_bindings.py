"""Binding tests.  Run from the project root:

    make shared
    PYTHONPATH=bindings/python python3 -m unittest discover -s bindings/python/tests -v
"""
from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path

import hpmarket as hp

ROOT = Path(__file__).resolve().parents[3]


class TestDiscreteHMM(unittest.TestCase):
    def test_uniform_likelihood(self):
        h = hp.DiscreteHMM(2, 2)
        _, ll = h.forward([0, 1, 0])
        self.assertAlmostEqual(ll, math.log(1.0 / 8.0), places=12)
        path, lp = h.viterbi([0, 1, 0])
        self.assertEqual(len(path), 3)
        self.assertTrue(math.isfinite(lp))

    def test_train_and_roundtrip(self):
        h = hp.DiscreteHMM(2, 3)
        obs = ([0] * 120 + [2] * 120) * 2
        rep = h.train(obs, max_iter=80, restarts=2, seed=7)
        self.assertTrue(math.isfinite(rep["loglik"]))
        self.assertEqual(h.n_states, 2)
        h.save(ROOT / "build" / "test_tmp" / "py_hmm.txt")
        h2 = hp.DiscreteHMM.load(ROOT / "build" / "test_tmp" / "py_hmm.txt")
        _, ll = h2.forward(obs)
        self.assertAlmostEqual(ll, rep["loglik"], places=8)


class TestGaussianHMM(unittest.TestCase):
    def test_recover_two_regimes(self):
        rng = hp.RNG(11, 0)
        x = []
        for _ in range(400):
            x.append([-2.0 + 0.6 * rng.normal()])
        for _ in range(400):
            x.append([2.0 + 0.6 * rng.normal()])
        m = hp.GaussianHMM(2, 1)
        rep = m.train(x, max_iter=100, restarts=3, seed=3)
        self.assertLess(m.state_mean(0), -1.5)
        self.assertGreater(m.state_mean(1), 1.5)
        self.assertTrue(rep["converged"])
        post, ll = m.posterior(x)
        self.assertEqual(len(post), len(x))
        self.assertAlmostEqual(ll, rep["loglik"], places=6)
        m.save(ROOT / "build" / "test_tmp" / "py_ghmm2.txt")
        m2 = hp.GaussianHMM.load(ROOT / "build" / "test_tmp" / "py_ghmm2.txt")
        _, ll2 = m2.forward(x)
        self.assertAlmostEqual(ll2, ll, places=8)


class TestPanel(unittest.TestCase):
    def test_load_synthetic(self):
        p = hp.Panel.load(ROOT / "data" / "synthetic_panel.csv")
        self.assertEqual(p.T, 3000)
        self.assertEqual(p.N, 10)
        self.assertEqual(p.symbols[0], "A1")
        self.assertEqual(p.dates[0], "2012-01-02")
        self.assertTrue(all(c > 0 for c in p.close_prices[:100]))
        ew = p.ew_returns()
        self.assertTrue(all(math.isfinite(v) for v in ew[1:]))
        expected = math.log(p.close_prices[10] / p.close_prices[0])
        self.assertAlmostEqual(p.return_series(0)[1], expected, places=12)


class TestPercolation(unittest.TestCase):
    def test_lattice(self):
        lat = hp.Lattice(48, "bond")
        rng = hp.RNG(5, 0)
        lat.bernoulli(0.6, rng)
        st = lat.analyze()
        self.assertEqual(st["n_sites"], 48 * 48)
        self.assertGreater(st["s_max"], 0.5)
        self.assertTrue(st["spanning"])
        self.assertEqual(sum(lat.cluster_sizes()), 48 * 48)

    def test_crossing(self):
        pc, sp = hp.perc_crossing(32, reps=40, seed=9, iters=10)
        self.assertLess(abs(pc - 0.5), 0.08)

    def test_meanfield(self):
        self.assertAlmostEqual(hp.er_giant(2.0), 0.7968121300200199, places=6)
        self.assertAlmostEqual(hp.bond_threshold(2.0), 0.5, places=12)
        pk = [math.exp(-2) * 2 ** k / math.factorial(k) for k in range(25)]
        s = sum(pk)
        pk = [v / s for v in pk]
        r = hp.config_giant(pk, 0.6)
        self.assertAlmostEqual(r["r0"], 1.2, places=3)
        self.assertGreater(r["giant"], 0.0)


class TestCorrBacktest(unittest.TestCase):
    def test_corr_block(self):
        n = 20
        corr = [[1.0 if i == j else (0.8 if (i < 10) == (j < 10) else 0.0)
                 for j in range(n)] for i in range(n)]
        crit = hp.corr_critical(corr)
        self.assertAlmostEqual(crit["theta_half"], 0.8, delta=0.03)
        g = hp.graph_giant(corr, 0.5)
        self.assertEqual(g["n_components"], 2)
        self.assertAlmostEqual(g["giant_fraction"], 0.5, places=12)

    def test_backtest_timing(self):
        ret = [math.nan] + [0.0] * 4 + [0.01, 0.0, 0.0]
        sig = [0, 0, 0, 0, 1, 0, 0, 0]
        out = hp.backtest(ret, sig, cost_bps=0.0)
        self.assertAlmostEqual(out["equity"][5], 1.01, places=12)
        self.assertAlmostEqual(out["total_return"], 0.01, places=12)
        no_future = hp.backtest(ret, [0, 0, 0, 0, 0, 1, 0, 0], cost_bps=0.0)
        self.assertAlmostEqual(no_future["equity"][7], 1.0, places=12)


class TestAgainstCpp(unittest.TestCase):
    """Recompute the C++ pipeline's training likelihood from Python."""

    def test_model_loglik(self):
        out = ROOT / "out" / "regime"
        metrics_path = out / "metrics.json"
        model_path = out / "model_train.txt"
        if not metrics_path.exists() or not model_path.exists():
            self.skipTest("run ./build/hp_regime first")
        meta = json.loads(metrics_path.read_text())
        p = hp.Panel.load(ROOT / "data" / "synthetic_panel.csv")
        ret = p.ew_returns()
        w = int(meta["config"]["vol_window"])
        start = int(meta["data"]["start"])
        train_end = int(meta["data"]["train_end"])
        x = []
        for t in range(start, len(ret)):
            window = ret[t - w + 1:t + 1]
            mean = sum(window) / w
            var = sum((v - mean) ** 2 for v in window) / (w - 1)
            vol = max(math.sqrt(var * 252.0), 1e-4)
            x.append([ret[t] * 100.0, math.log(vol * 100.0)])
        mean = meta["features"]["mean"]
        sd = meta["features"]["sd"]
        z = [[(row[d] - mean[d]) / sd[d] for d in range(2)] for row in x]
        z_train = z[: train_end - start]
        m = hp.GaussianHMM.load(model_path)
        _, ll = m.forward(z_train)
        self.assertAlmostEqual(ll, meta["model"]["loglik"], delta=2e-5)

    def test_connectedness_tail(self):
        csv_path = ROOT / "out" / "connectedness" / "connectedness.csv"
        if not csv_path.exists():
            self.skipTest("run ./build/hp_connectedness first")
        with csv_path.open() as f:
            rows = list(csv.DictReader(f))
        p = hp.Panel.load(ROOT / "data" / "synthetic_panel.csv")
        rets = [p.return_series(j) for j in range(p.N)]
        T = p.T
        w = 60
        window = [[rets[j][t] for j in range(p.N)] for t in range(T - w, T)]
        corr = hp.corr_matrix(window)
        crit = hp.corr_critical(corr)
        self.assertAlmostEqual(float(rows[-1]["theta_c"]), crit["theta_half"], places=10)
        self.assertAlmostEqual(float(rows[-1]["theta_chi"]), crit["theta_peak"], places=10)


if __name__ == "__main__":
    unittest.main()
