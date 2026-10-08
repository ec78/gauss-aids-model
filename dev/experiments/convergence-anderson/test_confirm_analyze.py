"""test_confirm_analyze.py -- unit + integration tests for confirm_analyze.py.

Run from this directory:   python -m unittest -v test_confirm_analyze
No solver; synthetic tests write small cell/truth files into temp dirs; the
integration test reads only development seeds 1-30 (pilot_out/, e0_out/).
"""
import filecmp
import json
import math
import os
import shutil
import tempfile
import unittest

import numpy as np

import confirm_analyze as ca

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = ca.load_config()
META = ca.read_meta_names(os.path.join(HERE, "pilot_out"))

# synthetic layout: model 0 (AIDS), nIntFull 2, n1 2 -> ng = 2+2+1+1 = 6, n = 3
NIF, N1, NG, N = 2, 2, 6, 3


def truth_params(seed):
    """Synthetic truth satisfying homogeneity (n-th gamma row = -sum of the
    n1 gamma rows) and adding-up (n-th column = -row sums, +1 on the
    constant), like _quaidsSyntheticDGP's trueParams, so that an endpoint
    equal to the truth has pfErr exactly 0."""
    rng = np.random.default_rng(seed)
    full = np.zeros((NG + 1, N))
    red = np.round(rng.normal(size=(NG, N1)), 3)          # reduced raw layout
    full[:NIF, :N1] = red[:NIF]
    full[NIF:NIF + N1, :N1] = red[NIF:NIF + N1]
    full[NIF + N1, :N1] = -red[NIF:NIF + N1].sum(axis=0)
    full[NIF + N1 + 1:, :N1] = red[NIF + N1:]
    full[:, N1] = -full[:, :N1].sum(axis=1)
    full[0, N1] += 1.0
    return full


def raw_from_truth(tp, shift=0.0):
    tr = ca.true_raw(tp, NIF, N1)[:, :N1].copy()
    return tr + shift


def write_dataset(d, cells, truths, model=0):
    """cells: list of dict(seed, method, start, conv, crit, raw, [bS], [wall], [bsOk]).
    Stored error columns are computed exactly as the harness would, so the
    truth gate passes unless a test deliberately corrupts them."""
    os.makedirs(d, exist_ok=True)
    shutil.copy(os.path.join(HERE, "pilot_out", "columns.txt"), os.path.join(d, "columns.txt"))
    with open(os.path.join(d, f"cells_m{model}.csv"), "w") as fh:
        for c in cells:
            tp = truths[c["seed"]]
            meta = {k: float("nan") for k in META}
            raw = c["raw"]
            bS = np.array(c.get("bS", np.zeros((NG + 1, N))))
            meta.update(model=model, seed=c["seed"], method=c["method"], start=c["start"],
                        conv=c["conv"], critExact=c["crit"], critTrue=-5.0, ng=NG, n1=N1,
                        bsOk=c.get("bsOk", 1), wallTime=c.get("wall", 1.0), status=0)
            meta["rawStructErr"] = ca.raw_struct_err(raw, tp, NIF, model)
            meta["bsErr"] = float(np.abs(bS - tp).max())
            vals = [meta[k] for k in META] + list(raw.ravel(order="F")) + list(bS.ravel(order="F"))
            fh.write(",".join("." if (isinstance(v, float) and math.isnan(v)) else repr(float(v)) for v in vals) + "\n")
    with open(os.path.join(d, f"truth_m{model}.csv"), "w") as fh:
        for s in sorted(truths):
            vals = [model, s, NG, N1, -5.0, -9.0] + list(truths[s].ravel(order="F"))
            fh.write(",".join(repr(float(v)) for v in vals) + "\n")


def cell(seed, method, start, conv=1, crit=-5.0, raw=None, **kw):
    d = dict(seed=seed, method=method, start=start, conv=conv, crit=crit, raw=raw)
    d.update(kw)
    return d


def full_grid(seeds, truths, override=None):
    """Every (seed, method 1-3, start 1-16) cell; default endpoint = truth
    (correct), objective -6; override(seed, method, start) -> dict updates."""
    out = []
    for s in seeds:
        for m in (1, 2, 3):
            for st in range(1, 17):
                c = cell(s, m, st, conv=1, crit=-6.0, raw=raw_from_truth(truths[s]))
                if override:
                    c.update(override(s, m, st) or {})
                out.append(c)
    return out


def cfg_with_seeds(lo, hi, tmp):
    cfg = dict(CFG)
    cfg["dev_seeds"] = [lo, hi]
    p = os.path.join(tmp, "cfg.json")
    with open(p, "w") as fh:
        json.dump(cfg, fh)
    return p


class TestSelectionUnits(unittest.TestCase):
    def mk(self, method, start, conv=1, crit=-5.0, raw=None):
        return dict(method=method, start=start, conv=conv, critExact=crit,
                    raw=np.zeros((NG, N1)) if raw is None else raw)

    def test_eligibility(self):
        self.assertTrue(ca.is_eligible(self.mk(2, 1)))
        self.assertFalse(ca.is_eligible(self.mk(2, 1, conv=0)), "non-converged must be ineligible")
        self.assertFalse(ca.is_eligible(self.mk(2, 1, crit=float("nan"))), "NaN objective ineligible")
        self.assertFalse(ca.is_eligible(self.mk(2, 1, crit=float("inf"))), "infinite objective ineligible")
        bad = np.zeros((NG, N1)); bad[3, 1] = np.inf
        self.assertFalse(ca.is_eligible(self.mk(2, 1, raw=bad)), "non-finite endpoint ineligible")

    def test_objective_selection_direction(self):
        pool = [self.mk(2, 1, crit=-7.84), self.mk(2, 2, crit=-5.47), self.mk(3, 1, crit=-6.0)]
        best, _ = ca.select(pool)
        self.assertEqual((best["method"], best["start"]), (2, 2), "largest critExact (-ln det) must win")

    def test_nonconverged_better_objective_excluded(self):
        pool = [self.mk(2, 1, conv=0, crit=-1.0), self.mk(3, 4, crit=-5.0)]
        best, elig = ca.select(pool)
        self.assertEqual((best["method"], best["start"]), (3, 4))
        self.assertEqual(len(elig), 1)

    def test_exact_ties(self):
        pool = [self.mk(3, 1, crit=-5.0), self.mk(2, 9, crit=-5.0), self.mk(2, 4, crit=-5.0)]
        best, _ = ca.select(pool)
        self.assertEqual((best["method"], best["start"]), (2, 4), "ties: cold LM before continuation, then lowest start")
        pool = [self.mk(3, 7, crit=-5.0), self.mk(3, 2, crit=-5.0)]
        best, _ = ca.select(pool)
        self.assertEqual(best["start"], 2)

    def test_nonfinite_endpoint_never_selected(self):
        bad = np.zeros((NG, N1)); bad[0, 0] = np.nan
        pool = [self.mk(2, 1, crit=-1.0, raw=bad), self.mk(2, 2, crit=-9.0)]
        best, _ = ca.select(pool)
        self.assertEqual(best["start"], 2)

    def test_no_result(self):
        best, elig = ca.select([self.mk(2, 1, conv=0), self.mk(3, 1, crit=float("nan"))])
        self.assertIsNone(best)
        self.assertEqual(elig, [])
        self.assertEqual(ca.classify(None, 1.0), "none")

    def test_classification_boundaries(self):
        self.assertEqual(ca.classify(1.0, 1.0), "correct", "pfErr == cutoff counts as correct")
        self.assertEqual(ca.classify(np.nextafter(1.0, 2.0), 1.0), "wrong")
        self.assertEqual(ca.classify(0.0, 1.0), "correct")

    def test_pferr_formula_hand_example(self):
        tp = truth_params(1)
        self.assertEqual(ca.pf_err(raw_from_truth(tp), tp, NIF), 0.0, "endpoint == truth must give pfErr 0")
        # uniform +3 shift of every reduced coefficient (hand-derived):
        #   n-th gamma row cols 1-2 = -(sum of the two gamma rows) -> shift -6;
        #   n-th column entries = -(row sum of two shifted entries) -> shift -6;
        #   their intersection = -(sum of the two -6-shifted entries) -> shift +12  => max error 12
        self.assertAlmostEqual(ca.pf_err(raw_from_truth(tp, 3.0), tp, NIF), 12.0, places=12)
        # explicit construction: truth satisfying homogeneity/adding-up exactly
        full = np.zeros((NG + 1, N))
        full[:NIF, :N1] = [[0.3, 0.2], [0.1, -0.1]]
        full[NIF:NIF + N1, :N1] = [[0.05, -0.02], [-0.02, 0.04]]
        full[NIF + N1, :N1] = -full[NIF:NIF + N1, :N1].sum(axis=0)
        full[NIF + N1 + 1:, :N1] = [[0.07, -0.03], [0.5, 0.6]]
        full[:, N1] = -full[:, :N1].sum(axis=1)
        full[0, N1] += 1.0
        rawx = ca.true_raw(full, NIF, N1)[:, :N1]
        self.assertTrue(np.array_equal(ca.unprojected_full(rawx, NIF), full), "recovery must invert exactly")
        self.assertEqual(ca.pf_err(rawx, full, NIF), 0.0)
        r2 = rawx.copy(); r2[-1, 0] += 50.0          # u row: excluded
        self.assertEqual(ca.pf_err(r2, full, NIF), 0.0)
        r3 = rawx.copy(); r3[NIF, 0] += 0.25         # gamma: enters gamma row, n-th row and n-th column
        self.assertAlmostEqual(ca.pf_err(r3, full, NIF), 0.25, places=12)

    def test_layout_assertion(self):
        self.assertEqual(ca.n_int_full(8, 4, 0), 2)
        self.assertEqual(ca.n_int_full(9, 4, 1), 2)
        with self.assertRaises(ca.IntegrityError):
            ca.n_int_full(6, 4, 1)

    def test_a2_interval_thresholds(self):
        lo161, _ = ca.cp_two_sided(161, 170)
        lo160, _ = ca.cp_two_sided(160, 170)
        lo154, _ = ca.cp_two_sided(154, 170)
        lo153, _ = ca.cp_two_sided(153, 170)
        self.assertGreaterEqual(lo161, 0.90); self.assertLess(lo160, 0.90)
        self.assertGreaterEqual(lo154, 0.85); self.assertLess(lo153, 0.85)
        self.assertAlmostEqual(lo161, 0.9019, places=4)
        self.assertAlmostEqual(lo154, 0.8517, places=4)
        self.assertEqual(ca.cp_two_sided(0, 30)[0], 0.0)
        self.assertEqual(ca.cp_two_sided(30, 30)[1], 1.0)
        self.assertAlmostEqual(ca.cp_two_sided(1, 170)[1], 0.0323, places=4)

    def test_a3_descriptive_pvalue(self):
        self.assertAlmostEqual(ca.mcnemar_one_sided(9, 0), 0.5 ** 9)
        self.assertAlmostEqual(ca.mcnemar_one_sided(17, 0), 0.5 ** 17)
        self.assertEqual(ca.mcnemar_one_sided(0, 0), 1.0)


class TestEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="confirm_test_")
        self.seeds = [1, 2, 3, 4]
        self.truths = {s: truth_params(s) for s in self.seeds}
        self.cfgp = cfg_with_seeds(1, 4, self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def build(self, override, name="d"):
        d = os.path.join(self.tmp, name)
        write_dataset(d, full_grid(self.seeds, self.truths, override), self.truths)
        return d

    def test_a4_counting_and_classes(self):
        # seed 1: all correct; seed 2: cold LM selects a WRONG endpoint (best objective);
        # seed 3: nothing eligible in L16 (no result); seed 4: continuation-only correct basin wins in L16+C16
        def ov(s, m, st):
            t = self.truths[s]
            if s == 2 and m == 2 and st == 5:
                return dict(crit=-1.0, raw=raw_from_truth(t, 3.0))
            if s == 3 and m == 2:
                return dict(conv=0)
            if s == 4 and m == 2:
                return dict(crit=-8.0, raw=raw_from_truth(t, 5.0))
            if s == 4 and m == 3 and st == 7:
                return dict(crit=-5.5)
            return None
        d = self.build(ov)
        res, rows = ca.run("dev", d, d, os.path.join(self.tmp, "out"), cfg_path=self.cfgp)
        L16 = res["models"]["0"]["configurations"]["L16"]
        self.assertEqual(L16["classesPf"], {"none": 1, "wrong": 2, "correct": 1})
        self.assertEqual(L16["A4"]["W"], 2, "no-result seed 3 must NOT count toward W")
        self.assertEqual(L16["A4"]["wrongSeeds"], [2, 4])
        self.assertEqual(L16["noResultSeeds"], [3])
        pooled = res["models"]["0"]["configurations"]["L16+C16"]
        self.assertEqual(pooled["classesPf"], {"none": 0, "wrong": 1, "correct": 3})
        r4 = [r for r in rows if r["seed"] == 4 and r["config"] == "L16+C16"][0]
        self.assertEqual((r4["selMethod"], r4["selStart"]), (3, 7))
        g = res["models"]["0"]["gates"]
        self.assertFalse(g["A4"]["pass"]); self.assertEqual(g["A4"]["W"], 2)

    def test_deterministic_output(self):
        d = self.build(lambda s, m, st: dict(crit=-6.0 + 0.001 * st) if m == 3 else None)
        o1, o2 = os.path.join(self.tmp, "o1"), os.path.join(self.tmp, "o2")
        ca.run("dev", d, d, o1, cfg_path=self.cfgp)
        ca.run("dev", d, d, o2, cfg_path=self.cfgp)
        for f in ("results.json", "selected.csv", "report.md"):
            self.assertTrue(filecmp.cmp(os.path.join(o1, f), os.path.join(o2, f), shallow=False), f)

    def test_dev_guard_refuses_out_of_range_seed(self):
        d = self.build(None)
        with open(os.path.join(d, "cells_m0.csv")) as fh:
            first = fh.readline()
        parts = first.split(",")
        parts[1] = "31"
        with open(os.path.join(d, "cells_m0.csv"), "a") as fh:
            fh.write(",".join(parts))
        with self.assertRaises(ca.GuardError):
            ca.run("dev", d, d, os.path.join(self.tmp, "out"), cfg_path=self.cfgp)

    def test_confirm_mode_requires_valid_manifest(self):
        d = self.build(None)
        with self.assertRaises(ca.GuardError):
            ca.run("confirm", d, d, os.path.join(self.tmp, "out"), cfg_path=self.cfgp)
        man = os.path.join(self.tmp, "m.json")
        with open(man, "w") as fh:
            json.dump({"repo_root_rel": ".", "files": {"cfg.json": "0" * 64}}, fh)
        with self.assertRaises(ca.GuardError):
            ca.run("confirm", d, d, os.path.join(self.tmp, "out"), manifest=man, cfg_path=self.cfgp)

    def test_empty_input_is_an_integrity_failure(self):
        d = os.path.join(self.tmp, "empty")
        os.makedirs(d)
        shutil.copy(os.path.join(HERE, "pilot_out", "columns.txt"), os.path.join(d, "columns.txt"))
        with self.assertRaises(ca.IntegrityError):
            ca.run("dev", d, d, os.path.join(self.tmp, "out"), cfg_path=self.cfgp)

    def test_truth_gate_detects_mismatch(self):
        d = self.build(None)
        p = os.path.join(d, "cells_m0.csv")
        with open(p) as fh:
            lines = fh.read().splitlines()
        parts = lines[0].split(",")
        parts[META.index("rawStructErr")] = "0.123"
        lines[0] = ",".join(parts)
        with open(p, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        with self.assertRaises(ca.IntegrityError):
            ca.run("dev", d, d, os.path.join(self.tmp, "out"), cfg_path=self.cfgp)

    def test_incomplete_grid_detected(self):
        d = self.build(None)
        p = os.path.join(d, "cells_m0.csv")
        with open(p) as fh:
            lines = fh.read().splitlines()[:-1]
        with open(p, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        with self.assertRaises(ca.IntegrityError):
            ca.run("dev", d, d, os.path.join(self.tmp, "out"), cfg_path=self.cfgp)


class TestReproducesE0(unittest.TestCase):
    """Integration: development seeds 1-30 only (pilot_out/, e0_out/)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="confirm_e0_")
        cls.res, cls.rows = ca.run("dev", os.path.join(HERE, "pilot_out"), os.path.join(HERE, "e0_out"),
                                   os.path.join(cls.tmp, "out"))
        cls.cells = ca.load_cells(os.path.join(HERE, "pilot_out"), (1, 30))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cell_level_pferr_exact(self):
        import csv
        ref = {}
        with open(os.path.join(HERE, "e0_out", "e0_cells_pferr.csv")) as fh:
            for r in csv.DictReader(fh):
                k = (int(r["model"]), int(r["seed"]), int(r["method"]), int(r["start"]))
                ref[k] = float(r["pfErr"]) if r["pfErr"] not in ("", "nan") else float("nan")
        truth = ca.load_truth(os.path.join(HERE, "e0_out"), (1, 30))
        n = 0
        for c in self.cells:
            k = (c["model"], c["seed"], c["method"], c["start"])
            if np.all(np.isfinite(c["raw"])):
                v = ca.pf_err(c["raw"], truth[(c["model"], c["seed"])]["trueParams"],
                              ca.n_int_full(c["ng"], c["n1"], c["model"]))
                self.assertEqual(v, ref[k], k)
            else:
                self.assertTrue(math.isnan(ref[k]), k)
            n += 1
        self.assertEqual(n, 2880)

    def test_class_tables_match_e0_report(self):
        # protocol section 11.3 (pfErr / raw / bS as none/wrong/correct) and 11.4 (W)
        expected = {
            ("0", "A16"): ("2/7/21", "2/6/22", "2/9/19"), ("0", "L1"): ("0/0/30", "0/0/30", "0/4/26"),
            ("0", "C1"): ("0/0/30", "0/0/30", "0/4/26"), ("0", "L1+C1"): ("0/0/30", "0/0/30", "0/3/27"),
            ("0", "L16"): ("0/0/30", "0/0/30", "0/3/27"), ("0", "C16"): ("0/0/30", "0/0/30", "0/4/26"),
            ("0", "L16+C1"): ("0/0/30", "0/0/30", "0/3/27"), ("0", "L16+C16"): ("0/0/30", "0/0/30", "0/3/27"),
            ("1", "A16"): ("1/16/13", "1/16/13", "1/19/10"), ("1", "L1"): ("4/3/23", "4/3/23", "4/9/17"),
            ("1", "C1"): ("0/3/27", "0/3/27", "0/9/21"), ("1", "L1+C1"): ("0/1/29", "0/1/29", "0/7/23"),
            ("1", "L16"): ("0/1/29", "0/1/29", "0/7/23"), ("1", "C16"): ("0/1/29", "0/1/29", "0/7/23"),
            ("1", "L16+C1"): ("0/1/29", "0/1/29", "0/7/23"), ("1", "L16+C16"): ("0/0/30", "0/0/30", "0/6/24"),
        }
        f = lambda d: f"{d['none']}/{d['wrong']}/{d['correct']}"
        for (m, name), (pf, raw, bs) in expected.items():
            cr = self.res["models"][m]["configurations"][name]
            self.assertEqual((f(cr["classesPf"]), f(cr["classesRaw"]), f(cr["classesBs"])), (pf, raw, bs), (m, name))
        W = {("0", "A16"): 7, ("0", "L16"): 0, ("0", "C16"): 0, ("0", "L16+C1"): 0, ("0", "L16+C16"): 0,
             ("1", "A16"): 16, ("1", "L16"): 1, ("1", "C16"): 1, ("1", "L16+C1"): 1, ("1", "L16+C16"): 0}
        for (m, name), w in W.items():
            self.assertEqual(self.res["models"][m]["configurations"][name]["A4"]["W"], w, (m, name))

    def test_selection_matches_e0_rescore(self):
        """Same selected (method, start) as e0_rescore.py's selector, every seed/config."""
        import pandas as pd
        e0 = pd.read_csv(os.path.join(HERE, "e0_out", "e0_cells_pferr.csv"))
        for r in self.rows:
            g = e0[(e0.model == r["model"]) & (e0.seed == r["seed"])]
            conf = CFG["configurations"][r["config"]]
            mask = np.zeros(len(g), bool)
            for meth, (a, b) in conf.items():
                mask |= ((g.method == int(meth)) & (g.start >= a) & (g.start <= b)).values
            c = g[mask & (g.conv == True).values].sort_values(["method", "start"])
            if not len(c):
                self.assertIsNone(r["selMethod"]); continue
            best = c.loc[c.critExact.idxmax()]
            self.assertEqual((r["selMethod"], r["selStart"]), (int(best.method), int(best.start)), (r["model"], r["seed"], r["config"]))

    def test_dev_a3_descriptive_and_reconciliation(self):
        a3 = self.res["models"]["0"]["A3_descriptive"]
        self.assertEqual((a3["b"], a3["c"], a3["ties"]), (9, 0, 21))
        self.assertEqual(a3["improvedSeeds"], [5, 6, 10, 13, 15, 16, 17, 24, 29])
        q = self.res["models"]["1"]["A3_descriptive"]
        self.assertEqual((q["b"], q["c"], q["ties"]), (17, 0, 13))

    def test_dev_gates(self):
        self.assertFalse(self.res["gatesBinding"], "dev-mode gate verdicts must be marked non-binding")
        for m in ("0", "1"):
            g = self.res["models"][m]["gates"]
            self.assertEqual(g["A2"]["k"], 30); self.assertEqual(g["A4"]["W"], 0)

    def test_no_selected_nonfinite_and_ties(self):
        self.assertTrue(all(r["pfErr"] is not None for r in self.rows if r["selMethod"] is not None))
        self.assertEqual(self.res["models"]["1"]["configurations"]["C16"]["maxExactTies"], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
