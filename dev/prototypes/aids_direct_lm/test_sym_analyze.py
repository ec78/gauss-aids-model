"""test_sym_analyze.py -- tests for sym_analyze.py (SYMMETRY_PROTOCOL.md v1).

I6: the pfErr/rawStructErr used here equal confirm_analyze.pf_err /
raw_struct_err on the stored E0 rows (seeds 1-30, AIDS), read only. The
frozen module is imported with bytecode writing disabled, so nothing is
written under dev/experiments/convergence-anderson/.

The other tests run the analyzer on small hand-built run directories in a
temporary directory: selection recomputation, INVALID on a harness
disagreement, the seed guard, determinism compare, byte-identical
reruns, interim mode without truth, the real-data G6 path, the no-overwrite
and freeze-precondition guards, an incomplete determinism re-solve, Clopper-Pearson
thresholds and the section 10 decision table. No solver; no protocol seed.

Run from this directory:  python -m unittest -v test_sym_analyze
"""
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sym_analyze as sa  # noqa: E402

FROZEN = os.path.normpath(os.path.join(HERE, "..", "..", "experiments", "convergence-anderson"))


def build_run(d, block, seeds, ng=8, n1=4, nobs=3000, truth=True, polish=False):
    """A consistent harness-format run: every start lands on the symmetric
    truth (one basin); start 3 has retcode 2; H's G is slightly asymmetric."""
    nif = ng - n1 - 2
    os.makedirs(d, exist_ok=True)
    cells = ["block,seed,arm,start,retcode,eligible,critExact,gnorm,iterations,rejected,covInvFail,"
             "stepInvFail,infeasTrial,rhoReject,lmAccepted,lambdaFinal,wallSec"]
    raws = ["block,seed,arm,start,row,col,value"]
    seedl = [",".join(["block", "seed", "arm", "ok", "selStart", "nEligible", "crit", "pfErr", "rawStructErr",
                       "classPf", "asym", "addUpErr", "homogErr", "exactTies", "margin", "nearTie", "iterations",
                       "wallSec", "N1", "N2", "LR", "condGN", "critAbsErr"])]
    trl = ["block,seed,row,col,value"]
    f = sa.fmt
    for seed in seeds:
        rng = np.random.default_rng(seed + 17)
        R = np.round(rng.normal(size=(ng, n1)), 1)
        G = R[nif:nif + n1]
        R[nif:nif + n1] = np.triu(G) + np.triu(G, 1).T
        tp = sa.unprojected_full(R, nif)
        rawArm = {"H": R.copy(), "S": R.copy(), "P": R.copy()}
        rawArm["H"][nif, 1] += 1e-3
        crits = {"H": [-5 + 1e-3 * s for s in range(1, 17)], "S": [-5.001 + 1e-3 * s for s in range(1, 17)]}
        for arm in ("H", "S"):
            for s in range(1, 17):
                rc = 2 if s == 3 else 0
                cells.append(",".join([block, f(seed), arm, f(s), f(rc), f(int(rc == 0)), f(crits[arm][s - 1]),
                                       "0.0001", "7", "0", "0", "0", "0", "0", "6", "0.001", "."]))
                for r in range(ng):
                    for c in range(n1):
                        raws.append(",".join([block, f(seed), arm, f(s), f(r + 1), f(c + 1), f(rawArm[arm][r, c])]))
        for r in range(ng):
            for c in range(n1):
                raws.append(",".join([block, f(seed), "P", "0", f(r + 1), f(c + 1), f(rawArm["P"][r, c])]))
        ch, cs = crits["H"][15], crits["S"][15]
        cp = cs + 1e-9
        asymH = float(np.abs(rawArm["H"][nif:nif + n1] - rawArm["H"][nif:nif + n1].T).max())
        common = ["0", "0"]
        seedl.append(",".join([block, f(seed), "H", "1", "16", "15", f(ch), ".", ".", ".", f(asymH)] + common +
                              ["0", ".", "0", "7", "10", ".", ".", ".", ".", "1e-12"]))
        seedl.append(",".join([block, f(seed), "S", "1", "16", "15", f(cs), ".", ".", ".", "0"] + common +
                              ["0", ".", "0", "7", "11", f(cs - ch), f(cp - cs), f(nobs * (ch - cs)), "100",
                               "1e-12"]))
        seedl.append(",".join([block, f(seed), "S1", "1", "1", ".", f(crits["S"][0]), ".", ".", ".", "0"] + common +
                              [".", ".", ".", "7", ".", ".", ".", ".", ".", "1e-12"]))
        seedl.append(",".join([block, f(seed), "P", "1", "0", ".", f(cp), ".", ".", ".", "0"] + common +
                              [".", ".", ".", ".", "1", ".", ".", ".", ".", "1e-12"]))
        for r in range(tp.shape[0]):
            for c in range(tp.shape[1]):
                trl.append(",".join([block, f(seed), f(r + 1), f(c + 1), f(tp[r, c])]))
    for name, lines in (("cells.csv", cells), ("raw.csv", raws), ("seeds.csv", seedl)):
        with open(os.path.join(d, name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(lines) + "\n")
    if truth:
        with open(os.path.join(d, "truth.csv"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(trl) + "\n")
    if polish:
        with open(os.path.join(d, "polish.csv"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("block,seed,retcode,iterations,critDiff\n" + f"{block},0,0,1,0\n")


class TestI6(unittest.TestCase):
    def test_pf_err_matches_frozen_definition_on_e0_rows(self):
        sys.path.insert(0, FROZEN)
        try:
            import confirm_analyze as ca
        finally:
            sys.path.remove(FROZEN)
        cells = ca.load_cells(os.path.join(FROZEN, "pilot_out"), (1, 30), models=(0,))
        truth = ca.load_truth(os.path.join(FROZEN, "e0_out"), (1, 30), models=(0,))
        n = 0
        for c in cells:
            nif = ca.n_int_full(c["ng"], c["n1"], 0)
            tp = truth[(0, c["seed"])]["trueParams"]
            a, b = sa.pf_err(c["raw"], tp, nif), ca.pf_err(c["raw"], tp, nif)
            self.assertTrue(sa.same_num(a, b), (c["seed"], c["method"], c["start"], a, b))
            a, b = sa.raw_struct_err(c["raw"], tp, nif), ca.raw_struct_err(c["raw"], tp, nif, 0)
            self.assertTrue(sa.same_num(a, b))
            n += 1
        self.assertGreater(n, 1000)


class TestStats(unittest.TestCase):
    def test_clopper_pearson_thresholds(self):
        self.assertAlmostEqual(sa.cp_lower(161, 170), 0.9019, places=4)
        self.assertAlmostEqual(sa.cp_lower(160, 170), 0.8945, places=4)
        self.assertGreaterEqual(sa.cp_lower(161, 170), 0.90)
        self.assertLess(sa.cp_lower(160, 170), 0.90)

    def test_decision_table(self):
        v = sa.verdict
        self.assertEqual(v(False, True, True, True, True, True, True, True), "INVALID")
        self.assertEqual(v(True, True, True, True, True, True, True, False), "PROCEED")
        self.assertEqual(v(True, False, True, True, True, True, True, True), "REVISE")
        self.assertEqual(v(True, True, False, True, True, True, True, False), "REJECT")
        self.assertEqual(v(True, True, True, False, True, True, True, False), "REVISE")
        self.assertEqual(v(True, True, True, True, True, True, False, False), "REVISE")
        self.assertEqual(v(True, True, True, True, True, True, None, False), "PENDING_G6")
        self.assertEqual(v(True, False, False, False, False, False, None, False), "REJECT")


class TestAnalyzer(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="symtest_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def p(self, *a):
        return os.path.join(self.tmp, *a)

    def test_fixture_end_to_end(self):
        build_run(self.p("run"), "fixture", [602])
        res, code = sa.analyze("fixture", self.p("run"), self.p("out"))
        self.assertEqual(code, 0, res["integrity"])
        self.assertEqual(res["integrity"]["problems"], [])
        self.assertEqual(res["gates"]["C2"]["statistic"], 0)
        self.assertEqual(res["descriptive"]["support"]["602"], 15)
        with open(self.p("out", "seeds.csv"), encoding="utf-8") as fh:
            rows = fh.read().splitlines()
        self.assertEqual(len(rows), 5)
        s_row = dict(zip(rows[0].split(","), rows[2].split(",")))
        self.assertEqual((s_row["arm"], s_row["pfErr"], s_row["classPf"]), ("S", "0", "correct"))

    def test_selection_disagreement_is_invalid(self):
        build_run(self.p("run"), "fixture", [602])
        path = self.p("run", "seeds.csv")
        with open(path, encoding="utf-8") as fh:
            txt = fh.read()
        txt = txt.replace("fixture,602,S,1,16,", "fixture,602,S,1,15,")
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(txt)
        res, code = sa.analyze("fixture", self.p("run"), self.p("out"))
        self.assertEqual(code, 2)
        self.assertEqual(res["verdict"], "INVALID")

    def test_seed_guard_refuses_before_parsing(self):
        build_run(self.p("run"), "fixture", [602])
        with open(self.p("run", "cells.csv"), "a", encoding="utf-8", newline="\n") as fh:
            fh.write("fixture,2001,H,not-a-number,garbage\n")
        with self.assertRaises(sa.GuardError):
            sa.analyze("fixture", self.p("run"), self.p("out"))
        self.assertEqual(sa.main(["analyze", "--stage", "fixture", "--run", self.p("run"), "--out", self.p("o2")]), 3)

    def test_block_label_guard(self):
        build_run(self.p("run"), "dev", [602])
        with self.assertRaises(sa.GuardError):
            sa.analyze("fixture", self.p("run"), self.p("out"))

    def test_compare_and_byte_identical_reruns(self):
        build_run(self.p("a"), "fixture", [602])
        shutil.copytree(self.p("a"), self.p("b"))
        path = self.p("b", "seeds.csv")
        with open(path, encoding="utf-8") as fh:
            txt = fh.read()
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(txt.replace(",10,.,.,.,.,1e-12", ",99,.,.,.,.,1e-12"))
        self.assertEqual(sa.compare_dirs(self.p("a"), self.p("b")), [])
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(txt.replace("fixture,602,P,1,0,.,", "fixture,602,P,0,0,.,"))
        self.assertEqual(len(sa.compare_dirs(self.p("a"), self.p("b"))), 1)
        sa.analyze("fixture", self.p("a"), self.p("o1"))
        sa.analyze("fixture", self.p("a"), self.p("o2"))
        for name in ("results.json", "report.md", "seeds.csv"):
            with open(self.p("o1", name), "rb") as f1, open(self.p("o2", name), "rb") as f2:
                self.assertEqual(f1.read(), f2.read(), name)

    def test_interim_reads_no_truth_and_reports_no_accuracy(self):
        build_run(self.p("conf"), "conf", list(range(2001, 2021)), truth=False)
        res, code = sa.analyze("interim", self.p("conf"), self.p("out"))
        self.assertEqual(code, 0, res["integrity"])
        self.assertNotIn("gates", res)
        self.assertNotIn("descriptive", res)
        self.assertFalse(os.path.exists(self.p("out", "seeds.csv")))
        self.assertIn("projected_stage_c_s", res["runtime"])

    def test_dev_with_determinism(self):
        seeds = list(range(1001, 1031))
        build_run(self.p("dev"), "dev", seeds)
        build_run(self.p("dev_det"), "dev_det", [1001, 1030])
        res, code = sa.analyze("dev", self.p("dev"), self.p("out"), det_dir=self.p("dev_det"))
        self.assertEqual(code, 0, res["integrity"])
        self.assertEqual(res["integrity"]["determinismDiffs"], 0)
        self.assertFalse(res["devStop"]["stopRevise"])
        self.assertEqual(res["gates"]["C1"]["statistic"], 30)
        res2, code2 = sa.analyze("dev", self.p("dev"), self.p("out2"))
        self.assertEqual(code2, 2)

    def test_partial_determinism_resolve_is_invalid(self):
        seeds = list(range(1001, 1031))
        build_run(self.p("dev"), "dev", seeds)
        build_run(self.p("dev_det"), "dev_det", [1001])
        res, code = sa.analyze("dev", self.p("dev"), self.p("out"), det_dir=self.p("dev_det"))
        self.assertEqual(code, 2)
        self.assertEqual(res["verdict"], "INVALID")
        self.assertTrue(res["devStop"]["stopRevise"])
        os.remove(self.p("dev_det", "seeds.csv"))
        self.assertTrue(sa.compare_dirs(self.p("dev"), self.p("dev_det")))

    def test_outputs_and_manifest_are_never_overwritten(self):
        build_run(self.p("run"), "fixture", [602])
        sa.analyze("fixture", self.p("run"), self.p("out"))
        with self.assertRaises(sa.GuardError):
            sa.analyze("fixture", self.p("run"), self.p("out"))
        build_run(self.p("real"), "real", [0], ng=6, n1=3, nobs=32, truth=False, polish=True)
        sa.analyze_real(self.p("real"), self.p("rout"))
        with self.assertRaises(sa.GuardError):
            sa.analyze_real(self.p("real"), self.p("rout"))
        cfg = sa.load_config()
        man, dev_res = self.p("m", "SYM_FREEZE_MANIFEST.json"), self.p("dev_out", "results.json")
        with self.assertRaises(sa.GuardError):  # no stage D analysis yet
            sa.make_manifest(cfg, "test", path=man, dev_results=dev_res)
        build_run(self.p("dev"), "dev", list(range(1001, 1031)))
        build_run(self.p("dev_det"), "dev_det", [1001, 1030])
        sa.analyze("dev", self.p("dev"), self.p("dev_out"), det_dir=self.p("dev_det"))
        sa.make_manifest(cfg, "test", path=man, dev_results=dev_res)
        self.assertTrue(os.path.exists(man) and os.path.exists(man[:-5] + ".sha256"))
        with self.assertRaises(sa.GuardError):
            sa.make_manifest(cfg, "test", path=man, dev_results=dev_res)

    def test_real_g6_path(self):
        build_run(self.p("real"), "real", [0], ng=6, n1=3, nobs=32, truth=False, polish=True)
        rd, code = sa.analyze_real(self.p("real"), self.p("out"))
        self.assertEqual(code, 0, rd["integrity"])
        self.assertTrue(rd["G6"]["pass"], rd["G6"])
        self.assertEqual(len(rd["R8"]["ownPriceGammaS"]), 3)


if __name__ == "__main__":
    unittest.main()
