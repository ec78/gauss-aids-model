"""confirm_analyze.py -- confirmation analysis for CORRECTNESS_VALIDATION_PROTOCOL.md.

EXPERIMENT ONLY. Implements protocol section 4 (selection), section 5
(metrics/classes) and section 7 (A1 analysis-level integrity; gates A2 and
A4; descriptive A3, A5, A6, A7) on harness cell files. No solver.

Modes
  --mode dev      Seeds restricted to confirm_config.json "dev_seeds" (1-30).
                  The seed of every row is checked while streaming the file,
                  BEFORE any value of that row is parsed; a row outside the
                  range aborts the run. Used to test the script against E0.
  --mode confirm  Seeds restricted to "confirm_seeds" (31-200). Refuses to run
                  unless --manifest verifies (every listed file's SHA-256
                  matches), so the frozen rules cannot drift after freeze.

Usage (from this directory):
  python confirm_analyze.py --mode dev --cells pilot_out --truth e0_out --out confirm_dev_out
  python confirm_analyze.py --mode confirm --cells e2_out --truth e2_out --out confirm_e2_out --manifest FREEZE_MANIFEST.json

Inputs
  <cells>/columns.txt, <cells>/cells_m{0,1}.csv   harness output (pilot format)
  <truth>/truth_m{0,1}.csv                         truth dump (e0_truth_dump.e format)
Outputs (deterministic: sorted keys, repr floats, no timestamps)
  <out>/results.json, <out>/report.md, <out>/selected.csv
Exit codes: 0 ok; 2 integrity failure (A1 analysis checks); 3 guard refusal.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import sys

import numpy as np
from scipy.stats import beta, binomtest

HERE = os.path.dirname(os.path.abspath(__file__))


class GuardError(Exception):
    pass


class IntegrityError(Exception):
    pass


# ----------------------------------------------------------------------------
# configuration and manifest
# ----------------------------------------------------------------------------
def load_config(path=None):
    with open(path or os.path.join(HERE, "confirm_config.json"), encoding="utf-8") as fh:
        return json.load(fh)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_manifest(manifest_path):
    """Every file listed in the manifest must exist with the recorded SHA-256.
    Paths are relative to the repository root recorded in the manifest's
    location convention (manifest lives in this directory; paths are relative
    to the repository root)."""
    with open(manifest_path, encoding="utf-8") as fh:
        man = json.load(fh)
    root = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(manifest_path)), man["repo_root_rel"]))
    bad = []
    for rel, digest in sorted(man["files"].items()):
        p = os.path.join(root, rel)
        if not os.path.exists(p):
            bad.append(f"missing {rel}")
        elif sha256_file(p) != digest:
            bad.append(f"hash mismatch {rel}")
    return bad


# ----------------------------------------------------------------------------
# loading (with the seed guard applied before a row's values are parsed)
# ----------------------------------------------------------------------------
def read_meta_names(cells_dir):
    with open(os.path.join(cells_dir, "columns.txt"), encoding="utf-8") as fh:
        header = fh.readline().strip().split(",")
    idx = next(i for i, h in enumerate(header) if h.startswith("raw["))
    return header[:idx]


def _num(tok):
    return float("nan") if tok == "." else float(tok)


def load_cells(cells_dir, seed_range, models=(0, 1)):
    """Returns a list of dict cells. Aborts (GuardError) on the first row whose
    seed lies outside seed_range, before parsing the rest of that row."""
    names = read_meta_names(cells_dir)
    nmeta = len(names)
    lo, hi = seed_range
    cells = []
    for m in models:
        path = os.path.join(cells_dir, f"cells_m{m}.csv")
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as fh:
            for ln, text in enumerate(fh, 1):
                text = text.strip()
                if not text:
                    continue
                head = text.split(",", 2)
                seed = int(float(head[1]))
                if not (lo <= seed <= hi):
                    raise GuardError(f"{path}:{ln}: seed {seed} outside allowed range {lo}-{hi}; refusing to read it")
                vals = [_num(t) for t in text.split(",")]
                c = dict(zip(names, vals[:nmeta]))
                ng, n1 = int(c["ng"]), int(c["n1"])
                n = n1 + 1
                k, nbs = ng * n1, (ng + 1) * n
                if len(vals) != nmeta + k + nbs:
                    raise IntegrityError(f"{path}:{ln}: row width {len(vals)} != {nmeta + k + nbs}")
                for key in ("model", "seed", "method", "start", "ng", "n1"):
                    c[key] = int(c[key])
                if c["model"] != m:
                    raise IntegrityError(f"{path}:{ln}: model column {c['model']} != file model {m}")
                c["raw"] = np.array(vals[nmeta:nmeta + k]).reshape((ng, n1), order="F")
                c["bS"] = np.array(vals[nmeta + k:]).reshape((ng + 1, n), order="F")
                cells.append(c)
    return cells


def load_truth(truth_dir, seed_range, models=(0, 1)):
    lo, hi = seed_range
    out = {}
    for m in models:
        path = os.path.join(truth_dir, f"truth_m{m}.csv")
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as fh:
            for ln, text in enumerate(fh, 1):
                text = text.strip()
                if not text:
                    continue
                head = text.split(",", 2)
                seed = int(float(head[1]))
                if not (lo <= seed <= hi):
                    raise GuardError(f"{path}:{ln}: truth seed {seed} outside allowed range {lo}-{hi}")
                v = [_num(t) for t in text.split(",")]
                model, ng, n1 = int(v[0]), int(v[2]), int(v[3])
                n = n1 + 1
                if len(v) != 6 + (ng + 1) * n:
                    raise IntegrityError(f"{path}:{ln}: truth row width")
                out[(model, seed)] = dict(ng=ng, n1=n1, critTrue=v[4], critStone=v[5],
                                          trueParams=np.array(v[6:]).reshape((ng + 1, n), order="F"))
    return out


# ----------------------------------------------------------------------------
# layout, errors, metric
# ----------------------------------------------------------------------------
def n_int_full(ng, n1, model):
    """Intercept-row count from the stored layout: alpha (nIntFull) | gamma
    (n1) | beta | [lambda, QUAIDS] | u  =>  ng = nIntFull + n1 + 2 + model."""
    nif = ng - n1 - 2 - model
    if nif < 1:
        raise IntegrityError(f"layout: ng={ng}, n1={n1}, model={model} gives nIntFull={nif} < 1")
    return nif


def true_raw(tp, nif, n1):
    """trueParams (bS layout, ng+1 rows) -> reduced raw layout (ng rows):
    drops the n-th (recovered) gamma row, as _directTrueRaw does."""
    n = n1 + 1
    return np.vstack([tp[:nif], tp[nif:nif + n1], tp[nif + n:]])


def raw_struct_err(raw, tp, nif, model):
    """Mirror of _pilotErrs structErr: max |raw - true| over all rows except
    u, columns 1:n1 (lambda included for QUAIDS)."""
    n1 = raw.shape[1]
    tr = true_raw(tp, nif, n1)[:, :n1]
    nstruct = nif + n1 + 1 + model
    return float(np.abs(raw[:nstruct] - tr[:nstruct]).max())


def unprojected_full(raw, nif):
    """Projection-free recovery: homogeneity gives the n-th price row
    (-sum of the n1 gamma rows); adding-up gives the n-th column (-row
    sums, +1 on the constant). Same formula as pilot_bs_diag.e _diagUnproj."""
    n1 = raw.shape[1]
    gam = raw[nif:nif + n1]
    full = np.vstack([raw[:nif], gam, -gam.sum(axis=0, keepdims=True), raw[nif + n1:]])
    full = np.hstack([full, -full.sum(axis=1, keepdims=True)])
    full[0, -1] += 1.0
    return full


def pf_err(raw, tp, nif):
    """pfErr: max |recovered - true| over every row except the last (u)."""
    full = unprojected_full(raw, nif)
    if full.shape != tp.shape:
        raise IntegrityError(f"pfErr shape {full.shape} != truth {tp.shape}")
    return float(np.abs(full[:-1] - tp[:-1]).max())


# ----------------------------------------------------------------------------
# section 4: eligibility and selection
# ----------------------------------------------------------------------------
def is_eligible(c):
    """Optimizer convergence (stored conv, which is the method's own
    predicate AND exact-objective feasibility) AND finite critExact AND every
    raw coefficient finite. Convergence is an ELIGIBILITY rule, not a
    reliability signal (protocol section 6)."""
    crit = c["critExact"]
    return (c["conv"] == 1) and (crit is not None) and math.isfinite(crit) and bool(np.all(np.isfinite(c["raw"])))


def in_config(c, config):
    spec = config.get(str(c["method"]))
    return spec is not None and spec[0] <= c["start"] <= spec[1]


def select(pool):
    """Strict maximum of critExact among eligible candidates, no tolerance.
    Candidates are visited in (method, start) ascending order and replaced
    only by a strictly larger value, so exact ties go to the earliest.
    Returns (selected cell or None, eligible list)."""
    elig = sorted((c for c in pool if is_eligible(c)), key=lambda c: (c["method"], c["start"]))
    best = None
    for c in elig:
        if best is None or c["critExact"] > best["critExact"]:
            best = c
    return best, elig


def classify(err, cutoff):
    if err is None:
        return "none"
    return "correct" if err <= cutoff else "wrong"


def basin_dist(a, b):
    return float(np.abs(a[:-1, :] - b[:-1, :]).max())


# ----------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------
def cp_two_sided(k, n, alpha=0.05):
    """Exact Clopper-Pearson two-sided (1-alpha) interval for k/n."""
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def mcnemar_one_sided(b, c):
    """Exact conditional McNemar: P(B >= b | B ~ Bin(b+c, 1/2))."""
    if b + c == 0:
        return 1.0
    return float(binomtest(b, b + c, 0.5, alternative="greater").pvalue)


def quantiles(v):
    v = sorted(v)
    if not v:
        return None
    a = np.array(v)
    return {"min": float(a.min()), "median": float(np.median(a)), "p90": float(np.quantile(a, 0.9)), "max": float(a.max())}


# ----------------------------------------------------------------------------
# A1 analysis-level integrity
# ----------------------------------------------------------------------------
def integrity(cells, truth, cfg, seed_range):
    issues = []
    lo, hi = seed_range
    s0, s1 = cfg["starts"]
    methods = sorted(int(m) for m in cfg["methods"])
    keys = {}
    for c in cells:
        k = (c["model"], c["seed"], c["method"], c["start"])
        if k in keys:
            issues.append(f"duplicate key {k}")
        keys[k] = c
    models_present = sorted({c["model"] for c in cells})
    for m in models_present:
        for s in range(lo, hi + 1):
            if (m, s) not in truth:
                issues.append(f"missing truth for model {m} seed {s}")
            for meth in methods:
                for st in range(s0, s1 + 1):
                    if (m, s, meth, st) not in keys:
                        issues.append(f"missing cell {(m, s, meth, st)}")
    nchk = 0
    for c in cells:
        T = truth.get((c["model"], c["seed"]))
        if T is None:
            continue
        if (T["ng"], T["n1"]) != (c["ng"], c["n1"]):
            issues.append(f"layout mismatch cell/truth {(c['model'], c['seed'])}")
            continue
        nif = n_int_full(c["ng"], c["n1"], c["model"])
        tp = T["trueParams"]
        checks = [("rawStructErr", raw_struct_err(c["raw"], tp, nif, c["model"])), ("critTrue", T["critTrue"])]
        if c["bsOk"] == 1:
            checks.append(("bsErr", float(np.abs(c["bS"] - tp).max())))
        for name, val in checks:
            st = c[name]
            same = (val == st) or (math.isnan(val) and math.isnan(st))
            nchk += 1
            if not same:
                issues.append(f"truth gate {name} {(c['model'], c['seed'], c['method'], c['start'])}: stored {st!r} regenerated {val!r}")
    return issues, nchk


# ----------------------------------------------------------------------------
# analysis
# ----------------------------------------------------------------------------
def analyze(cells, truth, cfg, seed_range):
    issues, nchk = integrity(cells, truth, cfg, seed_range)
    if issues:
        raise IntegrityError(f"{len(issues)} integrity issue(s); first: " + "; ".join(issues[:5]))

    for c in cells:
        T = truth[(c["model"], c["seed"])]
        nif = n_int_full(c["ng"], c["n1"], c["model"])
        c["pfErr"] = pf_err(c["raw"], T["trueParams"], nif) if np.all(np.isfinite(c["raw"])) else None

    by_seed = {}
    for c in cells:
        by_seed.setdefault((c["model"], c["seed"]), []).append(c)

    cut = cfg["primary_cutoff"]
    sel_rows = []
    per = {}
    for (m, s) in sorted(by_seed):
        group = by_seed[(m, s)]
        for name in sorted(cfg["configurations"]):
            conf = cfg["configurations"][name]
            pool = [c for c in group if in_config(c, conf)]
            best, elig = select(pool)
            row = {"model": m, "seed": s, "config": name, "nPool": len(pool), "nEligible": len(elig)}
            if best is None:
                row.update(selMethod=None, selStart=None, pfErr=None, rawStructErr=None, bsErr=None,
                           critExact=None, classPf="none", classRaw="none", classBs="none",
                           exactTies=0, marginOtherBasin=None, nearTie=False)
            else:
                ties = sum(1 for c in elig if c["critExact"] == best["critExact"]) - 1
                others = [c["critExact"] for c in elig if basin_dist(c["raw"], best["raw"]) > cfg["basin_tol"]]
                margin = (best["critExact"] - max(others)) if others else None
                bs = best["bsErr"] if best["bsOk"] == 1 and math.isfinite(best["bsErr"]) else None
                row.update(selMethod=best["method"], selStart=best["start"], pfErr=best["pfErr"],
                           rawStructErr=best["rawStructErr"], bsErr=bs, critExact=best["critExact"],
                           classPf=classify(best["pfErr"], cut), classRaw=classify(best["rawStructErr"], cut),
                           classBs=classify(bs, cut) if bs is not None else "none",
                           exactTies=ties, marginOtherBasin=margin,
                           nearTie=(margin is not None and margin < cfg["near_tie_eps"]))
                if best["pfErr"] is None:
                    raise IntegrityError(f"selected endpoint with non-finite pfErr {(m, s, name)}")
            row["wallTotal"] = float(sum(c["wallTime"] for c in pool))
            sel_rows.append(row)
            per.setdefault((m, name), []).append(row)

    results = {"seedRange": list(seed_range), "integrity": {"truthValuesReproduced": nchk, "issues": 0},
               "models": {}}
    n_seeds = seed_range[1] - seed_range[0] + 1
    for m in sorted({k[0] for k in per}):
        mname = cfg["models"][str(m)]
        mres = {"configurations": {}}
        for name in sorted(cfg["configurations"]):
            rows = per[(m, name)]
            n = len(rows)
            cnt = {k: sum(1 for r in rows if r["classPf"] == k) for k in ("none", "wrong", "correct")}
            k_ok = cnt["correct"]
            W = cnt["wrong"]
            lo, hi = cp_two_sided(k_ok, n, cfg["a2"]["alpha_two_sided"])
            wlo, whi = cp_two_sided(W, n, cfg["a4"]["alpha_two_sided"])
            def thr(col):
                return [sum(1 for r in rows if r[col] is not None and r[col] <= t) for t in cfg["sensitivity_cutoffs"]]
            mres["configurations"][name] = {
                "n": n,
                "classesPf": cnt,
                "classesRaw": {k: sum(1 for r in rows if r["classRaw"] == k) for k in ("none", "wrong", "correct")},
                "classesBs": {k: sum(1 for r in rows if r["classBs"] == k) for k in ("none", "wrong", "correct")},
                "thresholdCountsPf": thr("pfErr"),
                "thresholdCountsRaw": thr("rawStructErr"),
                "A2": {"k": k_ok, "cpLower": lo, "cpUpper": hi},
                "A4": {"W": W, "cpUpper": whi, "wrongSeeds": [r["seed"] for r in rows if r["classPf"] == "wrong"]},
                "noResultSeeds": [r["seed"] for r in rows if r["classPf"] == "none"],
                "borderline": [[r["seed"], r["pfErr"]] for r in rows if r["pfErr"] is not None
                               and cfg["borderline_band"][0] < r["pfErr"] <= cfg["borderline_band"][1]],
                "nearTieSeeds": [r["seed"] for r in rows if r["nearTie"]],
                "maxExactTies": max(r["exactTies"] for r in rows),
                "A7wallPerSeed": quantiles([r["wallTotal"] for r in rows]),
            }
        prim = cfg["primary"][str(m)]
        pc = mres["configurations"][prim]
        target = cfg["a2"]["target_lower_limit"][str(m)]
        mres["primary"] = prim
        mres["gates"] = {
            "A1_analysis": True,
            "A2": {"pass": pc["A2"]["cpLower"] >= target, "target": target, "k": pc["A2"]["k"],
                   "n": pc["n"], "cpLower": pc["A2"]["cpLower"]},
            "A4": {"pass": pc["A4"]["W"] <= cfg["a4"]["max_returned_wrong"], "W": pc["A4"]["W"],
                   "max": cfg["a4"]["max_returned_wrong"], "cpUpper": pc["A4"]["cpUpper"]},
        }
        mres["gatesPassAnalysisLevel"] = mres["gates"]["A2"]["pass"] and mres["gates"]["A4"]["pass"]
        # descriptive A3 / A5 (reference arm is NOT blind: protocol section 9)
        ref = {r["seed"]: r for r in per[(m, cfg["reference"])]}
        pri = {r["seed"]: r for r in per[(m, prim)]}
        b = [s for s in sorted(pri) if pri[s]["classPf"] == "correct" and ref[s]["classPf"] != "correct"]
        cc = [s for s in sorted(pri) if ref[s]["classPf"] == "correct" and pri[s]["classPf"] != "correct"]
        mres["A3_descriptive"] = {"gating": False, "b": len(b), "c": len(cc), "ties": len(pri) - len(b) - len(cc),
                                  "oneSidedExactP": mcnemar_one_sided(len(b), len(cc)),
                                  "improvedSeeds": b, "regressedSeeds": cc}
        mres["A5_regressionsToExplain"] = cc
        # A6 robustness (descriptive): primary classes under raw and at .5 / 2
        mres["A6_descriptive"] = {
            "rawCorrect": pc["classesRaw"]["correct"],
            "pfCorrectAt0.5": pc["thresholdCountsPf"][cfg["sensitivity_cutoffs"].index(0.5)],
            "pfCorrectAt2": pc["thresholdCountsPf"][cfg["sensitivity_cutoffs"].index(2.0)],
        }
        a7 = pc["A7wallPerSeed"]
        slow = sorted([c["seed"], c["method"], c["start"], c["wallTime"]] for c in cells
                      if c["model"] == m and c["wallTime"] > cfg["a7"]["slow_cell_s"])
        mres["A7_descriptive"] = {"wallPerSeed": a7, "medianTarget": cfg["a7"]["median_wall_target_s"][str(m)],
                                  "slowCells": slow}
        mres["name"] = mname
        if n_seeds != pc["n"]:
            raise IntegrityError(f"model {m}: {pc['n']} seeds analysed, expected {n_seeds}")
        results["models"][str(m)] = mres
    return results, sel_rows


# ----------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------
def _clean(o):
    if isinstance(o, float):
        return None if math.isnan(o) else o
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return _clean(float(o))
    return o


def write_outputs(results, sel_rows, cfg, out_dir, mode):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(_clean(results), fh, indent=1, sort_keys=True)
        fh.write("\n")
    cols = ["model", "seed", "config", "nPool", "nEligible", "selMethod", "selStart", "critExact", "pfErr",
            "rawStructErr", "bsErr", "classPf", "classRaw", "classBs", "exactTies", "marginOtherBasin",
            "nearTie", "wallTotal"]
    with open(os.path.join(out_dir, "selected.csv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for r in sel_rows:
            w.writerow(["" if r[c] is None else (repr(r[c]) if isinstance(r[c], float) else r[c]) for c in cols])
    L = [f"# Confirmation analysis ({mode} mode, seeds {results['seedRange'][0]}-{results['seedRange'][1]})", "",
         "Generated by `confirm_analyze.py`. Gates: A2 and A4 (plus A1). A3 is DESCRIPTIVE ONLY "
         "(reference arm not blind; protocol section 9).", "",
         ("Gate verdicts are BINDING (confirm mode, n = 170)." if results.get("gatesBinding") else
          "Gate verdicts are INFORMATIONAL ONLY in dev mode: the A2 thresholds are sized for n = 170, and "
          "with n = 30 even 30/30 cannot reach the AIDS target."), "",
         f"A1 analysis-level integrity: PASS ({results['integrity']['truthValuesReproduced']} truth-gate values reproduced exactly).", ""]
    for m, mres in sorted(results["models"].items()):
        g = mres["gates"]
        L.append(f"## {mres['name']} (primary `{mres['primary']}`)")
        L.append("")
        L.append(f"- A2: k = {g['A2']['k']}/{g['A2']['n']}, CP two-sided 95% lower = {g['A2']['cpLower']:.4f}, "
                 f"target {g['A2']['target']} -> {'PASS' if g['A2']['pass'] else 'FAIL'}")
        L.append(f"- A4: W = {g['A4']['W']} (max {g['A4']['max']}), CP upper = {g['A4']['cpUpper']:.4f} -> "
                 f"{'PASS' if g['A4']['pass'] else 'FAIL'}")
        a3 = mres["A3_descriptive"]
        L.append(f"- A3 (descriptive, not a gate): b = {a3['b']}, c = {a3['c']}, ties = {a3['ties']}, "
                 f"one-sided exact p = {a3['oneSidedExactP']:.3g}; improved {a3['improvedSeeds']}; regressed {a3['regressedSeeds']}")
        L.append(f"- A5 regressions needing a written mechanism: {mres['A5_regressionsToExplain'] or 'none'}")
        L.append(f"- A6 (descriptive): raw-correct {mres['A6_descriptive']['rawCorrect']}, pfErr<=0.5 "
                 f"{mres['A6_descriptive']['pfCorrectAt0.5']}, pfErr<=2 {mres['A6_descriptive']['pfCorrectAt2']}")
        L.append("")
        L.append("| config | pfErr none/wrong/correct | raw | bS | W | no-result seeds | borderline | near-tie seeds |")
        L.append("|---|---|---|---|---|---|---|---|")
        for name, cr in sorted(mres["configurations"].items()):
            f = lambda d: f"{d['none']}/{d['wrong']}/{d['correct']}"
            L.append(f"| {name} | {f(cr['classesPf'])} | {f(cr['classesRaw'])} | {f(cr['classesBs'])} | "
                     f"{cr['A4']['W']} | {cr['noResultSeeds']} | {[s for s, _ in cr['borderline']]} | {cr['nearTieSeeds']} |")
        L.append("")
    with open(os.path.join(out_dir, "report.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")


def run(mode, cells_dir, truth_dir, out_dir, manifest=None, cfg_path=None):
    cfg = load_config(cfg_path)
    if mode == "dev":
        seed_range = tuple(cfg["dev_seeds"])
    elif mode == "confirm":
        seed_range = tuple(cfg["confirm_seeds"])
        if not manifest:
            raise GuardError("confirm mode requires --manifest")
        bad = verify_manifest(manifest)
        if bad:
            raise GuardError("manifest verification failed: " + "; ".join(bad[:5]))
    else:
        raise GuardError(f"unknown mode {mode}")
    cells = load_cells(cells_dir, seed_range)
    truth = load_truth(truth_dir, seed_range)
    if not cells:
        raise IntegrityError(f"no cells found in {cells_dir}")
    present = sorted({c["model"] for c in cells})
    if mode == "confirm" and present != sorted(int(m) for m in cfg["models"]):
        raise IntegrityError(f"confirm mode requires every model; present {present}")
    results, sel_rows = analyze(cells, truth, cfg, seed_range)
    results["mode"] = mode
    # Gate verdicts are binding only in confirm mode (n = 170). In dev mode
    # (n = 30) they are informational: even 30/30 has a two-sided 95% CP
    # lower limit of 0.884, below the AIDS target, by construction.
    results["gatesBinding"] = (mode == "confirm")
    write_outputs(results, sel_rows, cfg, out_dir, mode)
    return results, sel_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["dev", "confirm"])
    ap.add_argument("--cells", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--manifest")
    ap.add_argument("--config")
    a = ap.parse_args()
    try:
        results, _ = run(a.mode, a.cells, a.truth, a.out, a.manifest, a.config)
    except GuardError as e:
        print("GUARD REFUSAL:", e)
        sys.exit(3)
    except IntegrityError as e:
        print("INTEGRITY FAILURE:", e)
        sys.exit(2)
    tag = "" if results["gatesBinding"] else " [informational: dev mode, n=30]"
    for m, mres in sorted(results["models"].items()):
        g = mres["gates"]
        print(f"{mres['name']}{tag}: A2 {'PASS' if g['A2']['pass'] else 'FAIL'} (k={g['A2']['k']}/{g['A2']['n']}), "
              f"A4 {'PASS' if g['A4']['pass'] else 'FAIL'} (W={g['A4']['W']})")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
