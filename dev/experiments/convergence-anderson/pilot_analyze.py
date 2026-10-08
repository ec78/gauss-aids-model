"""pilot_analyze.py -- analysis for the matched AIDS/QUAIDS pilot.

EXPERIMENT ONLY. Reads the per-cell CSVs written by
pilot_harness_core.src and applies the frozen rules in PILOT_PROTOCOL.md.
Selection, oracle, basin deduplication and every aggregate are computed
here, from disk, so they can be re-run without re-solving anything.

Usage (from this directory):
    python pilot_analyze.py --selftest
    python pilot_analyze.py --validate pilot_validation
    python pilot_analyze.py --report pilot_out [--seeds 1-30]

--selftest  runs V5 (selection direction / eligibility unit tests) and the
            basin-dedup unit test only.
--validate  runs V3b (benchmark reproduction), V3c (smoke-gate
            reproduction), V4 (exact-objective consistency), V5, V7 (keys),
            and prints PASS/FAIL lines plus a final summary.
--report    writes <dir>/pilot_report.md and <dir>/pilot_seed_table.csv.
"""

import argparse
import math
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = {1: "Anderson", 2: "cold LM", 3: "continuation"}
MODELS = {0: "AIDS", 1: "QUAIDS"}
THRESHOLDS = [0.5, 1.0, 2.0, 5.0]
PRIMARY_CUT = 1.0
BASIN_TOL = 0.05          # raw structural (excl. u), max abs diff
BENCH_DISTINCT_TOL = 1e-4  # the benchmark's own vec(bS) dedup
IDENT_TOL = 1e-6
NUM_STARTS = 16


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def load_cells(outdir):
    with open(os.path.join(outdir, "columns.txt")) as fh:
        header = fh.readline().strip().split(",")
    meta = header[:next(i for i, h in enumerate(header) if h.startswith("raw["))]
    nmeta = len(meta)
    frames = []
    for model in (0, 1):
        path = os.path.join(outdir, f"cells_m{model}.csv")
        if not os.path.exists(path):
            continue
        rows = []
        with open(path) as fh:
            for ln, text in enumerate(fh, 1):
                text = text.strip()
                if not text:
                    continue
                vals = [float("nan") if t == "." else float(t) for t in text.split(",")]
                rows.append((ln, vals))
        recs = []
        for ln, vals in rows:
            m = dict(zip(meta, vals[:nmeta]))
            ng, n1 = int(m["ng"]), int(m["n1"])
            n = n1 + 1
            k = ng * n1
            nbs = (ng + 1) * n
            if len(vals) != nmeta + k + nbs:
                raise ValueError(f"{path}:{ln}: row width {len(vals)} != {nmeta + k + nbs}")
            m["raw"] = np.array(vals[nmeta:nmeta + k]).reshape((ng, n1), order="F")
            m["bS"] = np.array(vals[nmeta + k:]).reshape((ng + 1, n), order="F")
            m["_line"] = ln
            recs.append(m)
        frames.append(pd.DataFrame(recs))
    if not frames:
        raise SystemExit(f"no cells found in {outdir}")
    df = pd.concat(frames, ignore_index=True)
    for c in ("model", "seed", "method", "start", "status", "ng", "n1"):
        df[c] = df[c].astype(int)
    df["conv"] = df["conv"] == 1
    return df


# --------------------------------------------------------------------------
# frozen rules
# --------------------------------------------------------------------------
def select_operational(cands):
    """cands: iterable of (start, conv, critExact). Returns the start with
    the LARGEST critExact among genuinely converged candidates (ties ->
    lowest start), or None if none converged. critExact = -ln det(Sigma):
    larger is better."""
    best = None
    for start, conv, crit in sorted(cands, key=lambda c: c[0]):
        if not conv or crit is None or (isinstance(crit, float) and math.isnan(crit)):
            continue
        if best is None or crit > best[1]:
            best = (start, crit)
    return None if best is None else best[0]


def select_oracle(cands):
    """cands: (start, conv, err). Smallest error among converged; ties ->
    lowest start. DIAGNOSTIC ONLY."""
    best = None
    for start, conv, err in sorted(cands, key=lambda c: c[0]):
        if not conv or err is None or math.isnan(err):
            continue
        if best is None or err < best[1]:
            best = (start, err)
    return None if best is None else best[0]


def struct_rows(raw):
    """All reduced rows except the last (u) row."""
    return raw[:-1, :]


def dedup(endpoints, tol, key=lambda e: e):
    """Leader clustering in the given order: a new endpoint opens a new
    basin if its max-abs distance to every representative exceeds tol.
    Returns (representatives' indices, label per endpoint)."""
    reps, labels = [], []
    for i, e in enumerate(endpoints):
        v = key(e)
        lab = None
        for j, r in enumerate(reps):
            if np.max(np.abs(key(endpoints[r]) - v)) <= tol:
                lab = j
                break
        if lab is None:
            reps.append(i)
            lab = len(reps) - 1
        labels.append(lab)
    return reps, labels


def classify(err, cut=PRIMARY_CUT):
    if err is None or (isinstance(err, float) and math.isnan(err)):
        return "never"
    return "correct" if err <= cut else "wrong"


# --------------------------------------------------------------------------
# V5 selftest
# --------------------------------------------------------------------------
def selftest(report):
    ok = True

    def chk(name, cond):
        nonlocal ok
        report.append(("PASS" if cond else "FAIL") + " V5 " + name)
        ok = ok and cond

    # larger critExact (= -ln det) wins
    c = [(1, True, -7.84), (2, True, -5.47), (3, True, -6.0)]
    chk("largest critExact selected (-5.47 over -6.0, -7.84)", select_operational(c) == 2)
    # sign-reversal guard: selecting by min would pick start 1; must not
    chk("sign reversal would pick the worst basin; rule does not", select_operational(c) != 1)
    # non-converged better objective excluded
    c = [(1, False, -1.0), (2, True, -5.0), (3, True, -9.0)]
    chk("non-converged better objective excluded", select_operational(c) == 2)
    # ties -> lowest start
    c = [(4, True, -5.0), (2, True, -5.0), (9, True, -6.0)]
    chk("tie -> lowest start", select_operational(c) == 2)
    # none converged
    c = [(1, False, -5.0), (2, False, -4.0)]
    chk("no converged start -> never (None)", select_operational(c) is None)
    # NaN objective on a converged flag is never selected
    c = [(1, True, float("nan")), (2, True, -8.0)]
    chk("NaN objective skipped", select_operational(c) == 2)
    # oracle picks smallest error among converged only
    c = [(1, False, 0.01), (2, True, 0.9), (3, True, 0.3)]
    chk("oracle smallest error among converged", select_oracle(c) == 3)
    # basin dedup
    a = np.zeros((3, 2))
    eps = [a, a + 0.04, a + 0.06, a + 1.0, a + 1.03]
    reps, labels = dedup(eps, BASIN_TOL)
    chk("basin dedup leader clustering (0,.04 | .06 | 1,1.03)", labels == [0, 0, 1, 2, 2])
    return ok


# --------------------------------------------------------------------------
# per (model, seed, method) portfolio summary
# --------------------------------------------------------------------------
def summarize_portfolio(g):
    g = g.sort_values("start")
    conv = g[g["conv"]]
    out = {
        "nAttempted": len(g),
        "nConv": int(g["conv"].sum()),
    }
    sel = select_operational(zip(g["start"], g["conv"], g["critExact"]))
    orb = select_oracle(zip(g["start"], g["conv"], g["bsErr"]))
    orr = select_oracle(zip(g["start"], g["conv"], g["rawStructErr"]))
    row = lambda s: g[g["start"] == s].iloc[0]
    if sel is None:
        out.update(selStart=np.nan, selBsErr=np.nan, selRawErr=np.nan, selCrit=np.nan,
                   selBsErrExclU=np.nan, selUErr=np.nan)
    else:
        r = row(sel)
        out.update(selStart=sel, selBsErr=r["bsErr"], selRawErr=r["rawStructErr"],
                   selCrit=r["critExact"], selBsErrExclU=r["bsErrExclU"], selUErr=r["rawUErr"])
    out["oraBsStart"] = np.nan if orb is None else orb
    out["oraBsErr"] = np.nan if orb is None else row(orb)["bsErr"]
    out["oraRawStart"] = np.nan if orr is None else orr
    out["oraRawErr"] = np.nan if orr is None else row(orr)["rawStructErr"]
    # Stone-only
    s1 = g[g["start"] == 1]
    if len(s1):
        r1 = s1.iloc[0]
        out.update(stoneConv=bool(r1["conv"]),
                   stoneBsErr=r1["bsErr"] if r1["conv"] else np.nan,
                   stoneRawErr=r1["rawStructErr"] if r1["conv"] else np.nan,
                   stoneWall=r1["wallTime"], stoneCrit=r1["critExact"])
    # basins
    eps = list(conv["raw"])
    reps, labels = dedup(eps, BASIN_TOL, key=struct_rows)
    out["nBasins"] = len(reps)
    reps2, _ = dedup(list(conv["bS"]), BENCH_DISTINCT_TOL, key=lambda b: b.ravel(order="F"))
    out["nDistinctBench"] = len(reps2)
    # basin of the selected endpoint, best-basin crit gap
    out["critTrue"] = g["critTrue"].iloc[0]
    # counters and costs over the whole portfolio
    out["wallTotal"] = g["wallTime"].sum()
    out["iterTotal"] = g["iterations"].sum()
    out["evalsTotal"] = g["approxEvals"].sum()
    for c in ("covInvFail", "stepInvFail", "infeasTrial", "rhoReject",
              "cFailRc3", "cFailRc5", "cFailRc2", "cFailRc4", "cStepReductions",
              "cMinStepExhausted", "cAnyStageFailed"):
        out[c] = g[c].sum(skipna=True)
    out["statusCounts"] = dict(g["status"].value_counts().sort_index())
    out["bsFail"] = int((g["bsOk"] != 1).sum())
    out["symInvalid"] = int((g["bsSymValid"] != 1).sum())
    return out


def build_seed_table(df):
    recs = []
    for (model, seed, method), g in df.groupby(["model", "seed", "method"]):
        s = summarize_portfolio(g)
        s.update(model=model, seed=seed, method=method)
        recs.append(s)
    t = pd.DataFrame(recs)
    t["classBs"] = [classify(e) for e in t["selBsErr"]]
    t["classRaw"] = [classify(e) for e in t["selRawErr"]]
    t["oraClassBs"] = [classify(e) for e in t["oraBsErr"]]
    t["stoneClassBs"] = [classify(e) for e in t["stoneBsErr"]]
    t["stoneClassRaw"] = [classify(e) for e in t["stoneRawErr"]]
    return t


# --------------------------------------------------------------------------
# validation (V3b, V3c, V4, V7)
# --------------------------------------------------------------------------
def parse_benchmark():
    path = os.path.join(HERE, "anderson_full_multistart_benchmark_output.txt")
    res = {}
    num = r"([-+0-9.eE]+)"
    pat = re.compile(r"^seed\s+" + num + r"\s+model(Iterated AIDS \(linear\)|QUAIDS)(.*)$")
    with open(path) as fh:
        for text in fh:
            m = pat.match(text.strip())
            if not m:
                continue
            seed = int(float(m.group(1)))
            model = 0 if m.group(2).startswith("Iterated") else 1
            rest = m.group(3)
            if "no safeguarded start converged" in rest:
                res[(model, seed)] = None
                continue
            fields = dict(re.findall(r"([A-Za-z]+)\s+" + num, rest))
            res[(model, seed)] = {k: float(v) for k, v in fields.items()}
    return res


SMOKE_REF = {
    # (model, seed, method): (critExact, rawStructErr, iterations, rejected)
    (0, 6, 2): (-5.4987359, 0.058349684, 4, 0),
    (0, 6, 3): (-5.4987232, 0.058367908, 33, 0),
    (1, 9, 2): (-11.344578, 4.0194095, 44, 8),
    (1, 9, 3): (-5.5772838, 0.081498821, 82, 8),
}


def validate(outdir):
    report = []
    ok_all = True

    def chk(name, cond, detail=""):
        nonlocal ok_all
        report.append(("PASS " if cond else "FAIL ") + name + (("  " + detail) if detail else ""))
        ok_all = ok_all and bool(cond)

    df = load_cells(outdir)
    # V7 keys
    keys = df[["model", "seed", "method", "start"]]
    chk("V7 keys unique", not keys.duplicated().any(), f"{len(df)} rows")
    seeds = sorted(df["seed"].unique())
    expected = 2 * len(seeds) * 3 * NUM_STARTS
    chk("V7 complete grid", len(df) == expected, f"{len(df)} of {expected}")
    for c in ("covInvFail", "stepInvFail", "infeasTrial", "rhoReject", "wallTime",
              "iterations", "status"):
        sub = df[df["method"].isin([2, 3])] if c in ("covInvFail", "stepInvFail",
                                                     "infeasTrial", "rhoReject") else df
        chk(f"V7 counter column '{c}' present and non-missing", sub[c].notna().all())
    cont = df[df["method"] == 3]
    for c in ("cStage0Ok", "cStepReductions", "cMinStepExhausted", "cFailRc2",
              "cFailRc3", "cFailRc4", "cFailRc5", "cFinalGnorm"):
        chk(f"V7 continuation column '{c}' present and non-missing", cont[c].notna().all())

    # V4 exact objective consistency
    lm = df[(df["method"] == 2) & df["critMethod"].notna()]
    d = (lm["critExact"] - lm["critMethod"]).abs().max()
    chk("V4 cold LM: critExact == returned objective", d == 0, f"max |diff| {d:.3g}")
    cr = df[(df["method"] == 3) & (df["finalHomotopy"] == 1) & df["critMethod"].notna()]
    d = (cr["critExact"] - cr["critMethod"]).abs().max()
    chk("V4 continuation (reached 1): critExact == returned objective", d == 0, f"max |diff| {d:.3g}")
    an = df[(df["method"] == 1) & df["conv"]]
    d = (an["critExact"] - an["critMethod"]).abs()
    chk("V4 Anderson converged: |critExact - qOut.homogCrit| small (<1e-2)",
        d.max() < 1e-2, f"max {d.max():.3g}, median {d.median():.3g}")
    d = an["bsHelperDiff"]
    report.append(f"INFO V1 Anderson converged endpoints: max|bS(b,b) - qOut.bS| median {d.median():.3g}, "
                  f"90% {d.quantile(.9):.3g}, max {d.max():.3g}; "
                  f"max |bsHelperErr - bsErr| {(an['bsHelperErr'] - an['bsErr']).abs().max():.3g}")

    # V3b benchmark reproduction for Anderson
    bench = parse_benchmark()
    for (model, seed), g in df[df["method"] == 1].groupby(["model", "seed"]):
        ref = bench.get((model, seed), "missing")
        s = summarize_portfolio(g)
        # the benchmark selected by qOut.homogCrit (critMethod), not critExact
        selBench = select_operational(zip(g["start"], g["conv"], g["critMethod"]))
        tag = f"V3b {MODELS[model]} seed {seed}"
        if ref == "missing":
            chk(tag + " benchmark line found", False)
            continue
        if ref is None:
            chk(tag + " no converged start (benchmark: none)", s["nConv"] == 0, f"pilot nConv {s['nConv']}")
            continue
        rowB = g[g["start"] == selBench].iloc[0]
        chk(tag + " converged starts", s["nConv"] == int(ref["convergedStarts"]),
            f"{s['nConv']} vs {int(ref['convergedStarts'])}")
        chk(tag + " distinct endpoints (bS, 1e-4)", s["nDistinctBench"] == int(ref["distinct"]),
            f"{s['nDistinctBench']} vs {int(ref['distinct'])}")
        chk(tag + " selected start (by homogCrit)", selBench == int(ref["selectedStart"]),
            f"{selBench} vs {int(ref['selectedStart'])}")
        chk(tag + " selected recErr", abs(rowB["bsErr"] - ref["selectedRecErr"]) <= 1e-7 * max(1.0, abs(ref["selectedRecErr"])),
            f"{rowB['bsErr']:.8g} vs {ref['selectedRecErr']:.8g}")
        chk(tag + " oracle recErr", abs(s["oraBsErr"] - ref["oracleRecErr"]) <= 1e-7 * max(1.0, abs(ref["oracleRecErr"])),
            f"{s['oraBsErr']:.8g} vs {ref['oracleRecErr']:.8g}")
        chk(tag + " truth criterion", abs(g["critTrue"].iloc[0] - ref["trueCrit"]) < 1e-3,
            f"{g['critTrue'].iloc[0]:.8g} vs {ref['trueCrit']:.8g} (direct objective clips exp; replica tol)")
        report.append(f"INFO {tag}: exact-objective selection picks start {s['selStart']} "
                      f"(benchmark homogCrit selection {selBench})")

    # V3c smoke-gate reproduction
    for (model, seed, method), (crit, err, it, rej) in SMOKE_REF.items():
        r = df[(df["model"] == model) & (df["seed"] == seed) & (df["method"] == method) & (df["start"] == 1)]
        tag = f"V3c smoke gate {MODELS[model]} seed {seed} {METHODS[method]} start 1"
        if not len(r):
            report.append(f"SKIP {tag} (cell not in this run)")
            continue
        r = r.iloc[0]
        # references are printed to 8 significant digits: relative tol 1e-7
        chk(tag, abs(r["critMethod"] - crit) <= 1e-7 * abs(crit)
            and abs(r["rawStructErr"] - err) <= 1e-7 * abs(err)
            and int(r["iterations"]) == it and int(r["rejectedSteps"]) == rej,
            f"crit {r['critMethod']:.8g} vs {crit}; structErr {r['rawStructErr']:.8g} vs {err}; "
            f"iter {int(r['iterations'])} vs {it}; rejected {int(r['rejectedSteps'])} vs {rej}")

    ok5 = selftest(report)
    ok_all = ok_all and ok5
    nfail = sum(1 for x in report if x.startswith("FAIL"))
    npass = sum(1 for x in report if x.startswith("PASS"))
    for x in report:
        print(x)
    print()
    if ok_all:
        print(f"VALIDATION: ALL {npass} CHECKS PASSED")
    else:
        print(f"VALIDATION: {nfail} CHECKS FAILED (of {npass + nfail})")
    return ok_all


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "-"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{x:.{nd}g}" if abs(x) < 1e4 else f"{x:.3e}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def quant(s):
    s = pd.Series(s).dropna()
    if not len(s):
        return "-"
    q = s.quantile([0, .25, .5, .75, .9, 1]).values
    return " / ".join(fmt(v) for v in q)


def write_report(outdir, seed_filter=None):
    df = load_cells(outdir)
    if seed_filter:
        df = df[df["seed"].isin(seed_filter)]
    t = build_seed_table(df)
    seeds = sorted(df["seed"].unique())
    models = sorted(df["model"].unique())
    L = []
    L.append(f"# Matched pilot report ({outdir})")
    L.append("")
    L.append(f"Generated by `pilot_analyze.py --report {outdir}` from the on-disk cells. "
             f"Seeds: {seeds[0]}-{seeds[-1]} ({len(seeds)} seeds); cells: {len(df)}. "
             "Primary class = selected bsErr <= 1 (benchmark-compatible); secondary = "
             "selected rawStructErr <= 1. Oracle columns are DIAGNOSTIC ONLY.")
    L.append("")
    keys = df[["model", "seed", "method", "start"]]
    L.append(f"Integrity: duplicate keys {int(keys.duplicated().sum())}; expected cells "
             f"{len(models)*len(seeds)*3*NUM_STARTS}, present {len(df)}.")
    L.append("")

    # ---- 1. counts by model x method
    L.append("## 1. Classification counts (of %d seeds)" % len(seeds))
    L.append("")
    rows = []
    for model in models:
        for method in (1, 2, 3):
            s = t[(t.model == model) & (t.method == method)]
            def cnt(col):
                v = s[col].value_counts()
                return f"{v.get('never',0)} / {v.get('wrong',0)} / {v.get('correct',0)}"
            rows.append([MODELS[model], METHODS[method], cnt("stoneClassBs"), cnt("classBs"),
                         cnt("oraClassBs"), cnt("stoneClassRaw"), cnt("classRaw"),
                         int(s["nConv"].sum()), f"{len(s)*NUM_STARTS}"])
    L.append(md_table(["model", "method", "Stone-only bS n/w/c", "16-start selected bS n/w/c",
                       "oracle bS n/w/c (diag.)", "Stone-only raw n/w/c", "selected raw n/w/c",
                       "converged starts", "attempted"], rows))
    L.append("")

    # ---- 2. threshold sensitivity
    L.append("## 2. Threshold sensitivity (seeds with selected/oracle error <= cutoff)")
    L.append("")
    rows = []
    for model in models:
        for method in (1, 2, 3):
            s = t[(t.model == model) & (t.method == method)]
            def thr(col):
                return " / ".join(str(int((s[col] <= c).sum())) for c in THRESHOLDS)
            rows.append([MODELS[model], METHODS[method], thr("stoneBsErr"), thr("selBsErr"),
                         thr("oraBsErr"), thr("stoneRawErr"), thr("selRawErr"), thr("oraRawErr"),
                         thr("selBsErrExclU")])
    L.append(md_table(["model", "method", "Stone bS .5/1/2/5", "selected bS .5/1/2/5",
                       "oracle bS (diag.)", "Stone raw .5/1/2/5", "selected raw .5/1/2/5",
                       "oracle raw (diag.)", "selected bS excl. u row .5/1/2/5 (diag.)"], rows))
    L.append("")

    # ---- 3. paired transitions
    L.append("## 3. Paired seed-level transitions (primary bS class, 16-start selected)")
    L.append("")
    piv = t.pivot_table(index=["model", "seed"], columns="method", values="classBs", aggfunc="first")
    rows = []
    for model in models:
        p = piv.loc[model]
        def trans(a, b):
            out = {}
            for _, r in p.iterrows():
                k = f"{r[a]}->{r[b]}"
                out[k] = out.get(k, 0) + 1
            return out
        for a, b in ((1, 3), (2, 3), (1, 2)):
            tr = trans(a, b)
            imp = sum(v for k, v in tr.items() if k.split("->")[0] != "correct" and k.endswith("correct"))
            reg = sum(v for k, v in tr.items() if k.startswith("correct") and not k.endswith("correct"))
            same = sum(v for k, v in tr.items() if k.split("->")[0] == k.split("->")[1])
            other = len(p) - imp - reg - same
            detail = ", ".join(f"{k}: {v}" for k, v in sorted(tr.items()))
            rows.append([MODELS[model], f"{METHODS[a]} -> {METHODS[b]}", imp, reg, same, other, detail])
    L.append(md_table(["model", "pair", "improved to correct", "regressed from correct",
                       "unchanged class", "other change", "all transitions"], rows))
    L.append("")
    # explicit seed lists
    for model in models:
        p = piv.loc[model]
        a2c = [s for s, r in p.iterrows() if r[1] != "correct" and r[3] == "correct"]
        c2a = [s for s, r in p.iterrows() if r[1] == "correct" and r[3] != "correct"]
        l2c = [s for s, r in p.iterrows() if r[2] != "correct" and r[3] == "correct"]
        c2l = [s for s, r in p.iterrows() if r[2] == "correct" and r[3] != "correct"]
        L.append(f"- {MODELS[model]}: Anderson wrong/never -> continuation correct: seeds {a2c or 'none'}; "
                 f"continuation regressions vs Anderson: {c2a or 'none'}; cold-LM wrong/never -> "
                 f"continuation correct: {l2c or 'none'}; continuation regressions vs cold LM: {c2l or 'none'}.")
    L.append("")

    # ---- 4. selection vs oracle
    L.append("## 4. Objective-selected vs oracle-best (solver failure vs selection failure)")
    L.append("")
    rows = []
    for model in models:
        for method in (1, 2, 3):
            s = t[(t.model == model) & (t.method == method)]
            never = int((s.classBs == "never").sum())
            solverFail = int(((s.oraClassBs == "wrong")).sum())
            selFail = int(((s.oraClassBs == "correct") & (s.classBs == "wrong")).sum())
            gap = (s["selBsErr"] - s["oraBsErr"])
            rows.append([MODELS[model], METHODS[method], never, solverFail, selFail,
                         int((gap > 1e-9).sum()), fmt(gap.median()), fmt(gap.max())])
    L.append(md_table(["model", "method", "never converged", "oracle wrong (solver: no correct endpoint)",
                       "oracle correct but selected wrong (selection failure)",
                       "seeds with selected != oracle error", "median gap", "max gap"], rows))
    L.append("")

    # ---- 5. cost
    L.append("## 5. Cost: Stone-only vs 16-start portfolio, wall-time distributions")
    L.append("")
    L.append("Quantiles shown as min / 25% / median / 75% / 90% / max, seconds.")
    L.append("")
    rows = []
    for model in models:
        for method in (1, 2, 3):
            c = df[(df.model == model) & (df.method == method)]
            s = t[(t.model == model) & (t.method == method)]
            rows.append([MODELS[model], METHODS[method], quant(c["wallTime"]),
                         quant(s["stoneWall"]), quant(s["wallTotal"]), fmt(s["wallTotal"].sum()),
                         fmt(c["iterations"].median()), fmt(c["approxEvals"].median())])
    L.append(md_table(["model", "method", "per-cell wall", "Stone-only wall per seed",
                       "16-start wall per seed", "total wall (s)", "median iterations per cell",
                       "median approx evals per cell"], rows))
    L.append("")
    L.append("Approx evals = iterations x (2k+1) + rejected steps for the direct methods (k = ng*n1); "
             "GLS map evaluations (= iterations) for Anderson. Not an instrumented count.")
    L.append("")

    # ---- 6. failure counters
    L.append("## 6. Failure and safeguard counters (summed over all cells)")
    L.append("")
    rows = []
    for model in models:
        for method in (2, 3):
            c = df[(df.model == model) & (df.method == method)]
            st = dict(c["status"].value_counts().sort_index())
            rows.append([MODELS[model], METHODS[method], ", ".join(f"{k}:{v}" for k, v in st.items()),
                         int(c["covInvFail"].sum()), int(c["stepInvFail"].sum()),
                         int(c["infeasTrial"].sum()), int(c["rhoReject"].sum()),
                         int(c["cAnyStageFailed"].sum()) if method == 3 else "-",
                         int(c["cStepReductions"].sum()) if method == 3 else "-",
                         int(c["cMinStepExhausted"].sum()) if method == 3 else "-",
                         (f"rc3 {int(c['cFailRc3'].sum())}, rc5 {int(c['cFailRc5'].sum())}, "
                          f"rc2 {int(c['cFailRc2'].sum())}, rc4 {int(c['cFailRc4'].sum())}") if method == 3 else "-",
                         int((c["bsOk"] != 1).sum()), int((c["bsSymValid"] != 1).sum())])
        c = df[(df.model == model) & (df.method == 1)]
        rows.append([MODELS[model], "Anderson", ", ".join(f"{k}:{v}" for k, v in
                     dict(c["status"].value_counts().sort_index()).items()),
                     "n/a", "n/a", "n/a", "n/a", "-", "-", "-", "-",
                     int((c["bsOk"] != 1).sum()), int((c["bsSymValid"] != 1).sum())])
    L.append(md_table(["model", "method", "status counts", "cov-inversion fails", "step-inversion fails",
                       "infeasible trials", "rho rejections", "cells with a failed stage",
                       "step reductions", "min-step exhausted", "failed stages by retcode",
                       "bS chain failures", "symmetry projection fallback"], rows))
    L.append("")
    L.append("Status codes: Anderson 0 converged / 2 not converged in 100 iterations. Cold LM: LM retcode "
             "(0 converged, 2 maxIter, 3 infeasible start, 4 stalled, 5 covariance inversion). Continuation: "
             "0 genuine, 1 stage zero failed, 2 min step exhausted, 3 other predicate failed.")
    L.append("")

    # ---- 7. basins
    L.append("## 7. Converged basins per seed (raw structural, tol 0.05)")
    L.append("")
    rows = []
    for model in models:
        for method in (1, 2, 3):
            s = t[(t.model == model) & (t.method == method)]
            rows.append([MODELS[model], METHODS[method], quant(s["nBasins"]), int((s["nBasins"] > 1).sum()),
                         quant(s["nDistinctBench"])])
    L.append(md_table(["model", "method", "basins min/25/50/75/90/max", "seeds with >1 basin",
                       "benchmark-style distinct (bS 1e-4)"], rows))
    L.append("")

    # ---- 8. seed table
    L.append("## 8. Seed-level table (16-start selected; oracle = diagnostic)")
    L.append("")
    rows = []
    for _, r in t.sort_values(["model", "seed", "method"]).iterrows():
        rows.append([MODELS[r.model], r.seed, METHODS[r.method], f"{r.nConv}/{r.nAttempted}",
                     fmt(r.selStart), fmt(r.selBsErr), fmt(r.selRawErr), fmt(r.selCrit, 6),
                     fmt(r.critTrue, 6), fmt(r.oraBsStart), fmt(r.oraBsErr), r.nBasins,
                     r.classBs, fmt(r.stoneBsErr), fmt(r.wallTotal)])
    L.append(md_table(["model", "seed", "method", "conv", "sel start", "sel bsErr", "sel raw",
                       "sel critExact", "truth crit", "oracle start", "oracle bsErr", "basins",
                       "class", "Stone bsErr", "wall (s)"], rows))
    L.append("")

    # ---- 9. anomalies
    L.append("## 9. Anomalies")
    L.append("")
    an = []
    for _, r in t.iterrows():
        tag = f"{MODELS[r.model]} seed {r.seed} {METHODS[r.method]}"
        if r.classBs != r.classRaw:
            an.append((tag, f"primary/secondary disagree: bsErr {fmt(r.selBsErr)} ({r.classBs}) vs raw "
                            f"{fmt(r.selRawErr)} ({r.classRaw}); bS excl. u {fmt(r.selBsErrExclU)}, raw u-row "
                            f"{fmt(r.selUErr)}; selected start {fmt(r.selStart)}"))
        if r.oraClassBs == "correct" and r.classBs == "wrong":
            an.append((tag, f"selection failure: selected start {fmt(r.selStart)} bsErr {fmt(r.selBsErr)} "
                            f"vs oracle start {fmt(r.oraBsStart)} bsErr {fmt(r.oraBsErr)}"))
        if r.classBs != "never" and not math.isnan(r.critTrue) and r.selCrit > r.critTrue + 1e-9 and r.classBs == "wrong":
            an.append((tag, f"wrong endpoint beats truth on the objective: selCrit {fmt(r.selCrit,6)} > "
                            f"truth {fmt(r.critTrue,6)} (bsErr {fmt(r.selBsErr)})"))
        if r.bsFail or r.symInvalid:
            an.append((tag, f"bS chain failures {r.bsFail}, symmetry fallback {r.symInvalid}"))
    # continuation vs cold-LM selected solution disagreement on the same seed
    for (model, seed), g in df[df.method.isin([2, 3])].groupby(["model", "seed"]):
        tt = t[(t.model == model) & (t.seed == seed)].set_index("method")
        if 2 in tt.index and 3 in tt.index and not math.isnan(tt.loc[2, "selStart"]) and not math.isnan(tt.loc[3, "selStart"]):
            e2 = g[(g.method == 2) & (g.start == int(tt.loc[2, "selStart"]))].iloc[0]["raw"]
            e3 = g[(g.method == 3) & (g.start == int(tt.loc[3, "selStart"]))].iloc[0]["raw"]
            d = np.max(np.abs(struct_rows(e2) - struct_rows(e3)))
            if d > BASIN_TOL:
                an.append((f"{MODELS[model]} seed {seed}",
                           f"cold LM and continuation selected different solutions (max raw struct diff {fmt(d)}; "
                           f"critExact {fmt(tt.loc[2,'selCrit'],6)} vs {fmt(tt.loc[3,'selCrit'],6)}; "
                           f"bsErr {fmt(tt.loc[2,'selBsErr'])} vs {fmt(tt.loc[3,'selBsErr'])})"))
    # continuation cells that failed
    for _, r in df[(df.method == 3) & (~df.conv)].iterrows():
        an.append((f"{MODELS[r.model]} seed {r.seed} continuation start {r.start}",
                   f"not genuinely converged: status {r.status}, stage0Ok {fmt(r.cStage0Ok)}, finalHomotopy "
                   f"{fmt(r.finalHomotopy)}, failed stages rc3/5/2/4 {fmt(r.cFailRc3)}/{fmt(r.cFailRc5)}/"
                   f"{fmt(r.cFailRc2)}/{fmt(r.cFailRc4)}, step reductions {fmt(r.cStepReductions)}"))
    if an:
        L.append(md_table(["case", "anomaly"], an))
    else:
        L.append("None.")
    L.append("")
    L.append("Reproduce any cell: `_quaidsSyntheticDGP(3000, seed, model, 1)`, Stone fit, "
             "`_pilotAllStarts(...)`, then the method call in `pilot_harness_core.src` "
             "(or rerun the harness with pilotSeeds/pilotStarts/pilotMethods set to that cell).")

    path = os.path.join(outdir, "pilot_report.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    tc = t.drop(columns=["statusCounts"]).copy()
    tc["model"] = tc["model"].map(MODELS)
    tc["method"] = tc["method"].map(METHODS)
    tc.to_csv(os.path.join(outdir, "pilot_seed_table.csv"), index=False)
    print(f"wrote {path} and pilot_seed_table.csv ({len(t)} seed/model/method rows)")
    return t, df


def compare_runs(dir_a, dir_b):
    """Determinism: every cell present in BOTH runs must agree in every
    metadata column except wallTime, and in the stored endpoint vectors."""
    a = load_cells(dir_a).set_index(["model", "seed", "method", "start"])
    b = load_cells(dir_b).set_index(["model", "seed", "method", "start"])
    common = a.index.intersection(b.index)
    skip = {"wallTime", "raw", "bS", "_line"}
    nbad = 0
    for key in common:
        ra, rb = a.loc[key], b.loc[key]
        for c in a.columns:
            if c in skip:
                continue
            x, y = ra[c], rb[c]
            if not ((pd.isna(x) and pd.isna(y)) or x == y):
                nbad += 1
                print("DIFF", key, c, x, y)
        if not (np.array_equal(ra["raw"], rb["raw"], equal_nan=True)
                and np.array_equal(ra["bS"], rb["bS"], equal_nan=True)):
            nbad += 1
            print("DIFF", key, "endpoint vectors")
    print(f"DETERMINISM: {len(common)} common cells, {nbad} differences (wallTime excluded)")
    return nbad == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--validate")
    ap.add_argument("--report")
    ap.add_argument("--seeds")
    ap.add_argument("--compare", nargs=2)
    a = ap.parse_args()
    rc = 0
    if a.selftest:
        rep = []
        ok = selftest(rep)
        print("\n".join(rep))
        print("SELFTEST:", "ALL PASSED" if ok else "FAILED")
        rc |= 0 if ok else 1
    if a.validate:
        rc |= 0 if validate(a.validate) else 1
    if a.compare:
        rc |= 0 if compare_runs(*a.compare) else 1
    if a.report:
        sf = None
        if a.seeds:
            lo, hi = a.seeds.split("-")
            sf = list(range(int(lo), int(hi) + 1))
        write_report(a.report, sf)
    sys.exit(rc)


if __name__ == "__main__":
    main()
