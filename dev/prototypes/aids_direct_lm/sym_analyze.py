"""sym_analyze.py -- analysis for SYMMETRY_PROTOCOL.md v1 (this directory).

PRIVATE EXPERIMENT ONLY. No solver. Reads the files sym_harness.src writes
(cells.csv, raw.csv, seeds.csv, truth.csv; polish.csv for the real stage),
completes the truth-dependent seeds.csv columns (pfErr, rawStructErr,
classPf), recomputes the selection independently from cells.csv/raw.csv and
fails INVALID on any disagreement with the harness (v1 artifact-schema
clarification), and evaluates the section 7 checks and gates.

Commands (run from this directory):
  python sym_analyze.py analyze --stage fixture --run DIR --out DIR
  python sym_analyze.py analyze --stage dev --run symmetry/dev --det symmetry/dev_det --out symmetry/dev_out
  python sym_analyze.py analyze --stage interim --run symmetry/conf --out symmetry/conf_interim_out
  python sym_analyze.py analyze --stage confirm --run symmetry/conf --det symmetry/conf_det
         --manifest symmetry/SYM_FREEZE_MANIFEST.json --out symmetry/conf_out
  python sym_analyze.py analyze --stage confirm ... --real symmetry/real_out/realdata.json --out symmetry/conf_g6_out
  python sym_analyze.py real --run symmetry/real --out symmetry/real_out
  python sym_analyze.py compare A B          (determinism: 0 diffs apart from wallSec)
  python sym_analyze.py manifest --gauss-version "26.1.4"   (freeze, section 8 step 4)
  python sym_analyze.py verify-manifest

Seed guard: every row's block label and seed are checked against the stage's
hard-coded seeds (sym_config.json) BEFORE any other field of that row is
parsed; a row outside them aborts the run (exit 3).

Interim mode never loads truth.csv and reports no accuracy statistic.
No output is overwritten: analyze/real refuse an --out that already holds
results.json/realdata.json, and manifest refuses an existing manifest and
requires a clean stage D analysis (symmetry/dev_out/results.json) first.
Outputs are deterministic: sorted JSON keys, %.17g floats, no timestamps.
Exit codes: 0 ok; 2 integrity failure (verdict INVALID); 3 guard refusal.
"""
import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(HERE, "sym_config.json")
REPO_ROOT_REL = "../../../.."
DEFAULT_MANIFEST = os.path.join(HERE, "symmetry", "SYM_FREEZE_MANIFEST.json")
DEFAULT_DEV_RESULTS = os.path.join(HERE, "symmetry", "dev_out", "results.json")
ARM_ORDER = {"H": 0, "S": 1, "S1": 2, "P": 3}

KEYS = {
    "cells.csv": ("seed", "arm", "start"),
    "raw.csv": ("seed", "arm", "start", "row", "col"),
    "seeds.csv": ("seed", "arm"),
    "truth.csv": ("seed", "row", "col"),
    "polish.csv": ("seed",),
}


class GuardError(Exception):
    pass


class IntegrityError(Exception):
    pass


def load_config():
    with open(CONFIG_PATH, encoding="utf-8") as fh:
        return json.load(fh)


def allowed_seeds(cfg, stage):
    st = cfg["stages"][stage]
    if "list" in st:
        return set(st["list"])
    lo, hi = st["seeds"]
    return set(range(lo, hi + 1))


def num(tok):
    return float("nan") if tok == "." else float(tok)


def fmt(x):
    if x is None:
        return "."
    if isinstance(x, str):
        return x
    x = float(x)
    if math.isnan(x):
        return "."
    return format(x, ".17g")


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return None if math.isnan(x) or math.isinf(x) else x
    return x


def dump_json(path, obj):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(jsonable(obj), fh, sort_keys=True, indent=1)
        fh.write("\n")


# ----------------------------------------------------------------------------
# reading (guard before parse)
# ----------------------------------------------------------------------------
def read_csv(path, block, allowed):
    """Rows as dicts of raw string tokens. The block label and seed of each
    row are checked BEFORE the rest of the row is split or parsed."""
    rows = []
    with open(path, encoding="utf-8", newline="") as fh:
        header = fh.readline().rstrip("\r\n").split(",")
        if header[:2] != ["block", "seed"]:
            raise IntegrityError(f"{path}: unexpected header {header[:2]}")
        for ln, text in enumerate(fh, 2):
            text = text.rstrip("\r\n")
            if not text:
                continue
            head = text.split(",", 2)
            if head[0] != block:
                raise GuardError(f"{path}:{ln}: block {head[0]!r} is not {block!r}; refusing to read it")
            try:
                seed = int(head[1])
            except ValueError:
                raise GuardError(f"{path}:{ln}: unreadable seed field; refusing to read it")
            if seed not in allowed:
                raise GuardError(f"{path}:{ln}: seed {seed} outside the stage's seeds; refusing to read it")
            toks = text.split(",")
            if len(toks) != len(header):
                raise IntegrityError(f"{path}:{ln}: {len(toks)} fields, header has {len(header)}")
            rows.append(dict(zip(header, toks)))
    return header, rows


def key_of(fname, r):
    return tuple(r[k] if k == "arm" else int(r[k]) for k in KEYS[fname])


def index_rows(fname, rows):
    out = {}
    for r in rows:
        k = key_of(fname, r)
        if k in out:
            raise IntegrityError(f"{fname}: duplicate key {k}")
        out[k] = r
    return out


class Run:
    """One run directory, guarded and indexed."""

    def __init__(self, run_dir, block, allowed, ng, n1, load_truth=True, polish=False):
        self.dir = run_dir
        self.ng, self.n1 = ng, n1
        _, cells = read_csv(os.path.join(run_dir, "cells.csv"), block, allowed)
        _, raws = read_csv(os.path.join(run_dir, "raw.csv"), block, allowed)
        self.seeds_header, seeds = read_csv(os.path.join(run_dir, "seeds.csv"), block, allowed)
        self.cells = index_rows("cells.csv", cells)
        self.seed_rows = index_rows("seeds.csv", seeds)
        self.raw_rows = index_rows("raw.csv", raws)
        self.raw = {}
        for (seed, arm, start, row, col), r in self.raw_rows.items():
            m = self.raw.setdefault((seed, arm, start), np.full((ng, n1), np.nan))
            if not (1 <= row <= ng and 1 <= col <= n1):
                raise IntegrityError(f"raw.csv: row/col out of layout {(seed, arm, start, row, col)}")
            m[row - 1, col - 1] = num(r["value"])
        self.truth = None
        if load_truth:
            _, tr = read_csv(os.path.join(run_dir, "truth.csv"), block, allowed)
            self.truth_rows = index_rows("truth.csv", tr)
            self.truth = {}
            for (seed, row, col), r in self.truth_rows.items():
                t = self.truth.setdefault(seed, {})
                t[(row, col)] = num(r["value"])
        self.polish = None
        if polish:
            _, pr = read_csv(os.path.join(run_dir, "polish.csv"), block, allowed)
            self.polish = index_rows("polish.csv", pr)
        self.seed_list = sorted({k[0] for k in self.seed_rows})

    def sval(self, seed, arm, col):
        return num(self.seed_rows[(seed, arm)][col])


# ----------------------------------------------------------------------------
# metrics (frozen definitions; I6 checks pf_err against confirm_analyze)
# ----------------------------------------------------------------------------
def true_raw(tp, nif, n1):
    n = n1 + 1
    return np.vstack([tp[:nif], tp[nif:nif + n1], tp[nif + n:]])


def unprojected_full(raw, nif):
    n1 = raw.shape[1]
    gam = raw[nif:nif + n1]
    full = np.vstack([raw[:nif], gam, -gam.sum(axis=0, keepdims=True), raw[nif + n1:]])
    full = np.hstack([full, -full.sum(axis=1, keepdims=True)])
    full[0, -1] += 1.0
    return full


def pf_err(raw, tp, nif):
    full = unprojected_full(raw, nif)
    if full.shape != tp.shape:
        raise IntegrityError(f"pfErr shape {full.shape} != truth {tp.shape}")
    return float(np.abs(full[:-1] - tp[:-1]).max())


def raw_struct_err(raw, tp, nif):
    n1 = raw.shape[1]
    tr = true_raw(tp, nif, n1)[:, :n1]
    nstruct = nif + n1 + 1
    return float(np.abs(raw[:nstruct] - tr[:nstruct]).max())


def classify(err, cutoff):
    if err is None or (isinstance(err, float) and math.isnan(err)):
        return "none"
    return "correct" if err <= cutoff else "wrong"


def asym_of(raw, nif):
    n1 = raw.shape[1]
    g = raw[nif:nif + n1]
    return float(np.abs(g - g.T).max())


def basin_dist(a, b):
    return float(np.abs(a[:-1, :] - b[:-1, :]).max())


def cp_lower(k, n, alpha=0.05):
    """Exact Clopper-Pearson two-sided (1 - alpha) lower limit."""
    if k == 0:
        return 0.0
    return float(stats.beta.ppf(alpha / 2, k, n - k + 1))


def same_num(a, b):
    a, b = float(a), float(b)
    return (math.isnan(a) and math.isnan(b)) or a == b


# ----------------------------------------------------------------------------
# selection, recomputed independently (protocol 4.2 / design section 5)
# ----------------------------------------------------------------------------
def select(run, seed, arm, cfg):
    ns = cfg["n_starts"]
    elig, crit = [], []
    for s in range(1, ns + 1):
        c = run.cells[(seed, arm, s)]
        raw = run.raw[(seed, arm, s)]
        cr = num(c["critExact"])
        e = num(c["retcode"]) == 0 and math.isfinite(cr) and bool(np.all(np.isfinite(raw)))
        if float(e) != num(c["eligible"]):
            raise IntegrityError(f"seed {seed} arm {arm} start {s}: eligible {c['eligible']} != recomputed {int(e)}")
        elig.append(e)
        crit.append(cr)
    best = None
    for s in range(ns):
        if elig[s] and (best is None or crit[s] > crit[best]):
            best = s
    if best is None:
        return dict(ok=0, selStart=float("nan"), nEligible=0, crit=float("nan"), exactTies=0,
                    margin=float("nan"), nearTie=0, iterations=float("nan"), support=0)
    braw = run.raw[(seed, arm, best + 1)]
    ties, other, support = 0, None, 1
    for s in range(ns):
        if not elig[s] or s == best:
            continue
        if crit[s] == crit[best]:
            ties += 1
        if basin_dist(run.raw[(seed, arm, s + 1)], braw) > cfg["basin_tol"]:
            if other is None or crit[s] > other:
                other = crit[s]
        else:
            support += 1
    margin = float("nan") if other is None else crit[best] - other
    near = 0 if other is None else int(margin < cfg["near_tie_eps"])
    return dict(ok=1, selStart=best + 1, nEligible=sum(elig), crit=crit[best], exactTies=ties,
                margin=margin, nearTie=near, iterations=num(run.cells[(seed, arm, best + 1)]["iterations"]),
                support=support)


def endpoint(run, seed, arm):
    """Selected (or S1 / P) reduced raw, or None for no result."""
    if run.sval(seed, arm, "ok") != 1:
        return None
    if arm in ("H", "S"):
        return run.raw[(seed, arm, int(run.sval(seed, arm, "selStart")))]
    if arm == "S1":
        return run.raw[(seed, "S", 1)]
    return run.raw[(seed, "P", 0)]


# ----------------------------------------------------------------------------
# integrity (C0 parts computable from files)
# ----------------------------------------------------------------------------
def check_grid(run, seeds, cfg, with_truth):
    ns, ng, n1 = cfg["n_starts"], run.ng, run.n1
    probs = []
    if sorted(run.seed_list) != sorted(seeds):
        probs.append(f"seeds present {run.seed_list} != expected {sorted(seeds)}")
    for seed in seeds:
        for arm in cfg["arms_seeds"]:
            if (seed, arm) not in run.seed_rows:
                probs.append(f"seeds.csv missing ({seed}, {arm})")
        for arm in cfg["arms_cells"]:
            for s in range(1, ns + 1):
                if (seed, arm, s) not in run.cells:
                    probs.append(f"cells.csv missing ({seed}, {arm}, {s})")
                if (seed, arm, s) not in run.raw:
                    probs.append(f"raw.csv missing ({seed}, {arm}, {s})")
        p_ok = (seed, "P") in run.seed_rows and run.sval(seed, "P", "ok") == 1
        if p_ok != ((seed, "P", 0) in run.raw):
            probs.append(f"raw.csv P rows inconsistent with P ok for seed {seed}")
    n_p = sum(1 for k in run.raw if k[1] == "P")
    exp_cells = len(seeds) * len(cfg["arms_cells"]) * ns
    exp_raw = (len(seeds) * len(cfg["arms_cells"]) * ns + n_p) * ng * n1
    exp_seeds = len(seeds) * len(cfg["arms_seeds"])
    if len(run.cells) != exp_cells:
        probs.append(f"cells.csv has {len(run.cells)} rows, expected {exp_cells}")
    if len(run.raw_rows) != exp_raw:
        probs.append(f"raw.csv has {len(run.raw_rows)} rows, expected {exp_raw}")
    if len(run.seed_rows) != exp_seeds:
        probs.append(f"seeds.csv has {len(run.seed_rows)} rows, expected {exp_seeds}")
    if with_truth:
        lay = cfg["layout"]["synthetic"]
        exp_truth = len(seeds) * lay["truth_rows"] * lay["truth_cols"]
        if len(run.truth_rows) != exp_truth:
            probs.append(f"truth.csv has {len(run.truth_rows)} rows, expected {exp_truth}")
    return probs


def check_selection(run, seed, cfg, nobs):
    """Harness seeds.csv vs the independent recomputation. Returns
    (problems, {arm: selection dict})."""
    probs, sel = [], {}
    for arm in ("H", "S"):
        r = select(run, seed, arm, cfg)
        sel[arm] = r
        for col in ("ok", "selStart", "nEligible", "crit", "exactTies", "margin", "nearTie", "iterations"):
            if not same_num(r[col], run.sval(seed, arm, col)):
                probs.append(f"seed {seed} arm {arm}: {col} harness {run.seed_rows[(seed, arm)][col]} != recomputed {fmt(r[col])}")
    c1 = run.cells[(seed, "S", 1)]
    ok1 = num(c1["eligible"])
    if not same_num(ok1, run.sval(seed, "S1", "ok")):
        probs.append(f"seed {seed} S1: ok mismatch")
    if ok1 == 1 and not same_num(num(c1["critExact"]), run.sval(seed, "S1", "crit")):
        probs.append(f"seed {seed} S1: crit mismatch")
    ch, cs, cp = (run.sval(seed, a, "crit") for a in ("H", "S", "P"))
    n1g, n2g, lr = cs - ch, cp - cs, nobs * (ch - cs)
    for col, val in (("N1", n1g), ("N2", n2g), ("LR", lr)):
        if not same_num(val, run.sval(seed, "S", col)):
            probs.append(f"seed {seed}: {col} harness {run.seed_rows[(seed, 'S')][col]} != recomputed {fmt(val)}")
    for arm in cfg["arms_seeds"]:
        e = endpoint(run, seed, arm)
        a = float("nan") if e is None else asym_of(e, run.ng - run.n1 - 2)
        if not same_num(a, run.sval(seed, arm, "asym")):
            probs.append(f"seed {seed} arm {arm}: asym harness {run.seed_rows[(seed, arm)]['asym']} != recomputed {fmt(a)}")
    return probs, sel


def c0_numeric(run, seeds, cfg):
    """asym == 0 (S, P), adding-up/homogeneity, critAbsErr (S)."""
    probs = []
    for seed in seeds:
        for arm in cfg["arms_seeds"]:
            if run.sval(seed, arm, "ok") != 1:
                continue
            if arm in ("S", "S1", "P") and run.sval(seed, arm, "asym") != 0:
                probs.append(f"seed {seed} arm {arm}: asym {run.seed_rows[(seed, arm)]['asym']} != 0")
            for col in ("addUpErr", "homogErr"):
                if not run.sval(seed, arm, col) <= cfg["tol_addup_homog"]:
                    probs.append(f"seed {seed} arm {arm}: {col} {run.seed_rows[(seed, arm)][col]} > tol")
        if run.sval(seed, "S", "ok") == 1 and not run.sval(seed, "S", "critAbsErr") <= cfg["tol_crit_abs"]:
            probs.append(f"seed {seed}: S critAbsErr {run.seed_rows[(seed, 'S')]['critAbsErr']} > {cfg['tol_crit_abs']}")
    return probs


def compare_dirs(a_dir, b_dir):
    """Determinism: every row of B (the re-solve) must equal the same-key
    row of A in every field except block and wallSec. Returns a list of
    differences (empty = identical)."""
    diffs = []
    for fname in ("cells.csv", "raw.csv", "seeds.csv", "truth.csv", "polish.csv"):
        pa, pb = os.path.join(a_dir, fname), os.path.join(b_dir, fname)
        if os.path.exists(pa) != os.path.exists(pb):
            diffs.append(f"{fname}: present in only one of the two runs")
            continue
        if not os.path.exists(pb):
            continue
        rows = {}
        for d in (a_dir, b_dir):
            with open(os.path.join(d, fname), encoding="utf-8", newline="") as fh:
                header = fh.readline().rstrip("\r\n").split(",")
                rr = {}
                for text in fh:
                    text = text.rstrip("\r\n")
                    if text:
                        r = dict(zip(header, text.split(",")))
                        rr[key_of(fname, r)] = r
                rows[d] = rr
        ra, rb = rows[a_dir], rows[b_dir]
        b_seeds = {k[0] for k in rb}
        a_sub = {k for k in ra if k[0] in b_seeds}
        if not rb:
            diffs.append(f"{fname}: the re-solve has no rows")
        if a_sub != set(rb):
            diffs.append(f"{fname}: key sets differ for seeds {sorted(b_seeds)}")
        for k, r in rb.items():
            if k not in ra:
                continue
            for col, v in r.items():
                if col in ("block", "wallSec"):
                    continue
                if ra[k].get(col) != v:
                    diffs.append(f"{fname} {k} {col}: {ra[k].get(col)} != {v}")
    return diffs


# ----------------------------------------------------------------------------
# manifest (section 8 step 4)
# ----------------------------------------------------------------------------
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def repo_root():
    return os.path.normpath(os.path.join(os.path.dirname(DEFAULT_MANIFEST), REPO_ROOT_REL))


def make_manifest(cfg, gauss_version, path=DEFAULT_MANIFEST, dev_results=DEFAULT_DEV_RESULTS):
    """Freeze (section 8 step 4). Refuses (GuardError) to overwrite an
    existing manifest or sidecar, and requires the stage D analysis to exist
    with no integrity problem and no stop (revise) rule triggered."""
    if os.path.exists(path) or os.path.exists(path[:-5] + ".sha256"):
        raise GuardError(f"{path} already exists; the freeze is written once and never overwritten")
    if not os.path.exists(dev_results):
        raise GuardError(f"no stage D analysis at {dev_results}; the freeze follows stage D")
    with open(dev_results, encoding="utf-8") as fh:
        dev = json.load(fh)
    if dev.get("mode") != "dev" or dev["integrity"]["problems"] or dev.get("devStop", {}).get("stopRevise", True):
        raise GuardError(f"{dev_results} is not a clean stage D pass (integrity problem or stop rule); no freeze")
    root = repo_root()
    files = {rel: sha256_file(os.path.join(root, rel)) for rel in cfg["freeze_files"]}
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True).stdout.strip()
    import scipy
    man = dict(protocol=cfg["protocol"], repo_root_rel=REPO_ROOT_REL, files=files, git_head=head,
               versions=dict(gauss=gauss_version, python=platform.python_version(),
                             numpy=np.__version__, scipy=scipy.__version__))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    dump_json(path, man)
    with open(path[:-5] + ".sha256", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"{sha256_file(path)}  {os.path.basename(path)}\n")
    return man


def verify_manifest(path=DEFAULT_MANIFEST):
    bad = []
    side = path[:-5] + ".sha256"
    if not os.path.exists(path) or not os.path.exists(side):
        return ["manifest or its .sha256 sidecar missing"]
    with open(side, encoding="utf-8") as fh:
        if fh.read().split()[0] != sha256_file(path):
            bad.append("manifest digest does not match its sidecar")
    with open(path, encoding="utf-8") as fh:
        man = json.load(fh)
    root = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), man["repo_root_rel"]))
    for rel, digest in sorted(man["files"].items()):
        p = os.path.join(root, rel)
        if not os.path.exists(p):
            bad.append(f"missing {rel}")
        elif sha256_file(p) != digest:
            bad.append(f"hash mismatch {rel}")
    return bad


# ----------------------------------------------------------------------------
# verdict (section 10; first matching row wins)
# ----------------------------------------------------------------------------
def verdict(c0, c1, c2, c3, c4, c5, g6, p_passes_c1_c2):
    if not c0:
        return "INVALID"
    if c1 and c2 and c3 and c4 and c5 and g6 is True:
        return "PROCEED"
    if not (c1 and c2):
        return "REVISE" if p_passes_c1_c2 else "REJECT"
    if not (c3 and c4 and c5):
        return "REVISE"
    if g6 is None:
        return "PENDING_G6"
    return "REVISE"


# ----------------------------------------------------------------------------
# analyze (fixture / dev / interim / confirm)
# ----------------------------------------------------------------------------
def analyze(stage, run_dir, out_dir, det_dir=None, manifest=None, real_json=None, cfg=None):
    cfg = cfg or load_config()
    st = cfg["stages"][stage]
    allowed = allowed_seeds(cfg, stage)
    lay = cfg["layout"]["synthetic"]
    ng, n1, nif, nobs = lay["ng"], lay["n1"], lay["n_int_full"], lay["nobs"]
    interim = stage == "interim"
    if os.path.exists(os.path.join(out_dir, "results.json")):
        raise GuardError(f"{out_dir} already holds results.json; analysis outputs are never overwritten")
    run = Run(run_dir, st["block"], allowed, ng, n1, load_truth=not interim)
    seeds = sorted(allowed)

    integrity = []
    manifest_ok = None
    if stage == "confirm":
        bad = verify_manifest(manifest or DEFAULT_MANIFEST)
        manifest_ok = not bad
        integrity += [f"manifest: {b}" for b in bad]
    integrity += check_grid(run, seeds, cfg, with_truth=not interim)
    sel = {}
    if not integrity:
        for seed in seeds:
            p, sel[seed] = check_selection(run, seed, cfg, nobs)
            integrity += p
        integrity += c0_numeric(run, seeds, cfg)
    det_diffs = None
    if det_dir is not None:
        det_stage = {"dev": "dev_det", "confirm": "conf_det"}.get(stage)
        if det_stage is not None:
            det_seeds = sorted(allowed_seeds(cfg, det_stage))
            det_run = Run(det_dir, cfg["stages"][det_stage]["block"], set(det_seeds), ng, n1)
            integrity += [f"determinism: re-solve {p}" for p in check_grid(det_run, det_seeds, cfg, with_truth=True)]
        det_diffs = compare_dirs(run_dir, det_dir)
        integrity += [f"determinism: {d}" for d in det_diffs]
    elif stage in ("dev", "confirm"):
        integrity.append("determinism re-solve directory (--det) not supplied")

    # runtime (all stages; never a decision input)
    walls = []
    slow = []
    for seed in run.seed_list:
        w = 0.0
        for arm in ("H", "S", "P"):
            v = run.sval(seed, arm, "wallSec")
            if math.isfinite(v):
                w += v
                if v > cfg["runtime"]["slow_cell_s"]:
                    slow.append([seed, arm, v])
        walls.append(w)
    runtime = dict(median_seed_s=float(np.median(walls)) if walls else None,
                   p90_seed_s=float(np.percentile(walls, 90)) if walls else None,
                   max_seed_s=float(max(walls)) if walls else None, slow_arms=slow)
    if walls:
        proj = runtime["median_seed_s"] * cfg["runtime"]["n_confirm_seeds"]
        runtime["projected_stage_c_s"] = proj
        runtime["exceeds_cap"] = proj > cfg["runtime"]["cap_hours"] * 3600

    os.makedirs(out_dir, exist_ok=True)
    res = dict(protocolVersion="v1", mode=stage, block=st["block"], nSeeds=len(run.seed_list),
               integrity=dict(problems=integrity, manifestOk=manifest_ok,
                              determinismDiffs=None if det_diffs is None else len(det_diffs)),
               runtime=runtime)
    if interim:
        res["verdict"] = None
        res["note"] = "interim look: integrity and runtime only; no accuracy statistic is computed"
        dump_json(os.path.join(out_dir, "results.json"), res)
        write_report(os.path.join(out_dir, "report.md"), res)
        return res, (2 if integrity else 0)
    if integrity:
        res["verdict"] = "INVALID"
        if stage == "dev":
            # section 8 step 3: asym != 0 and a determinism difference are
            # stop (revise) reasons, not only integrity failures
            stop = [p for p in integrity if p.startswith("determinism:") or "asym" in p]
            res["devStop"] = dict(stopRevise=bool(stop) or None, reasons=stop,
                                  note="INVALID run: the other stop rules are not evaluated")
        dump_json(os.path.join(out_dir, "results.json"), res)
        write_report(os.path.join(out_dir, "report.md"), res)
        return res, 2

    # completed seeds.csv and per-seed metrics
    cut = cfg["primary_cutoff"]
    per = {}
    lines = [",".join(run.seeds_header)]
    for seed in seeds:
        t = run.truth[seed]
        tp = np.array([[t[(r, c)] for c in range(1, lay["truth_cols"] + 1)] for r in range(1, lay["truth_rows"] + 1)])
        for arm in sorted(cfg["arms_seeds"], key=ARM_ORDER.get):
            e = endpoint(run, seed, arm)
            pf = None if e is None else pf_err(e, tp, nif)
            rs = None if e is None else raw_struct_err(e, tp, nif)
            per[(seed, arm)] = dict(pf=pf, cls=classify(pf, cut))
            r = dict(run.seed_rows[(seed, arm)])
            r["pfErr"], r["rawStructErr"], r["classPf"] = fmt(pf), fmt(rs), classify(pf, cut)
            lines.append(",".join(r[h] for h in run.seeds_header))
    with open(os.path.join(out_dir, "seeds.csv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    def cls(seed, arm, c=cut):
        pf = per[(seed, arm)]["pf"]
        return classify(pf, c)

    def count(arm, which, c=cut):
        return sum(1 for s in seeds if cls(s, arm, c) == which)

    n = len(seeds)
    g = cfg["gates"]
    k_s, w_s = count("S", "correct"), count("S", "wrong")
    n2_vals = [run.sval(s, "S", "N2") for s in seeds]
    n2_flags = [s for s, v in zip(seeds, n2_vals) if math.isfinite(v) and v > cfg["n2_flag"]]
    n2_missing = sum(1 for v in n2_vals if not math.isfinite(v))
    c_cnt = sum(1 for s in seeds if cls(s, "P") == "correct" and cls(s, "S") != "correct")
    b_cnt = sum(1 for s in seeds if cls(s, "S") == "correct" and cls(s, "P") != "correct")
    regress = [s for s in seeds if cls(s, "H") == "correct" and cls(s, "S") != "correct"]
    k_p, w_p = count("P", "correct"), count("P", "wrong")
    gates = {
        "C0": dict(statistic=len(integrity), threshold="no integrity problem", **{"pass": not integrity}),
        "C1": dict(statistic=k_s, cpLower=cp_lower(k_s, n, g["C1"]["alpha_two_sided"]),
                   threshold=g["C1"]["min_correct"], **{"pass": k_s >= g["C1"]["min_correct"]}),
        "C2": dict(statistic=w_s, threshold=g["C2"]["max_returned_wrong"], **{"pass": w_s <= g["C2"]["max_returned_wrong"]}),
        "C3": dict(statistic=len(n2_flags), seeds=n2_flags, n2Missing=n2_missing,
                   threshold=g["C3"]["max_n2_flags"], **{"pass": len(n2_flags) <= g["C3"]["max_n2_flags"]}),
        "C4": dict(statistic=c_cnt - b_cnt, c=c_cnt, b=b_cnt, threshold=g["C4"]["max_c_minus_b"],
                   **{"pass": c_cnt - b_cnt <= g["C4"]["max_c_minus_b"]}),
        "C5": dict(statistic=len(regress), seeds=regress, threshold=g["C5"]["max_regressions"],
                   mechanisms="a written mechanism is required for each listed seed (report.md)",
                   **{"pass": len(regress) <= g["C5"]["max_regressions"]}),
    }
    g6 = None
    if real_json is not None:
        with open(real_json, encoding="utf-8") as fh:
            g6 = bool(json.load(fh)["G6"]["pass"])
    gates["G6"] = dict(statistic=None, threshold="all G6 conditions", **{"pass": g6})
    p_pass = k_p >= g["C1"]["min_correct"] and w_p <= g["C2"]["max_returned_wrong"]

    # descriptive (section 7.3)
    pairs = [(per[(s, "S")]["pf"], per[(s, "H")]["pf"]) for s in seeds
             if per[(s, "S")]["pf"] is not None and per[(s, "H")]["pf"] is not None]
    d = np.array([a - b for a, b in pairs])
    wil = None
    if len(d) and np.any(d != 0):
        try:
            wil = float(stats.wilcoxon(d, zero_method="wilcox", method="auto", alternative="two-sided").pvalue)
        except ValueError:
            wil = None
    lr_crit = float(stats.chi2.ppf(0.95, cfg["lr_df"]["synthetic"]))
    lrs = [run.sval(s, "S", "LR") for s in seeds]
    lrs = [v for v in lrs if math.isfinite(v)]
    sens = {arm: {fmt(c): {k: count(arm, k, c) for k in ("correct", "wrong", "none")}
                  for c in cfg["sensitivity_cutoffs"]} for arm in cfg["arms_seeds"]}
    lo, hi = cfg["borderline_band"]
    border = [[s, a, per[(s, a)]["pf"]] for s in seeds for a in cfg["arms_seeds"]
              if per[(s, a)]["pf"] is not None and lo < per[(s, a)]["pf"] <= hi]
    descriptive = dict(
        k_P=k_p, W_P=w_p, k_H=count("H", "correct"), W_H=count("H", "wrong"),
        k_S1=count("S1", "correct"), W_S1=count("S1", "wrong"),
        none={a: count(a, "none") for a in cfg["arms_seeds"]},
        pairedPfErrSminusH=dict(n=len(d), median=float(np.median(d)) if len(d) else None, wilcoxonP=wil,
                                note="context only"),
        sensitivity=sens, borderline=border,
        n1Flags=[s for s in seeds if run.sval(s, "S", "N1") > cfg["n1_flag"]],
        lrExceedance=dict(chi2_df=cfg["lr_df"]["synthetic"], crit95=lr_crit,
                          rate=(sum(v > lr_crit for v in lrs) / len(lrs)) if lrs else None,
                          note="preview only; validates no test"),
        support={str(s): sel[s]["S"]["support"] for s in seeds},
    )
    v = verdict(True, gates["C1"]["pass"], gates["C2"]["pass"], gates["C3"]["pass"], gates["C4"]["pass"],
                gates["C5"]["pass"], g6, p_pass)
    res.update(gates=gates, descriptive=descriptive, pPassesC1C2=p_pass)
    if stage == "confirm":
        res["verdict"] = v
        res["manifestSha256"] = sha256_file(manifest or DEFAULT_MANIFEST)
    else:
        res["verdict"] = None
        res["informationalVerdict"] = v
    if stage == "dev":
        ds = cfg["dev_stop"]
        no_result = count("S", "none")
        stop = []
        if w_s > ds["max_returned_wrong"]:
            stop.append(f"S returned-wrong seeds: {w_s}")
        if no_result > ds["max_no_result"]:
            stop.append(f"S no-result seeds: {no_result}")
        if len(n2_flags) > ds["max_n2_flags"]:
            stop.append(f"N2 > {cfg['n2_flag']} on seeds {n2_flags}")
        res["devStop"] = dict(stopRevise=bool(stop), reasons=stop)
    dump_json(os.path.join(out_dir, "results.json"), res)
    write_report(os.path.join(out_dir, "report.md"), res)
    return res, 0


def write_report(path, res):
    out = [f"# Symmetry protocol v1 -- {res['mode']} analysis", "",
           f"- block: {res['block']}; seeds: {res['nSeeds']}",
           f"- verdict: {res.get('verdict')}" + (f" (informational: {res['informationalVerdict']})"
                                                 if res.get("informationalVerdict") else ""), ""]
    out += ["## Integrity", ""]
    probs = res["integrity"]["problems"]
    out += [f"- {p}" for p in probs] or ["- no problem found"]
    out += ["", "## Runtime (never a decision input)", "", "```", json.dumps(jsonable(res["runtime"]), sort_keys=True, indent=1), "```"]
    if "gates" in res:
        out += ["", "## Gates", "", "| Gate | Statistic | Threshold | Pass |", "|---|---|---|---|"]
        for k in sorted(res["gates"]):
            gt = res["gates"][k]
            out.append(f"| {k} | {gt['statistic']} | {gt['threshold']} | {gt['pass']} |")
        if res["gates"]["C5"]["seeds"]:
            out += ["", "## C5 mechanisms (to be written)", ""]
            out += [f"- seed {s}: MECHANISM REQUIRED" for s in res["gates"]["C5"]["seeds"]]
        out += ["", "## Descriptive (never pass/fail)", "", "```",
                json.dumps(jsonable(res["descriptive"]), sort_keys=True, indent=1), "```"]
    if "devStop" in res:
        out += ["", "## Stage D stop rule", "", f"- stop (revise): {res['devStop']['stopRevise']}"]
        out += [f"- {r}" for r in res["devStop"]["reasons"]]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")


# ----------------------------------------------------------------------------
# real data (section 6.2, gate G6)
# ----------------------------------------------------------------------------
def analyze_real(run_dir, out_dir, cfg=None):
    cfg = cfg or load_config()
    lay = cfg["layout"]["real"]
    ng, n1, nif, nobs = lay["ng"], lay["n1"], lay["n_int_full"], lay["nobs"]
    if os.path.exists(os.path.join(out_dir, "realdata.json")):
        raise GuardError(f"{out_dir} already holds realdata.json; analysis outputs are never overwritten")
    run = Run(run_dir, "real", {0}, ng, n1, load_truth=False, polish=True)
    probs = check_grid(run, [0], cfg, with_truth=False)
    sel = {}
    if not probs:
        probs, sel = check_selection(run, 0, cfg, nobs)
    sv = lambda arm, col: run.sval(0, arm, col)
    pol = run.polish[(0,)]
    rc, it, cd = num(pol["retcode"]), num(pol["iterations"]), num(pol["critDiff"])
    r5 = {arm: dict(asym=sv(arm, "asym"), addUpErr=sv(arm, "addUpErr"), homogErr=sv(arm, "homogErr"))
          for arm in ("S", "P")}
    eS, eH, eP = (endpoint(run, 0, a) for a in ("S", "H", "P"))
    r8 = dict(
        SminusH=None if eS is None or eH is None else (eS - eH).tolist(),
        PminusS=None if eP is None or eS is None else (eP - eS).tolist(),
        ownPriceGammaS=None if eS is None else [float(eS[nif + i, i]) for i in range(n1)])
    lr = sv("S", "LR")
    rd = dict(
        ok=dict(H=sv("H", "ok"), S=sv("S", "ok"), P=sv("P", "ok")),
        R1=sv("S", "N1"), R2=sv("S", "N2"),
        R3=dict(LR=lr, chi2_df=cfg["lr_df"]["real"],
                p=float(stats.chi2.sf(lr, cfg["lr_df"]["real"])) if math.isfinite(lr) else None,
                note="non-inferential: u is a generated regressor and no size was validated"),
        R4=dict(support=sel.get("S", {}).get("support"), margin=sv("S", "margin"), nearTie=sv("S", "nearTie"),
                exactTies=sv("S", "exactTies")),
        R5=r5, R6=dict(retcode=rc, iterations=it, critDiff=cd), R7=sv("S", "condGN"), R8=r8,
        integrity=probs, note="truth-free falsification checks only; none is a correctness claim")
    tol = cfg["tol_addup_homog"]
    r5_ok = all(r5[a]["asym"] == 0 and r5[a]["addUpErr"] <= tol and r5[a]["homogErr"] <= tol for a in ("S", "P"))
    r6_ok = rc == 0 and it == 1 and cd <= cfg["tol_polish"]
    conds = dict(S_ok=sv("S", "ok") == 1, P_ok=sv("P", "ok") == 1, R1=rd["R1"] <= cfg["n1_flag"],
                 R2=rd["R2"] <= cfg["n2_flag"], R5=r5_ok, R6=r6_ok, integrity=not probs)
    rd["G6"] = dict(conditions=conds, **{"pass": all(conds.values())})
    os.makedirs(out_dir, exist_ok=True)
    dump_json(os.path.join(out_dir, "realdata.json"), rd)
    return rd, (2 if probs else 0)


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("analyze")
    a.add_argument("--stage", required=True, choices=["fixture", "dev", "interim", "confirm"])
    a.add_argument("--run", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--det")
    a.add_argument("--manifest")
    a.add_argument("--real")
    r = sub.add_parser("real")
    r.add_argument("--run", required=True)
    r.add_argument("--out", required=True)
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    m = sub.add_parser("manifest")
    m.add_argument("--gauss-version", required=True)
    sub.add_parser("verify-manifest")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "analyze":
            if args.stage == "confirm" and verify_manifest(args.manifest or DEFAULT_MANIFEST):
                print("REFUSED: confirm mode needs a verified freeze manifest")
                for b in verify_manifest(args.manifest or DEFAULT_MANIFEST):
                    print("  " + b)
                return 3
            res, code = analyze(args.stage, args.run, args.out, args.det, args.manifest, args.real)
            print(f"{args.stage}: verdict {res.get('verdict')} informational {res.get('informationalVerdict')} "
                  f"integrity problems {len(res['integrity']['problems'])}")
            return code
        if args.cmd == "real":
            rd, code = analyze_real(args.run, args.out)
            print(f"real: G6 pass {rd['G6']['pass']}")
            return code
        if args.cmd == "compare":
            diffs = compare_dirs(args.a, args.b)
            for d in diffs:
                print(d)
            print(f"{len(diffs)} differences (wallSec and block excluded)")
            return 0 if not diffs else 2
        if args.cmd == "manifest":
            make_manifest(load_config(), args.gauss_version)
            print(f"wrote {DEFAULT_MANIFEST}")
            return 0
        if args.cmd == "verify-manifest":
            bad = verify_manifest()
            for b in bad:
                print(b)
            print("manifest OK" if not bad else "manifest FAILED")
            return 0 if not bad else 2
    except GuardError as e:
        print(f"GUARD REFUSAL: {e}")
        return 3
    except IntegrityError as e:
        print(f"INTEGRITY FAILURE: {e}")
        return 2
    return 1


if __name__ == "__main__":
    sys.exit(main())
