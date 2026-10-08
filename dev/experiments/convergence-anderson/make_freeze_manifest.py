"""make_freeze_manifest.py -- SHA-256 manifest of the frozen analysis inputs.

EXPERIMENT ONLY. See CORRECTNESS_VALIDATION_PROTOCOL.md section 14.5.

  python make_freeze_manifest.py --write    # (re)write FREEZE_MANIFEST.json
  python make_freeze_manifest.py --verify   # exit 0 and print MANIFEST OK if every hash matches

The manifest has no timestamp, so regenerating it from unchanged files gives
a byte-identical file. Paths are relative to the repository root.
confirm_analyze.py --mode confirm refuses to run unless this manifest
verifies.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_REL = "../../.."
ROOT = os.path.normpath(os.path.join(HERE, ROOT_REL))
EXP = "dev/experiments/convergence-anderson"
MANIFEST = os.path.join(HERE, "FREEZE_MANIFEST.json")

GROUPS = {
    "protocol": [f"{EXP}/CORRECTNESS_VALIDATION_PROTOCOL.md", f"{EXP}/PILOT_PROTOCOL.md"],
    "analysis": [f"{EXP}/confirm_analyze.py", f"{EXP}/confirm_config.json", f"{EXP}/test_confirm_analyze.py",
                 f"{EXP}/pilot_analyze.py", f"{EXP}/make_freeze_manifest.py"],
    "harness_executed_by_E2": [f"{EXP}/pilot_harness_core.src", f"{EXP}/pilot_lib.src", f"{EXP}/direct_objective.src",
                               f"{EXP}/quaidsfit_anderson_prototype.src",
                               "src/quaids.sdf", "src/quaidsutil.src", "src/quaidsiv.src", "src/quaidselas.src",
                               "src/quaidsslutzky.src", "src/quaids.src", "tests/quaidsfixtures.src"],
    "E2_entry_scripts": [f"{EXP}/e2_run_m0.e", f"{EXP}/e2_run_m1.e", f"{EXP}/e2_truth_dump.e",
                         f"{EXP}/e2_determinism.e", f"{EXP}/pilot_resume_test.e", f"{EXP}/pilot_resume_test.ps1"],
    "reference_outputs": [f"{EXP}/anderson_full_multistart_benchmark_output.txt"],
    "dev_analysis_inputs": [f"{EXP}/pilot_out/columns.txt", f"{EXP}/pilot_out/cells_m0.csv",
                            f"{EXP}/pilot_out/cells_m1.csv", f"{EXP}/e0_out/truth_m0.csv",
                            f"{EXP}/e0_out/truth_m1.csv", f"{EXP}/e0_out/e0_cells_pferr.csv"],
    "dev_canonical_outputs": [f"{EXP}/confirm_dev_out/results.json", f"{EXP}/confirm_dev_out/selected.csv",
                              f"{EXP}/confirm_dev_out/report.md"],
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def environment():
    import numpy, pandas, scipy  # noqa: E401
    try:
        head = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception as e:  # pragma: no cover
        head = f"unavailable ({e})"
    gauss = "unknown"
    banner = os.path.join(HERE, "pilot_validation", "pilot_validate_console.txt")
    if os.path.exists(banner):
        with open(banner, encoding="utf-8-sig", errors="replace") as fh:
            for text in fh:
                text = text.strip()
                if text:
                    if text.startswith("GAUSS"):
                        gauss = text
                    break
    return {"git_head": head, "python": sys.version.split()[0], "numpy": numpy.__version__,
            "scipy": scipy.__version__, "pandas": pandas.__version__, "gauss_banner_from_stage_A": gauss}


def build():
    files, groups = {}, {}
    for g, paths in sorted(GROUPS.items()):
        groups[g] = sorted(paths)
        for rel in paths:
            p = os.path.join(ROOT, rel)
            if not os.path.exists(p):
                raise SystemExit(f"missing file for manifest: {rel}")
            files[rel] = sha256(p)
    return {
        "manifest_version": 1,
        "status": ("PRE-FREEZE CANDIDATE until the protocol's section 15 is ticked; regenerate with --write "
                   "after the approval edit -- that regenerated manifest is the frozen record."),
        "repo_root_rel": ROOT_REL,
        "environment": environment(),
        "groups": groups,
        "files": dict(sorted(files.items())),
    }


def write():
    man = build()
    with open(MANIFEST, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(f"wrote {MANIFEST}: {len(man['files'])} files; manifest sha256 {sha256(MANIFEST)}")


def verify():
    with open(MANIFEST, encoding="utf-8") as fh:
        man = json.load(fh)
    bad = []
    for rel, digest in sorted(man["files"].items()):
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            bad.append(f"MISSING {rel}")
        elif sha256(p) != digest:
            bad.append(f"MISMATCH {rel}")
    if bad:
        print("MANIFEST FAILED:")
        print("\n".join("  " + b for b in bad))
        return 1
    print(f"MANIFEST OK: {len(man['files'])} files verified; manifest sha256 {sha256(MANIFEST)}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true")
    g.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    if a.write:
        write()
    else:
        sys.exit(verify())
