"""Reproduce Phase D0 from committed Phase C/Q0 CSV evidence."""

from __future__ import annotations

import csv
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CONV = ROOT / "dev/experiments/tvp_conv_study/results/categories.csv"
Q0 = ROOT / "dev/experiments/tvp_q_identification/results"
OUT = HERE / "results"


def read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = required - fields
        if missing:
            raise RuntimeError(f"{path}: missing columns {sorted(missing)}")
        return list(reader)


def num(value: str) -> float:
    return float(value)


def fmt(value: float | None, digits: int = 6) -> str:
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.{digits}g}"


def classify(
    *, all_rc0: bool, max_ll_gap: float, tol: float,
    n_boundary: int, n_flat_eig: int, n_floor_under_marker: int,
) -> str:
    """Apply the frozen precedence; thresholds are descriptive, not tests."""
    if not all_rc0:
        return "FIT-FAILED"
    if max_ll_gap > tol:
        return "START-SENSITIVE"
    if n_boundary > 0 or n_flat_eig > 0 or n_floor_under_marker > 0:
        return "GEOMETRY-CAUTION"
    return "NO-CAUTION-OBSERVED"


def main() -> None:
    conv = read_csv(CONV, {"cell_id", "category"})
    endpoints = read_csv(
        Q0 / "endpoint_summary.csv",
        {
            "endpoint", "k", "tol", "all_starts_rc0", "max_start_ll_gap",
            "max_start_dlog10", "n_boundary", "n_floor_under_marker",
            "loss_floor_all", "n_flat_eig",
        },
    )
    fits = read_csv(Q0 / "fits.csv", {"endpoint", "qmult", "rc"})
    hints = read_csv(Q0 / "hint.csv", {"endpoint", "step", "cond"})

    if len(endpoints) != 5:
        raise RuntimeError(f"expected 5 Q0 endpoints, found {len(endpoints)}")
    endpoint_names = {row["endpoint"] for row in endpoints}
    if {row["endpoint"] for row in fits} != endpoint_names:
        raise RuntimeError("fits.csv endpoint set does not match endpoint_summary.csv")
    fit_counts = {name: 0 for name in endpoint_names}
    for row in fits:
        fit_counts[row["endpoint"]] += 1
    if any(count != 3 for count in fit_counts.values()):
        raise RuntimeError(f"expected exactly 3 fits per endpoint, found {fit_counts}")

    conv_global = sum(row["category"] == "CONV-GLOBAL" for row in conv)
    if conv_global != len(conv):
        raise RuntimeError("Phase C categories input contains a non-CONV-GLOBAL row")

    cond = {
        row["endpoint"]: num(row["cond"])
        for row in hints
        if abs(num(row["step"]) - 1.0) < 1e-12
    }
    rows: list[dict[str, str]] = []
    for row in sorted(endpoints, key=lambda x: x["endpoint"]):
        name = row["endpoint"]
        all_rc0 = row["all_starts_rc0"].lower() == "true"
        label = classify(
            all_rc0=all_rc0,
            max_ll_gap=num(row["max_start_ll_gap"]),
            tol=num(row["tol"]),
            n_boundary=int(float(row["n_boundary"])),
            n_flat_eig=int(float(row["n_flat_eig"])),
            n_floor_under_marker=int(float(row["n_floor_under_marker"])),
        )
        rows.append(
            {
                "endpoint": name,
                "k": str(int(float(row["k"]))),
                "all_starts_rc0": str(all_rc0),
                "max_start_ll_gap": row["max_start_ll_gap"],
                "tol": row["tol"],
                "max_start_dlog10": row["max_start_dlog10"],
                "n_boundary": str(int(float(row["n_boundary"]))),
                "n_flat_eig": str(int(float(row["n_flat_eig"]))),
                "n_floor_under_marker": str(int(float(row["n_floor_under_marker"]))),
                "loss_floor_all": row["loss_floor_all"],
                "interior_condition": "" if name not in cond else repr(cond[name]),
                "label": label,
            }
        )

    OUT.mkdir(exist_ok=True)
    csv_path = OUT / "diagnostics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    caution_count = sum(row["label"] == "GEOMETRY-CAUTION" for row in rows)
    lines = [
        "# TVP-AIDS Phase D0: reliability-diagnostic feasibility results",
        "",
        "This is a dev-only retrospective analysis of committed Phase C/Q0 evidence.",
        "It adds no public diagnostic, inference claim, or shipped source change.",
        "",
        "## Outcome",
        "",
        f"All {len(rows)} endpoints were stable across the three starts tried, while ",
        f"all {caution_count} received `GEOMETRY-CAUTION` under the frozen descriptive rules. ",
        "The two facts coexist: optimizer endpoint agreement does not remove boundary/flat-geometry risk.",
        "The current public result can report only the fit return code and direct scale checks; ",
        "objective agreement and geometry need extra computation.",
        "",
        "## Endpoint table",
        "",
        "| endpoint | k | max start LL gap / tol | max dlog10 Q | boundary | flat eig | floor < 1.92 | all-floor loss | interior cond | label |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        ratio = num(row["max_start_ll_gap"]) / num(row["tol"])
        condition = None if not row["interior_condition"] else num(row["interior_condition"])
        lines.append(
            f"| {row['endpoint']} | {row['k']} | {fmt(ratio)} | "
            f"{fmt(num(row['max_start_dlog10']))} | {row['n_boundary']} | "
            f"{row['n_flat_eig']} | {row['n_floor_under_marker']} | "
            f"{fmt(num(row['loss_floor_all']))} | {fmt(condition)} | {row['label']} |"
        )
    lines += [
        "",
        "## Feasibility matrix",
        "",
        "| diagnostic | one public fit | extra starts | private objective work | recommendation |",
        "|---|---|---|---|---|",
        "| CMLMT completion | yes: `mleRetcode` | no | no | keep as the mandatory hard status |",
        "| Q positivity/finiteness and `Qfit/q0` scale | yes | no | no | cheap descriptive scale warning only; `q0` is caller-chosen |",
        "| Endpoint agreement | no | yes, plus comparable LL | objective exposure is currently private | opt-in multi-start reliability layer |",
        "| Boundary proximity against a defensible reference | only start-relative | sometimes | yes for Q0's truth/reference rule | do not call this identification from one fit |",
        "| Score/gradient or fixed-point residual | no | no | yes | not available in `quaidsTVPOut`; investigate separately |",
        "| Hessian rank/conditioning | no | no | yes, multiple evaluations | advanced diagnostic; boundary-aware interpretation required |",
        "| Profile/floor loss | no | no | yes, many reduced fits | research/advanced mode, not a default fit warning |",
        "",
        "## Findings",
        "",
        "1. A zero `mleRetcode` and agreement from 0.1x/1x/10x starts are useful but incomplete. Every endpoint passed both and every endpoint still showed at least one frozen geometry caution.",
        "2. The example-data anchor is the clearest case: all 12 Q elements are on Q0's boundary rule, 10 Hessian directions are flat, and the all-floor point loses only about 0.112 total LL, despite start agreement.",
        "3. A single-fit public warning can honestly cover only solver status and caller-start-relative scale collapse. It must not claim identification, standard errors, or global optimality.",
        "4. The next public-design candidate should be a separate result/procedure with explicit cost tiers: cheap single-fit checks, opt-in multi-start agreement, and advanced geometry diagnostics. Do not overload `mleRetcode`.",
        f"5. As a broader correctness cross-check, all {len(conv)} committed Phase C rows are `CONV-GLOBAL`. Phase S0 nevertheless found the available analytic gradient slower in 7 of 9 pilot cells, so advanced score/curvature work should not be assumed cheap.",
        "",
        "## Limits",
        "",
        "The set is five purposive endpoints, excludes `n1=4`, and reuses Q0's descriptive thresholds. It measures no false-positive/false-negative rate for warnings. The DGP reference Q is unavailable in applications, while a caller's `q0` is arbitrary. Therefore D0 supports an API feasibility decision, not calibrated warning thresholds for a release.",
        "",
        "Regenerate with `python dev/experiments/tvp_reliability_diagnostics/analyze_d0.py`.",
        "",
    ]
    (HERE / "RESULTS.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
