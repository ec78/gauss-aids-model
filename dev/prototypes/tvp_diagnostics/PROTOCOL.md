# TVP-AIDS Phase D1: private tiered diagnostic prototype

Status: **FROZEN v1 before live prototype measurement.**

## Scope

D1 implements a private, dev-only feasibility prototype following Phase D0.
Nothing under this directory is shipped, loaded by `library quaids`, listed in
`package.json`, or public API.

The prototype has two deliberately separate cost tiers:

1. `_quaidsTVPDiagSingle(tvOut, scaleRatio)` reads one existing
   `quaidsTVPOut`. It checks `mleRetcode`, Q/q0 shape and finite positivity,
   and counts fitted Q elements below `scaleRatio*q0`. The default candidate
   marker is `1e-3`. Because q0 is caller-chosen, this is named
   **start-relative scale collapse**, never a boundary or identification test.
2. `_quaidsTVPDiagMultiStart(...)` performs the fixed 0.1x, 1x and 10x q0
   starts, records each return code and total filter log-likelihood, and calls
   the starts stable only when every retcode is zero and every
   `bestLL - LL_i <= 1e-6*nobs`. Q-vector equality does not gate stability:
   equivalent likelihoods may lie along a flat boundary ridge.

Status precedence is `FIT-FAILED`, `INVALID-Q`, `SCALE-CAUTION`, `OK` for the
single-fit tier and `FIT-FAILED`, `START-SENSITIVE`, `START-STABLE-OBSERVED`
for the multi-start tier. The last label is not a global-optimum claim.

## Explicit exclusions

- no Hessian, rank, profile, score-threshold or fixed-point diagnostic;
- no Q standard errors, intervals, chi-square reference or identification
  claim;
- no public struct/procedure, source, docs, example, package or version edit;
- no SSLib change and no broad grid.

## Frozen checks

- `tvpdiag_unit_test.e`: synthetic structure and pure-summary tests covering
  every precedence branch, the inclusive likelihood tolerance, invalid Q,
  and the start-relative threshold boundary;
- `tvpdiag_equiv.e`: the Stage-6 n1=2 fixture, requiring the multi-start 1x
  column to equal `quaidsTVPFit()` exactly in Qfit and retcode and to use the
  same total filter likelihood;
- at most two descriptive live probes after the freeze: Phase C/Q0's
  n1=2/tobs=200/seed=1001/oracle-H endpoint and the P6 n1=3/tobs=500/seed=204
  static-H anchor.

All live runs use the repository's GAUSS config override and pinned SSLib.
Protocol, prototype, includes, unit/equivalence files and runner are hashed in
`FREEZE.sha256` and committed before equivalence/probe measurement. A later
public proposal requires a separate owner decision.
