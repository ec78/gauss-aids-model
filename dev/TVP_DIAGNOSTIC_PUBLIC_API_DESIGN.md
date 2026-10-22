# Public TVP diagnostic API design

Status: **D2 design candidate — implementation not yet started.**

This design promotes the two evidence-backed diagnostic tiers from the private
D1/P6 work without promoting its underspecified edges. The implementation is a
new public API surface and therefore targets package version **0.5.0**.

## 1. Scope and claims

The public feature answers two distinct questions:

1. Given one completed `quaidsTVPFit()` result, did the optimizer report
   success, are fitted/start Q values usable, and did fitted Q become small
   relative to the caller's own start?
2. Across the fixed 0.1x, 1x and 10x Q starts, did three complete fits produce
   valid results and agree in observed total filter likelihood within
   `1e-6*nobs`?

It does **not** claim global optimality, a statistical boundary, parameter
identification, Q standard errors, confidence intervals or valid inference.
The scale marker is explicitly start-relative because `q0` is caller-chosen.
Likelihood agreement never requires Q-vector agreement: the private studies
show that materially different Q values can lie on a flat likelihood ridge.

Excluded from 0.5.0: custom start grids, Hessian/rank/profile diagnostics,
score thresholds, fixed-point residuals, automatic retry, parallel fitting,
new estimation math and changes to `quaidsTVPFit()` itself.

## 2. Public procedures

### 2.1 Cheap tier

```gauss
proc (struct quaidsTVPDiagOut) = quaidsTVPDiagnostics(
    struct quaidsTVPOut tvOut, scaleRatio=1e-3);
```

This procedure performs no fit and no Kalman-filter call. It only summarizes
fields already returned by `quaidsTVPFit()`.

### 2.2 Opt-in multi-start tier

```gauss
proc (struct quaidsTVPMultiStartOut) = quaidsTVPMultiStartDiagnostics(
    w, prices, totexp, H, q0, struct quaidsTVPControl tvpCtl,
    scaleRatio=1e-3);
```

This procedure performs three sequential, complete Q-MLE fits from
`0.1*q0`, `q0` and `10*q0`. Only `tvpCtl.maxIters` and `tvpCtl.maxTime` are
used; `smooth` and `othnam` are irrelevant because no state path or report is
returned. Worst-case optimizer time is approximately three times the per-pass
cap, plus the documented possible one-iteration overruns.

The fixed start grid and fixed likelihood tolerance are part of the 0.5.0
contract. They are not configurable until separate evidence supports a wider
surface.

No public print procedures are added initially. The result structs are compact,
and avoiding two more printer APIs keeps the first public contract small.

## 3. Public result structures

Both structs contain only GAUSS matrix/string fields and therefore add no
compile-time SSLib type dependency to `quaids.sdf`.

```gauss
struct quaidsTVPDiagOut {
    matrix status;
    string statusLabel;
    matrix mleRetcode;
    matrix fitOk;
    matrix qValid;
    matrix scaleRatio;
    matrix scaleCaution;
    matrix nScaleCollapse;
    matrix minQOverQ0;
    matrix maxQOverQ0;
};
```

Field shapes are scalar except `statusLabel` (string). Both structs use the
shared status vocabulary defined below; the cheap tier uses codes 0, 1, 2 and
5:

| code | label | meaning |
|---:|---|---|
| 0 | `OK` | converged, Q valid, no start-relative scale caution |
| 1 | `FIT-FAILED` | `mleRetcode != 0` |
| 2 | `INVALID-Q` | Qfit/q0 shape, finiteness or strict positivity failed |
| 5 | `SCALE-CAUTION` | one or more `Qfit/q0 < scaleRatio` |

Precedence is `FIT-FAILED`, `INVALID-Q`, `SCALE-CAUTION`, `OK`. Independent
fields preserve coexisting facts; for example, `qValid=0` remains visible when
the overall status is `FIT-FAILED`. When Q is invalid, `scaleCaution`,
`nScaleCollapse`, `minQOverQ0` and `maxQOverQ0` are missing rather than zero.

```gauss
struct quaidsTVPMultiStartOut {
    matrix status;
    string statusLabel;
    matrix scaleRatio;
    matrix multipliers;
    matrix retcodes;
    matrix fitOk;
    matrix qValid;
    matrix llValid;
    matrix scaleCaution;
    matrix anyScaleCaution;
    matrix loglikelihoods;
    matrix qfits;
    matrix nScaleCollapse;
    matrix bestStartIndex;
    matrix bestMultiplier;
    matrix bestLL;
    matrix maxLLGap;
    matrix tolerance;
    matrix startStable;
};
```

Exact shapes are:

- scalar: `status`, `scaleRatio`, `anyScaleCaution`, `bestStartIndex`,
  `bestMultiplier`, `bestLL`, `maxLLGap`, `tolerance`, `startStable`;
- string: `statusLabel`;
- `3x1` columns in fixed 0.1/1/10 order: `multipliers`, `retcodes`, `fitOk`,
  `qValid`, `llValid`, `scaleCaution`, `loglikelihoods`, `nScaleCollapse`;
- `k_states x 3`: `qfits`, where column i corresponds to multiplier i.

Both result structs share one numeric vocabulary so generic status handling
cannot silently reinterpret a code. The multi-start tier uses all six codes:

| code | label | meaning |
|---:|---|---|
| 0 | `START-STABLE-OBSERVED` | all results valid/converged, LL agreement passes, no scale caution |
| 1 | `FIT-FAILED` | one or more optimizer return codes are nonzero |
| 2 | `INVALID-Q` | one or more fitted Q vectors are unusable |
| 3 | `INVALID-LOGLIKELIHOOD` | one or more common-scale LL values are unusable |
| 4 | `START-SENSITIVE` | valid converged fits exceed the LL tolerance |
| 5 | `SCALE-CAUTION` | starts agree in LL, but at least one arm has start-relative collapse |

The shared meanings are 0 = tier-specific success, 1 = fit failed, 2 = invalid
Q, 3 = invalid likelihood, 4 = start sensitive and 5 = scale caution. Cheap
diagnostics never return codes 3 or 4.

Precedence is `FIT-FAILED`, `INVALID-Q`, `INVALID-LOGLIKELIHOOD`,
`START-SENSITIVE`, `SCALE-CAUTION`, `START-STABLE-OBSERVED`. Retaining the
mandatory nonzero-retcode warning as the top-level hard status is consistent
with `quaidsTVPFit()`, D0 and the cheap tier. `fitOk`, `qValid`, `llValid`,
`startStable`, `scaleCaution` and
`anyScaleCaution` preserve every independent fact. Thus the P6 pattern is
reported honestly as `SCALE-CAUTION` with `startStable=1`, not as a claim that
the warning disappeared.

`bestStartIndex`, `bestMultiplier`, `bestLL` and `maxLLGap` are missing whenever
any arm has invalid Q or LL. `bestStartIndex` is the 1-based row index into the
fixed multiplier vector, not a Q value. No missing/nonfinite value is passed to
`maxindc()`. An exact LL tie selects the first arm, matching GAUSS `maxindc()`
and the frozen private behavior.

## 4. Validation and numerical rules

### 4.1 `scaleRatio`

A shared private resolver validates `scaleRatio` before any expensive work. It
must be numeric, scalar, nonmissing, non-NaN, finite and strictly
`0 < scaleRatio < 1`. Invalid values fail loudly with a procedure-specific
message. The collapse comparison remains strict (`<`), so equality at the
threshold is not counted.

### 4.2 Q validity

The cheap tier requires `Qfit` and `q0` to be nonempty conformable column
vectors whose entries are finite and strictly positive.

The multi-start tier validates the base `q0` shape/finiteness/positivity before
the first fit, in addition to existing TVP model guards. It then constructs all
three derived starts and validates each as conformable, finite and strictly
positive before any optimizer call. This catches overflow or underflow caused
by scaling an otherwise-valid extreme `q0` and prevents a partially executed
diagnostic. After each fit it
validates the returned Q vector **before assigning it** to a preallocated
`k_states x 3` result. The result matrix begins all-missing. Only a conformable
`k_states x 1`, finite, strictly positive vector is stored. This prevents a
wrong-shaped returned vector from causing a GAUSS assignment error. If an arm's
Q is invalid, the implementation does not call the Kalman filter for that arm;
its Q column, LL, collapse count and scale flag remain missing, while
`qValid=0` and `llValid=0`.

For a valid-Q arm, collapse is exactly
`Qfit_i < scaleRatio*(multipliers[i]*q0)`. Its collapse count and 0/1 scale flag
are populated even if the later LL is invalid. `anyScaleCaution=1` if any valid
arm cautions; it is 0 only when all three Q vectors are valid and none cautions;
otherwise it is missing.

### 4.3 Likelihood validity and comparison

Each valid-Q arm is re-evaluated through the production fixed-Q Kalman-filter
path, and its total LL is the sum of per-observation filter log likelihoods.
An LL is valid only if it is a finite, nonmissing scalar. A GAUSS/filter engine
error can still abort the call; status code 3 covers returned invalid numerical
values, not arbitrary engine exceptions that GAUSS cannot reliably trap. An
invalid returned LL is normalized to missing; its raw invalid value is not
stored in the public result.

When all three LLs are valid:

```text
tolerance = 1e-6*nobs
maxLLGap = max(bestLL - loglikelihoods)
startStable = all(retcodes == 0) and maxLLGap <= tolerance
```

The tolerance is inclusive. Nonzero-retcode arms may retain descriptive valid
Q/LL values and comparison fields, but `startStable` is zero and overall status
is `FIT-FAILED`. `startStable` is always a 0/1 scalar and is zero whenever any
retcode is nonzero, Q/LL validity is incomplete, or the LL gap exceeds the
tolerance; it is never missing.

## 5. Source and optional dependency placement

Implementation belongs at the end of existing optional adapter
`src/quaidstvpfit.src`, after the current public fit/print procedures. This is
the smallest loading change: that file is already last in the documented TVP
include chain and already has every private dependency required by the
multi-start tier.

The new structs belong in `src/quaids.sdf`. `package.json`'s `src` and `deps`
arrays remain unchanged: diagnostics are part of the existing optional TVP
module and must not make SSLib a core-package dependency. The documented setup
remains:

```gauss
library cmlmt, tsmt, sslib, quaids;
#include quaids.sdf
#include quaidstvp.src
#include quaidstvpkalman.src
#include quaidstvpmle.src
#include quaidstvpsmooth.src
#include quaidstvpelas.src
#include quaidstvpfit.src
```

No new source filename means no new intentionally-unlisted manifest entry.

## 6. Documentation and compatibility

Version 0.5.0 changes together:

- `package.json`, `docs/public-api.json`, `CITATION.cff` and release metadata;
- `CHANGELOG.md`;
- `docs/COMMAND_REFERENCE.md` plus command pages for both procedures;
- `README.md`, `docs/USAGE_GUIDE.md`, `docs/FEATURE_SUPPORT_MATRIX.md` and
  relevant TVP cross-links;
- `docs/public-api.json` supported procedures, TVP optional-module procedures,
  and public struct names.

Example 14 should demonstrate only the cheap call after its existing fit. The
three-fit diagnostic is documented with a bounded-code snippet but is not added
to routine example smoke execution because its cost is material.

This is additive relative to 0.4.0. Existing structs/procedures and
`quaidsTVPFit()` behavior are unchanged. Once released, names, fields, codes,
precedence and fixed-grid semantics become public compatibility commitments.

## 7. Implementation and verification phases

### D3 — semantics freeze and implementation

Files: `src/quaids.sdf`, `src/quaidstvpfit.src`, a focused new test, guard
scripts and a frozen implementation protocol under `dev/`.

Acceptance:

- pure classifier tests cover every precedence branch;
- scale threshold is strict and LL tolerance inclusive;
- exact ties choose the first arm;
- missing/NaN/Inf Q and LL are never passed to summary operations;
- bad scale ratios (zero, one, negative, missing, NaN, Inf, nonscalar) fail
  before a fit for both public procedures;
- malformed base q0 (empty, wrong rows/columns, missing, NaN, Inf, zero or
  negative) fails before the first fit;
- a derived 0.1x, 1x or 10x start that becomes nonfinite or nonpositive through
  scaling fails before any optimizer call;
- wrong-shaped/invalid returned Q is neither assigned nor filtered, invalid LL
  is normalized to missing, and every overlapping retcode/Q/LL precedence case
  is covered;
- `SCALE-CAUTION` code 5 can coexist with `startStable=1`, and each arm's
  collapse comparison uses its own scaled start;
- the Stage-6 fixture's 1x arm exactly matches `quaidsTVPFit()` in Q, retcode
  and common-scale LL;
- private D1 files remain byte-identical.

### D4 — public integration and documentation

Update public API metadata, documentation, version/changelog/release metadata
and the installed optional-TVP smoke test. Keep the existing include chain and
confirm the core package still loads without SSLib.

Acceptance:

- public API, manifest, docs consistency and docs quality checks pass;
- focused source and installed-artifact tests pass;
- example smoke remains bounded;
- the new API reproduces the frozen D1 n1=2 evidence.

### D5 — bounded evidence and release gate

Run the public implementation against the frozen P6 anchor outside routine CI,
then run the full 0.5.0 release gate.

Acceptance:

- P6 retains valid rc/LL/Q results, observed-start agreement and explicit
  `SCALE-CAUTION` with `startStable=1`;
- strict SSLib 2.1.2 pin passes;
- full source suite, artifact build/install verification, installed public API,
  examples and release checksum record all pass;
- an independent review finds no overclaim or optional-dependency regression.

## 8. Decision gate

Approve this contract before D3 modifies shipped source. The main owner choice
is whether the multi-start overall status should surface scale caution as code 5
while preserving `startStable=1` (this design), or retain the private prototype's
code 0 and expose scale caution only in secondary fields. This design recommends
shared code 5 because a public diagnostic should make the warning difficult to
miss.
