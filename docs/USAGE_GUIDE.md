# QUAIDS Usage Guide

This guide summarizes the main API choices, model-selection switches, and
output conventions for the GAUSS QUAIDS package.

## Choosing An API

Use [quaidsPreflight](command-reference/quaidsPreflight.md) before fitting
when you want a silent data/design screen for common applied problems:
dimension mismatches, non-finite values, share adding-up failures,
zero/negative shares, weak instruments, low variation, cluster counts, and
basic convergence-risk hints:

```gauss
pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, clusterId=householdId);

if not pOut.ok;
    call printQuaidsPreflight(pOut);
    end;
endif;
```

Use [quaidsFit](command-reference/quaidsFit.md) when you want a silent,
struct-returning call with no console output -- the right choice for
scripts, simulations, and anywhere you only need the returned structure:

```gauss
aCtl = quaidsControlCreate();

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

Use [quaids](command-reference/quaids.md) (the original, backward-compatible
call) when you want the full console report -- estimation tables,
elasticities at the mean/quartiles, descriptive statistics, and the
Slutzky diagnostic -- printed automatically:

```gauss
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);
```

Use [quaidsFull](command-reference/quaidsFull.md) instead of
[quaidsFit](command-reference/quaidsFit.md) when your data is already a
GAUSS dataframe and you'd rather select columns by name than assemble
matrices by hand:

```gauss
data = loadd("mydata.csv");
qOut = quaidsFull(data, shareVars, priceVars, "totexp", "instr", extraVars, aCtl);
```

Use [quaidsWorkflowFit](command-reference/quaidsWorkflowFit.md) when you
want the first applied workflow bundle: run the compact preflight summary,
estimate the system, evaluate predicted shares and elasticities at the
sample mean, and compute robust or cluster-robust standard errors in one
silent, struct-returning call:

```gauss

// clusterId omitted defaults to 0 -- heteroskedasticity-robust SE.
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);

// Or pass a Tx1 group-label vector for cluster-robust SE.
wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl, clusterId=householdId);
```

The workflow object is intentionally a composition layer, not a separate
estimator. Its fit fields come from `quaidsFit()`, and its post-estimation
fields match direct calls to `quaidsSharesFit()`, `quaidsElasFit()`, and
`quaidsRobustFit()` on the same sample-mean evaluation point, including
robust/cluster-robust standard errors propagated into shares and
elasticities via `quaidsRobustCovariance()`. It also includes
model/restriction summary fields such as `symPval`, `overidPvf`, and,
for unconstrained QUAIDS fits, `quadraticPval`. Workflow preflight fields
such as `preflightOk`, `preflightWarnings`, `preflightIVFstat`, and
`preflightNClusters` are diagnostic only; call `quaidsPreflight()` first
and stop on `pOut.ok == 0` when you want a hard pre-estimation guard.

Use [quaidsWorkflowScenarioFit](command-reference/quaidsWorkflowScenarioFit.md)
when the same workflow should also return exact CV/EV for a price-change
scenario:

```gauss
pricesPt1 = wfOut.evalPrices;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

wfScenario = quaidsWorkflowScenarioFit(w, intcpt, prices, totexp, instr, aCtl,
    wfOut.evalIntcpt, wfOut.evalPrices, pricesPt1, wfOut.evalTotexp);

if wfScenario.welfareValid;
    print wfScenario.cv wfScenario.seCV;
    print wfScenario.ev wfScenario.seEV;
endif;

if wfScenario.welfareRobustValid;
    print wfScenario.cv wfScenario.seCVRobust;
    print wfScenario.ev wfScenario.seEVRobust;
endif;
```

Use [quaidsSurveyWorkflowFit](command-reference/quaidsSurveyWorkflowFit.md)
when a household survey or other microdata source has sampling weights and
you want both a genuinely weighted estimator and population-representative
post-estimation summaries:

```gauss
wfSurvey = quaidsSurveyWorkflowFit(w, intcpt, prices, totexp, instr, aCtl,
    sampwt, clusterId=householdId);

print wfSurvey.weighted wfSurvey.weightSum wfSurvey.effN;
print wfSurvey.surveyWeightSum;
print wfSurvey.evalTotexp;
print wfSurvey.shares;
print wfSurvey.incomeElas;
```

Since Milestone 26, `sampwt` does double duty: it both fits `quaidsFit()`
as a genuine sampling-weight-adjusted estimator (via
[quaidsWorkflowFit](command-reference/quaidsWorkflowFit.md)'s own optional
`weight` argument) and changes the evaluation point used for predicted
shares and elasticities, with classical and robust/cluster-robust
delta-method standard errors recomputed at that point using the weighted
fit. If you only want the estimator-level weighting without also
reweighting the evaluation point, call
[quaidsWorkflowFit](command-reference/quaidsWorkflowFit.md) directly with
its own optional `weight` argument instead. Weighted point estimates and a
matching weighted/clustered sandwich SE are the current scope; formal
strata as a concept distinct from clustering, replicate-weight
(BRR/jackknife) variance, and finite-population correction remain roadmap
items.

There is no formula-string (`"y ~ x1 + x2"`) API -- AIDS/QUAIDS is a
multi-equation system (N budget shares against N parallel log prices),
which doesn't fit GAUSS's single-equation formula grammar. Column-name
lists, matched positionally between `shareVars` and `priceVars`, are the
natural fit instead.

## Choosing A Model: LA-AIDS vs. Iterated AIDS vs. QUAIDS

All three are the same estimator, `quaidsFit`, selected by
`aCtl.linear` and `aCtl.maxiter`:

| Model | `aCtl.linear` | `aCtl.maxiter` | Price index |
| --- | --- | --- | --- |
| LA-AIDS | either | `1` | Stone (linear approximation) |
| Iterated AIDS | `1` | `> 1` | Nonlinear translog, iterated |
| QUAIDS | `0` | `> 1` | Nonlinear translog, iterated, plus a quadratic log-expenditure term |

```gauss
aCtl = quaidsControlCreate();

// LA-AIDS: one-step, Stone price index.
aCtl.maxiter = 1;

// Iterated (linear) AIDS: nonlinear translog price index, no quadratic term.
aCtl.linear = 1;
aCtl.maxiter = 100;

// QUAIDS: nonlinear translog price index, plus the quadratic term.
aCtl.linear = 0;
aCtl.maxiter = 100;
```

`aCtl.maxiter == 1` always uses the Stone index regardless of
`aCtl.linear`'s value -- LA-AIDS never iterates past the starting value.
For `aCtl.maxiter > 1`, `qOut.model` reports which of `"AIDS"`/`"QUAIDS"`
was actually fit.

**Convergence is not guaranteed** for the iterated estimator. A real,
committed 200-seed sweep (`tests/quaids_convergence_sweep.e`, Milestone
12 -- run it yourself with `tests/run_convergence_sweep.ps1`) measured,
at default settings: iterated AIDS fails (never converges, or converges
to a self-consistent but wrong answer) 58% of the time; QUAIDS 76%. An
optional damping control can help:

```gauss
aCtl.relax = .75;   // default is 1 (no damping); .75 measurably reduced
                    // the failure rate in the Milestone 12 sweep
```

`aCtl.relax` under-relaxes the fixed-point update (`b_new = relax*b +
(1-relax)*b_old` each iteration); more aggressive damping (`.5`, `.3`)
did not help further in testing and often made things worse. This is a
modest, evidence-backed mitigation, not a guarantee -- always check
`qOut.converged` and `qOut.iterations` after any iterated fit before
trusting the result. See
[Feature Support Matrix](FEATURE_SUPPORT_MATRIX.md#notes) and
`GOLD_STANDARD_TODO.md`'s Milestone 12 section for the full breakdown.

## Instrumental Variables Are Always Required

`quaidsFit()`, `quaids()`, and `quaidsFull()` always treat log total
expenditure as endogenous and instrument it via a control-function
approach -- `instr` is a required argument, not optional, and there is no
"no IV" estimation mode. Provide at least one instrument column; provide
more than one (`ninst > nu`, where `nu` is the number of endogenous
total-expenditure regressors) to activate the overidentification test
(`qOut.overidValid`/`qOut.overidGamma`/`qOut.overidFstat`/`qOut.overidPvf`).

## Homogeneity, Symmetry, and Overidentification

Set `aCtl.homogenous = 1` (the default) to impose homogeneity via
minimum-distance FGLS and additionally test/report symmetry given
homogeneity (`qOut.symStat`/`qOut.symPval`) and fit a
symmetry-constrained system (`qOut.symcB`/`qOut.symcV`). Set
`aCtl.homogenous = 0` to fit unconstrained instead -- required before
calling the standalone Wald tests:

```gauss
aCtl = quaidsControlCreate();
aCtl.homogenous = 0;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ statH, pvalH, dfH } = quaidsHomogeneityTest(qOut);
{ statJ, pvalJ, dfJ } = quaidsJointTest(qOut);
```

Both [quaidsHomogeneityTest](command-reference/quaidsHomogeneityTest.md)
and [quaidsJointTest](command-reference/quaidsJointTest.md) error clearly
if called on a `qOut.homogenous == 1` fit.

## Is QUAIDS Needed? (`quaidsQuadraticTest`)

If `aCtl.linear = 0` (QUAIDS) and `aCtl.homogenous = 0` (same
unconstrained-fit requirement as above),
[quaidsQuadraticTest](command-reference/quaidsQuadraticTest.md) is a Wald
test of whether the quadratic log-expenditure term is actually needed --
a failure to reject means plain AIDS (`aCtl.linear = 1`, simpler, more
numerically stable, see "Choosing A Model" above) fits just as well:

```gauss
{ statQ, pvalQ, dfQ } = quaidsQuadraticTest(qOut);
if pvalQ > 0.05;
    print "Quadratic term not statistically justified -- consider AIDS.";
endif;
```

Only callable on a QUAIDS fit -- an AIDS fit (`aCtl.linear = 1`) never
estimates the quadratic term at all, so there is nothing to test.

## Elasticities At Any Point

`quaidsElas_()` (the low-level computation) always accepted an arbitrary
evaluation point -- the Milestone 5 generalization was giving that a
silent, struct-returning entry point:

```gauss
n = qOut.n;
nint = qOut.nint;

// At the sample mean:
m_ = meanc(qOut.intcptFull~prices~totexp~w~instr);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call printQuaidsElas(elasOut);

// At a synthetic counterfactual -- e.g. a hypothetical 20% price increase
// on good 1, evaluated at the sample mean otherwise:
pricesCf = pricesMean;
pricesCf[1] = pricesCf[1] + ln(1.20);
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesCf, totexpMean, aCtl);
```

Always evaluate elasticities/[quaidsSlutzky](command-reference/quaidsSlutzky.md)
against `qOut.bestB`/`qOut.bestV` -- "whichever is the most-constrained
estimate actually fit" (symmetric if homogeneity was imposed, else the
recovered unconstrained fit).

## Predicted Budget Shares At Any Point

[quaidsSharesFit](command-reference/quaidsSharesFit.md) (Milestone 16)
exposes the same model-implied share `quaidsElasFit`'s elasticities are
built on, directly, at any evaluation point -- useful for out-of-sample
prediction and policy simulation without hand-deriving the share
equation:

```gauss
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
call printQuaidsShares(sharesOut);

print "sum of predicted shares:" sumc(sharesOut.w);   // exactly 1, by construction
```

Reports budget shares only (`w_i`, the fraction of expenditure spent on
good `i`), not physical quantities demanded -- converting to quantities
would require assuming `prices`/`totexp` are logs of levels in mutually
consistent units, which nothing else in this library requires. `sharesOut.v`
is the full covariance of `w` (not just marginal standard errors), so
hypotheses spanning more than one good -- e.g. whether two goods' shares
differ significantly -- can be tested with the correct combined variance.

## Welfare Analysis (Compensating/Equivalent Variation)

[quaidsWelfareFit](command-reference/quaidsWelfareFit.md) computes exact
CV/EV for a hypothetical price change, holding nominal expenditure fixed
-- works for any model choice (LA-AIDS, iterated AIDS, QUAIDS), no extra
package required:

```gauss
n = qOut.n;
nint = qOut.nint;
intcptPt = meanc(qOut.intcptFull);
pricesPt0 = meanc(prices);
totexpPt0 = meanc(totexp);

// A hypothetical 5% price increase on good 1, all else unchanged:
pricesPt1 = pricesPt0;
pricesPt1[1] = pricesPt1[1] + ln(1.05);

wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
call printQuaidsWelfare(wOut);
```

`wOut.cv`/`wOut.ev` are positive when the price change reduces welfare,
negative when it improves welfare, and exactly zero when
`pricesPt1 == pricesPt0`. See
[Methodology Notes](METHODOLOGY_NOTES.md#welfare-measures) for the exact
formula and how it was verified.

## Imposing Curvature (Diewert-Wales)

`quaidsSlutzky()` always diagnoses curvature (Slutzky negative
semidefiniteness) but never imposes it.
[quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md) can impose
it locally, at the sample mean, for LA-AIDS/AIDS (`aCtl.linear = 1`) and,
since Milestone 13, QUAIDS (`aCtl.linear = 0`) too -- requires the
`optmt` package either way:

```gauss
library optmt, quaids;

aCtl = quaidsControlCreate();
aCtl.linear = 1;         // or 0 for QUAIDS
aCtl.maxiter = 100;
aCtl.homogenous = 1;    // required -- quaidsCurvatureFit needs a
                        // homogeneity+symmetry-constrained starting fit

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

aCtl.relax = .25;    // recommended for QUAIDS -- its curvature outer
                     // loop is measurably less stable than AIDS's own;
                     // undamped (relax=1, the default) runs can diverge.
                     // Not needed for AIDS.

cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
call printQuaidsCurvature(cOut);

print "Slutzky eigenvalues at the sample mean:" cOut.eigenvalues';  // all <= 0
```

For QUAIDS, `cOut.b`/`cOut.se` gain a trailing `lambda` row (no `u` row
either way -- see [quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md)).
See the [Limitations section](#limitations) below for the standard-error
caveat, which applies to both models.

## Zero Budget Shares (Corner Solutions)

Real survey/microdata routinely has corner solutions -- some households
report zero expenditure on some goods. Fitting [quaidsFit](command-reference/quaidsFit.md)
directly on such data is a known source of bias, since the linear/log-
linear share equation has no mechanism for a censored dependent variable.
[quaidsZeroFit](command-reference/quaidsZeroFit.md) corrects for this via
the Shonkwiler & Yen (1999) two-step procedure -- a per-good first-stage
probit (the probability of a non-zero share) followed by a corrected
second-stage GLS fit:

```gauss
aCtl = quaidsControlCreate();
aCtl.linear = 0;          // 1 for AIDS, 0 for QUAIDS
aCtl.maxiter = 100;
aCtl.homogenous = 1;      // 1 to impose homogeneity/symmetry; 0 unconstrained

zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOut);

print "fraction of zero shares per good:" zOut.shareZeroFrac';
```

`zOut.b`'s trailing `n x n` block (`delta`) is the estimated own-good
hazard coefficient, restricted to be diagonal (each good's correction only
depends on its own censoring probability) by construction. See
[Methodology Notes](METHODOLOGY_NOTES.md#zero-budget-share-correction-shonkwiler-yen)
for the full derivation and the [Limitations section](#limitations) below
for what is deliberately out of scope in this first pass (no homogeneity/
symmetry imposition, a simplified standard-error formula, and adding-up
not holding exactly for the corrected coefficients -- a real property of
the method itself, not a bug).

## Robust and Cluster-Robust Standard Errors

[quaidsRobustFit](command-reference/quaidsRobustFit.md) generalizes every
other covariance in this library (which are all built on a pooled,
homoskedastic `S = sse/nobs`) to a heteroskedasticity-robust or
cluster-robust sandwich, given an already-fitted `qOut`:

```gauss
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

// Heteroskedasticity-robust (clusterId omitted -- defaults to 0):
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl);
call printQuaidsRobust(rOut);

// Cluster-robust (householdId is a Tx1 vector of group labels):
rOutCluster = quaidsRobustFit(qOut, w, prices, totexp, aCtl, clusterId=householdId);
call printQuaidsRobust(rOutCluster);
```

Robust and cluster-robust are the same formula -- `clusterId = 0` means
heteroskedasticity-robust (every observation is its own cluster); a
supplied `Tx1` vector of group labels means cluster-robust with a CR1
small-sample correction.

A cluster-aware bootstrap is also available, resampling whole clusters
(or plain rows, when `clusterId` is unset) and refitting
[quaidsFit](command-reference/quaidsFit.md) on each resample:

```gauss
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, clusterId=householdId, seed=42);
call printQuaidsRobustBootstrap(rbOut);
```

For shares, elasticities, or welfare, expand the reduced robust covariance
back into `qOut.bestB`'s full coefficient basis before calling the
delta-method post-estimation procs:

```gauss
{ bR, vR } = quaidsRobustCovariance(qOut, rOutCluster, aCtl);

sharesR = quaidsSharesFit(bR, vR, intcptPt, pricesPt, totexpPt, aCtl);
elasR = quaidsElasFit(bR, vR, intcptPt, pricesPt, totexpPt, aCtl);
welfareR = quaidsWelfareFit(bR, vR, intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);

{ bB, vB } = quaidsRobustBootstrapCovariance(qOut, rbOut, aCtl);
sharesB = quaidsSharesFit(bB, vB, intcptPt, pricesPt, totexpPt, aCtl);
```

**A real, important caveat, found empirically**:
[quaidsRobustFit](command-reference/quaidsRobustFit.md)'s closed-form
sandwich uses a *simplified* bread (`inv(gg)`-based, not the full
nonlinear-price-index-feedback Jacobian correction
[quaidsFit](command-reference/quaidsFit.md) itself uses), which makes its
`se` dramatically more *conservative* (often 10-100x larger) than
`qOut.homogSE`/`symcSE` -- confirmed to be an expected consequence of
comparing a simple equation-by-equation sandwich against the full,
cross-equation-efficient FGLS system, not a bug in the formula itself
(verified against an independent hand-derivation using the same
regressors/residuals). `quaidsRobustBootstrapFit`'s `seBoot`, which
resamples the *actual* estimator, is typically much closer to `qOut`'s
own SE -- prefer it when this gap matters. See
[Methodology Notes](METHODOLOGY_NOTES.md#robust-and-cluster-robust-standard-errors)
for the full derivation.

The reduced coefficient table only covers the `n1` independently-estimated
equations (equation `n` is recovered via adding-up, matching every other
diagnostic in this library). Use `quaidsRobustCovariance()` or
`quaidsRobustBootstrapCovariance()` for the full-basis covariance needed
by shares, elasticities, and welfare.

## Replicate-Weight (Jackknife/BRR) Standard Errors

Use [quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md)
when your survey extract ships **pre-computed replicate weight columns**
(a common design for household-expenditure surveys) rather than requiring
you to implement your own resampling scheme:

```gauss
// replicateWeights: TxR matrix, one alternate Tx1 weight column per
// replicate, supplied by your survey's own documentation. scaleFactorJK1
// is that design's own prescribed scale factor -- e.g. (R-1)/R for a
// standard JK1 delete-one-PSU jackknife.
rOut = quaidsReplicateWeightFit(w, intcpt, prices, totexp, instr, aCtl,
    replicateWeights, scaleFactorJK1, weight=surveyWeight, method="JK1");

call printQuaidsReplicateWeight(rOut);
print "completed:" rOut.nCompleted "failed:" rOut.nFailed;
```

This is a genuinely different inference approach from
[quaidsRobustFit](command-reference/quaidsRobustFit.md)'s closed-form
sandwich or [quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md)'s
own resampling bootstrap -- it refits
[quaidsFit](command-reference/quaidsFit.md) once per **caller-supplied**
replicate weight column (not a random resample), and combines the results
using the general linearized replication-variance formula
`V = sum_r c_r * vec(b_r - b_full) * vec(b_r - b_full)'`. This proc does
not implement or auto-detect jackknife, BRR, or any other specific
design -- `replicateWeights` and `scaleFactor` are always required,
caller-supplied inputs, matching this library's own established "never
silently guess an inference-affecting parameter" discipline (the same
rule `quaidsCurvatureBootstrapFit`'s own `B` already follows).

A real, practical advantage over `quaidsRobustFit()`: `rOut.b`/`rOut.v`
are already in `quaidsFit()`'s own full `bestB` basis, so they feed
directly into `quaidsSharesFit()`/`quaidsElasFit()`/`quaidsWelfareFit()`
with no separate expansion step (`quaidsRobustCovariance()` has no
counterpart needed here). See
[Methodology Notes](METHODOLOGY_NOTES.md#replicate-weight-jackknifebrr-variance-estimation)
for the full derivation, including a real, non-trappable crash mode found
and guarded against while building this (a replicate weight concentrated
on too few effectively-weighted rows).

## Reporting (`pubtable`)

`src/pubtable_quaids.src` is an optional adapter onto the `pubtable`
package for LaTeX/Markdown/CSV/RTF/HTML/XLSX table export. It is **not**
loaded by `library quaids;` -- `#include` it directly after both `quaids`
and `pubtable` are available:

```gauss
library pubtable, quaids;
#include quaids.sdf
#include pubtable_quaids.src

coefTbl = ptFromQuaids(qOut);
call ptExport(coefTbl, "results.tex");

elasTbls = ptTablesFromQuaidsElas(elasOut);   // 3x1: income, uncompensated, compensated
call ptExport(elasTbls[1], "income_elasticities.md");

workflowTbls = ptTablesFromQuaidsWorkflow(wfOut);  // shares + elasticity tables, plus welfare if present
call ptExport(workflowTbls[1], "workflow_shares.md");
```

See the [command reference](COMMAND_REFERENCE.md#reporting-optional-requires-pubtable)
for each adapter proc, and `examples/pubtable_export_example.e` for a full
runnable example.

## Limitations

- Curvature **imposition** ([quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md))
  is available for LA-AIDS/AIDS and QUAIDS, at the sample mean, and its
  standard errors are a simplified delta-method approximation that is
  known to be unreliable when the estimated Cholesky factor has boundary
  (near-zero) entries -- see the [Methodology Notes](METHODOLOGY_NOTES.md#curvature-imposition-diewert-wales)
  for why this happens and why point estimates and the exact curvature
  property are unaffected. [quaidsCurvatureBootstrapFit](command-reference/quaidsCurvatureBootstrapFit.md)
  (Milestone 15) offers a bootstrap alternative that does not share this
  weakness, reported alongside (not replacing) the delta-method SE --
  see [printQuaidsCurvatureBootstrap](command-reference/printQuaidsCurvatureBootstrap.md).
  It has no default replication count: a single AIDS curvature fit takes
  under a second, but a single QUAIDS curvature fit takes several seconds,
  so a conventional replication count (200-1000) can mean anywhere from a
  few minutes to a couple of hours depending on the model -- choose `B`
  deliberately. For QUAIDS, `aCtl.relax` is effectively required (its
  curvature outer loop is measurably less stable than AIDS's own), and
  there is no known-curvature-consistent synthetic fixture to validate
  true-parameter recovery against (unlike AIDS) --
  `tests/quaids_curvature_test.e`'s QUAIDS checks validate convergence/
  exact NSD/shape instead, a real but weaker tier of evidence.
  [quaidsCurvatureBootstrapCI](command-reference/quaidsCurvatureBootstrapCI.md)
  (Milestone 18) computes percentile confidence intervals directly from
  `quaidsCurvatureBootstrapFit()`'s raw draws (`bootOut.bBoot`) -- no new
  resampling needed, but intervals from a small `B` are correspondingly
  crude.
- No guaranteed convergence for the iterated estimator (or the curvature-
  constrained outer iteration built on top of it) -- see "Choosing A
  Model" above. `aCtl.relax` (Milestone 12) is an evidence-backed, opt-in
  mitigation, not a fix.
- IV is mandatory; there is no exogenous-total-expenditure estimation mode.
- [quaidsFit](command-reference/quaidsFit.md) (Milestone 26) accepts an
  optional sampling-weight argument -- a genuine weighted point estimate,
  with a matching weighted/clustered sandwich SE via
  [quaidsRobustFit](command-reference/quaidsRobustFit.md), or
  replicate-weight (jackknife/BRR-style) SE via
  [quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md)
  (Milestone 27, always against caller-supplied replicate columns and
  scale factor -- no design is auto-detected). Formal strata as a concept
  distinct from clustering, and finite-population correction, are not
  implemented yet.
  [quaidsSurveyWorkflowFit](command-reference/quaidsSurveyWorkflowFit.md)
  wires the same base weight into both the estimator and the workflow's
  representative evaluation point.
- [quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md)
  (Milestone 27) does not retry a replicate that fails to converge (fixed,
  caller-supplied columns cannot be redrawn) -- a failed replicate is
  simply dropped from the variance sum, a documented simplification since
  the formal JK1/BRR literature does not define a missing-replicate
  adjustment this library implements. A replicate whose effective sample
  size falls below a defensive `2x`-design-columns heuristic is skipped
  before ever calling `quaidsFit()`, avoiding a real, non-trappable
  `error G0058` crash mode found while building this milestone -- see the
  [Methodology Notes](METHODOLOGY_NOTES.md#replicate-weight-jackknifebrr-variance-estimation).
- [quaidsZeroFit](command-reference/quaidsZeroFit.md) (Milestones 19/30)
  supports unconstrained estimation and homogeneity/symmetry imposition via
  `aCtl.homogenous`, but reports a simplified standard error that does not
  account for the nonlinear
  translog-price-index feedback or first-stage probit/IV generated-
  regressor uncertainty. Adding-up does not hold exactly for its corrected
  coefficients, a real property of the Shonkwiler-Yen method itself. Its
  first-stage probit relies on GAUSS's built-in `glm()`, which can
  hard-crash (not just fail to converge) on some degenerate inputs -- a
  known, non-trappable failure mode, not hardened against in this pass.
- [quaidsRobustFit](command-reference/quaidsRobustFit.md) (Milestone 20)
  uses a simplified bread, making its closed-form `se` dramatically more
  conservative than `qOut`'s own classical SE (a confirmed, expected
  property, not a bug -- see that proc's own Remarks).
  [quaidsRobustBootstrapFit](command-reference/quaidsRobustBootstrapFit.md)'s
  bootstrap SE does not share this weakness. Neither proc's covariance
  propagates automatically into `qOut.symcV` or into elasticities/shares/
  welfare's own delta-method SEs.
