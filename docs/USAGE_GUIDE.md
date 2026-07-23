# QUAIDS Usage Guide

This guide summarizes the main API choices, model-selection switches, and
output conventions for the GAUSS QUAIDS package.

## Choosing An API

Use [quaidsFit](command-reference/quaidsFit.md) when you want a silent,
struct-returning call with no console output -- the right choice for
scripts, simulations, and anywhere you only need the returned structure:

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();

struct quaidsOut qOut;
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
struct quaidsControl aCtl;
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
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.homogenous = 0;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ statH, pvalH, dfH } = quaidsHomogeneityTest(qOut);
{ statJ, pvalJ, dfJ } = quaidsJointTest(qOut);
```

Both [quaidsHomogeneityTest](command-reference/quaidsHomogeneityTest.md)
and [quaidsJointTest](command-reference/quaidsJointTest.md) error clearly
if called on a `qOut.homogenous == 1` fit.

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

struct quaidsElasOut elasOut;
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

struct quaidsWelfareOut wOut;
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
semidefiniteness) but never imposes it. For LA-AIDS/AIDS
(`aCtl.linear = 1`), [quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md)
can impose it locally, at the sample mean, requiring the `optmt` package:

```gauss
library optmt, quaids;

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;    // required -- quaidsCurvatureFit needs a
                        // homogeneity+symmetry-constrained starting fit

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
call printQuaidsCurvature(cOut);

print "Slutzky eigenvalues at the sample mean:" cOut.eigenvalues';  // all <= 0
```

QUAIDS (`aCtl.linear = 0`) is not yet supported -- `quaidsCurvatureFit`
errors clearly if called on a QUAIDS fit. See the
[Limitations section](#limitations) below for the standard-error caveat.

## Reporting (`pubtable`)

`src/pubtable_quaids.src` is an optional adapter onto the `pubtable`
package for LaTeX/Markdown/CSV/RTF/HTML/XLSX table export. It is **not**
loaded by `library quaids;` -- `#include` it directly after both `quaids`
and `pubtable` are available:

```gauss
library pubtable, quaids;
#include pubtable_quaids.src

struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);
call ptExport(coefTbl, "results.tex");

struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);   // 3x1: income, uncompensated, compensated
call ptExport(elasTbls[1], "income_elasticities.md");
```

See the [command reference](COMMAND_REFERENCE.md#reporting-optional-requires-pubtable)
for each adapter proc, and `examples/pubtable_export_example.e` for a full
runnable example.

## Limitations

- Curvature **imposition** ([quaidsCurvatureFit](command-reference/quaidsCurvatureFit.md))
  is only available for LA-AIDS/AIDS, at the sample mean, and its standard
  errors are a simplified delta-method approximation that is known to be
  unreliable when the estimated Cholesky factor has boundary (near-zero)
  entries -- see the [Methodology Notes](METHODOLOGY_NOTES.md#curvature-imposition-diewert-wales)
  for why this happens and why point estimates and the exact curvature
  property are unaffected. QUAIDS curvature imposition is deferred.
- No guaranteed convergence for the iterated estimator (or the curvature-
  constrained outer iteration built on top of it) -- see "Choosing A
  Model" above. `aCtl.relax` (Milestone 12) is an evidence-backed, opt-in
  mitigation, not a fix.
- IV is mandatory; there is no exogenous-total-expenditure estimation mode.
