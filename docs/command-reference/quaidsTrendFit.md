# quaidsTrendFit

## Purpose

A cheap, one-shot screening diagnostic for whether an AIDS demand
system's coefficients show a linear trend over time. Widens the shared
GLS design matrix with trend-interacted regressors and returns a joint
Wald test of whether the entire trend block is zero, along with the
individual trend-slope coefficients and standard errors. Silent, no
printing -- see [printQuaidsTrend](printQuaidsTrend.md).

**Not genuine time-varying-parameter estimation.** A planned
Kalman-filter-based random-walk-coefficient extension
(`quaidsTVPFit()`) covers that; this proc exists as the fast, near-free
first check of whether that much larger estimation effort is even
warranted on a given dataset.

## Format

```gauss
tOut = quaidsTrendFit(w, intcpt, prices, totexp, instr, aCtl);
```

## Parameters

- `w` (*TxN matrix*) - budget shares.
- `intcpt` (*TxK matrix, or `0`*) - extra intercept-shifter variables, or
  `0` for none. Not itself trended -- see Remarks.
- `prices` (*TxN matrix*) - absolute log prices.
- `totexp` (*Tx1 vector*) - log total expenditure (treated as endogenous,
  same IV handling as [quaidsFit](quaidsFit.md)).
- `instr` (*TxH matrix*) - instruments for log total expenditure.
- `aCtl` (*`quaidsControl` structure*) - only `aCtl.homogenous` and
  `aCtl.othnam` are read. **`aCtl.homogenous` must be `1`** -- errors
  clearly otherwise (see Remarks). `aCtl.linear`/`aCtl.maxiter` are not
  read; this proc always uses the LA-AIDS/Stone-index specification
  internally.

## Returns

`tOut` is a `quaidsTrendOut` structure:

- `n`, `nobs`, `nint`, `n1` - dimensions (`n1 = n - 1`, always, since
  homogeneity is required).
- `ng`, `ngTrend` - row counts of the level block (`1+nint+n1+1(lx)+nu`)
  and the trend block (`1(t)+n1+1(t*lx)`).
- `xnam`, `wnam`, `unam`, `intcptFull`, `u` - names and the first-stage IV
  residuals, same conventions as [quaidsFit](quaidsFit.md)'s `qOut`
  fields.
- `b0`, `se0`, `t0`, `pvt0` - the level-block coefficients (row order
  `intcpt | prices[.,1:n1] | lx | u`, the same convention
  [quaidsElasFit](quaidsElasFit.md)'s `b` docstring uses, minus the QUAIDS
  `lx2` row since this proc is always LA-AIDS), standard errors, t-stats,
  p-values.
- `b1`, `se1`, `t1`, `pvt1` - the trend-slope-block coefficients (row
  order `t | t.*prices[.,1:n1] | t.*lx`, a mean-centered linear trend
  interacted with each regressor), standard errors, t-stats, p-values.
- `trendStat`, `trendPval`, `trendDf` - the joint Wald test that the
  entire trend block is zero, across the `n-1` independently-estimated
  equations (`trendDf = ngTrend*(n-1)`).

## Remarks

**Design decisions, made explicit:**

- **Always LA-AIDS internally.** A screening diagnostic does not need the
  nonlinear translog price index's own iteration -- that iteration exists
  in [quaidsFit](quaidsFit.md) to get the best possible point estimate;
  here the question is only "is there a hint of drift worth investigating
  further," and the Stone index's well-documented approximation bias (see
  [quaidsFit](quaidsFit.md)'s own Milestone-3-era history) is a much
  smaller concern for a screen than for a reported final estimate.
- **`aCtl.homogenous == 1` is required.** The exact adding-up/homogeneity
  guarantee this proc relies on (below) only holds cleanly when the
  reference good's own price is dropped from the shared design, exactly
  as [quaidsFit](quaidsFit.md)'s own homogeneous branch does. An
  unconstrained version is a real, separate extension, not attempted
  here.
- **Only a bare linear trend is interacted with prices/log expenditure**;
  existing `intcpt` shifters are not separately trended. Whether a
  demographic shifter's own effect drifts over time is a bigger, separate
  question than this quick screen answers.
- **The IV control-function term `u` is not trended** -- total-expenditure
  endogeneity is treated once (one first-stage regression), not as itself
  evolving over time, matching the same design decision documented for
  the planned full TVP model.
- **Simplified covariance**: omits the generated-regressor correction for
  `u`'s own first-stage sampling variability (the term
  [quaidsFit](quaidsFit.md)'s full covariance includes) -- the same class
  of documented simplification this library already ships in
  [quaidsRobustFit](quaidsRobustFit.md)'s "simplified bread" and
  [quaidsZeroFit](quaidsZeroFit.md)'s own covariance. With that term
  omitted, the sandwich collapses exactly (verified by Kronecker algebra,
  not just asserted) to the textbook SUR-with-shared-regressors form
  `v = S[1:n-1,1:n-1].*.inv(gg)`.

**The central mathematical claim, and why it needs no separate
restriction-imposition step**: adding-up and homogeneity hold *exactly*
(to floating-point precision) on *both* the level (`b0`) and trend-slope
(`b1`) coefficient blocks, purely as a consequence of the shared-design-
matrix GLS mechanism every AIDS estimate in this library already uses.
Budget shares sum to 1 identically, and every equation is fit against the
*identical* regressor set, so `b`'s own column-sum for any regressor
other than the bare constant is forced to exactly 0 (and the constant's
own column-sum is forced to exactly 1) -- a pure linear-algebra
consequence of OLS/GLS on data satisfying a linear identity. This extends
unchanged to the added trend-interaction columns; nothing new has to be
derived or separately enforced for them.

**Verified directly, not just asserted**: confirmed this claim holds to
~1e-16 precision on synthetic data where adding-up is exact by
construction (`tests/quaidstrend_test.e`), and separately confirmed that
running against *real* published data (`Blanciforti86`, rounded to a few
decimals in its original 1986 source) gives only *approximate* (~1e-3)
adding-up -- traced directly to the real data's own rounding, not a bug
in this proc's math.

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl.homogenous = 1;   // required

tOut = quaidsTrendFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsTrend(tOut);

if tOut.trendPval < 0.05;
    print "Some evidence of coefficient drift -- a full TVP-AIDS fit may be worthwhile.";
endif;
```

## Source

`quaidstrend.src`

## See Also

[printQuaidsTrend](printQuaidsTrend.md), [quaidsFit](quaidsFit.md),
[quaidsQuadraticTest](quaidsQuadraticTest.md) (a similarly-scoped "run
this cheap test before committing to a more complex model" precedent)
