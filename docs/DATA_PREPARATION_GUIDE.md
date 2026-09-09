# Data Preparation Guide

Everything this library needs from your raw data, in one place, before you
call [quaidsFit](command-reference/quaidsFit.md)/[quaidsFull](command-reference/quaidsFull.md).
[examples/00_real_data_quickstart.e](../examples/00_real_data_quickstart.e)
is the worked, runnable counterpart to this guide -- it links back to the
relevant section here at each preparation decision it makes. Read this
guide once, then use the [checklist](#final-input-contract-checklist) at
the end as a pre-flight reference on your own data.

## 1. Budget Shares and Total Expenditure

`w` (`TxN`) is the fraction of total expenditure spent on each of `N`
goods -- `w_i = expenditure_i / totalExpenditure`. Every row must sum to
`1`: this is an accounting identity, not something the library estimates,
and [quaidsPreflight](command-reference/quaidsPreflight.md) treats a
violation as a hard failure (`shareAddOk`, tolerance `1e-6`).

**Real, rounded data will not sum to floating-point-exact `1`.** A
published or survey dataset reported to 2-3 decimal places routinely sums
to `1 +/- 0.001`, which is well outside preflight's tight tolerance --
this is not a data error, it is normal rounding, and it will hard-fail
preflight if not corrected first. The standard fix is to row-normalize
before fitting:

```gauss
w = w./sumc(w');
```

This forces exact adding-up while barely perturbing values that were
already within a fraction of a percent of summing to 1. See
[examples/00_real_data_quickstart.e](../examples/00_real_data_quickstart.e)
for this exact fix applied to real published data.

`totexp` (`Tx1`) must be **log** total expenditure, and it must be the
log of the *same* total that `w`'s denominator uses -- if `w_i =
x_i/X`, then `totexp = ln(X)`, not the log of some other aggregate
(household income, a different expenditure total, etc.). Using a
mismatched total silently breaks the model's adding-up identity even
though `w` itself still sums to 1.

## 2. Price and Expenditure Transformations and Unit Consistency

`prices` (`TxN`) must be **absolute log prices**, not relative and not
levels -- `quaidsFit()` converts to relative prices internally, so pass
raw log price levels, e.g. `prices = ln(priceLevels)` if your source data
has price indices or levels rather than logs.

`totexp` must likewise be **logged**, not a raw expenditure level (see
above).

**A common trap when using [quaidsFull](command-reference/quaidsFull.md)
with a dataframe**: `quaidsFull()` selects columns by name only -- it does
**not** transform them. If your loaded CSV has raw price/expenditure
levels, you must add already-logged columns to the dataframe yourself
(`data = dfaddcol(data, "lnP1", ln(data[.,"P1"]));`) before calling
`quaidsFull()`, and pass the logged column names as `priceVars`/
`totexpVar`. The matrix API (`quaidsFit()`) has the identical requirement,
just applied to a plain matrix instead of a named column.

**Unit consistency**: every price column should use the same base period/
units as every other price column, and consistently across all `T`
observations -- a price index rebased partway through a panel, or one
good's price in different units than the rest, will not error out
mechanically but will produce a meaningless fit. There is no automated
check for this; it is a data-construction responsibility.

## 3. Good/Category Ordering Across Shares and Prices

`shareVars[i]` and `priceVars[i]` (or matrix columns `w[.,i]`/
`prices[.,i]`) must refer to the **same good**. Matching is strictly by
position, never by name -- there is no name-matching magic in
`quaidsFull()`, and the matrix API has no names to match at all. Getting
this wrong does not error; it silently fits a demand system where good 1's
share is regressed against good 3's price, for example, with plausible-
looking but meaningless output. Double-check column order explicitly when
building `shareVars`/`priceVars`, especially when they come from two
separate `loadd()` calls or two independently-maintained lists.

The last good (`N`) is used internally as the reference good when forming
relative prices, and its own row/column carries a "reference" label in
some printed reports (e.g. [printQuaids](command-reference/printQuaids.md)'s
"Reference price (absolute)" block) -- this has no effect on which data
you supply or how to interpret coefficients for other goods, it's purely
an internal computational convention.

## 4. Missing Values, Invalid Observations, Zeros, and Corner Solutions

**Missing/non-finite values**: `quaidsFit()` has no built-in imputation.
[quaidsPreflight](command-reference/quaidsPreflight.md)'s `finiteOk`/
`nonFiniteCount` check is a hard failure on any `NaN`/missing/infinite
cell in `w`/`intcpt`/`prices`/`totexp`/`instr` -- clean or drop affected
rows before fitting.

**Negative shares**: a hard preflight failure
(`negativeShareCount`). Legitimate budget shares are never negative;
a negative value almost always indicates a data-construction error
(e.g. a net-purchase or refund convention that produced negative
expenditure) rather than something this library should absorb silently.

**Zero shares (corner solutions)**: real survey/microdata routinely has
households reporting zero expenditure on some goods -- a fundamentally
different situation from a data error. `quaidsFit()`'s linear/log-linear
share equation has no mechanism for a censored dependent variable, and
fitting it directly on data with many zeros is a known source of bias.
If zero shares are common in your data (`quaidsPreflight`'s
`zeroShareCount` warning is your signal), use
[quaidsZeroFit](command-reference/quaidsZeroFit.md) (the Shonkwiler-Yen
correction) instead of `quaidsFit()` -- see
[Zero Budget Shares](USAGE_GUIDE.md#zero-budget-shares-corner-solutions)
in the usage guide. A handful of incidental zeros in an otherwise
zero-share-free dataset is a warning, not necessarily a reason to switch
estimators; use judgment based on how common they are.

## 5. Instrument Selection and Weak-Instrument Diagnostics

Log total expenditure is **always** treated as endogenous -- `instr` is a
required argument with no "exogenous" mode. You need at least one
instrument column genuinely correlated with total expenditure but
plausibly uncorrelated with the demand-equation errors; classic choices
in the applied-demand literature are household income, total expenditure
across a broader set of categories than the system being estimated (the
approach [examples/00_real_data_quickstart.e](../examples/00_real_data_quickstart.e)
uses -- total food expenditure instrumented by total expenditure across
all commodity groups), or lagged expenditure in panel data.

Supplying more instrument columns than endogenous regressors
(`ninst > nu`) activates the overidentification test
(`qOut.overidValid`/`overidFstat`/`overidPvf`) automatically -- a useful
free diagnostic if you have more than one candidate instrument.

**Always check instrument strength before trusting a fit.**
[quaidsPreflight](command-reference/quaidsPreflight.md)'s `ivFstat`/
`weakIV` fields report the first-stage F statistic (the same one
`quaidsFit()`'s own first-stage regression produces) and flag `weakIV`
using the conventional `F < 10` screening rule. This is a screening
diagnostic, not a substitute for a full weak-IV analysis, but a low
first-stage F is a real warning sign that your instrument may not be
informative enough for reliable IV estimates.

## 6. Demographic Intercept Shifters

`intcpt` (`TxK`, or scalar `0` for none) holds extra variables that shift
each good's intercept -- household size, region indicators, a time trend,
or similar demographic/control variables. These are distinct from, and
should not be confused with, `clusterId` (a grouping label used only for
cluster-robust standard errors) or `weight` (a sampling weight used only
for weighting the point estimate/SE) -- three separate arguments with
three separate roles, all of which can be non-trivial (non-zero/non-`0`)
in the same call.

Each shifter gets its own coefficient row per good in the output (visible
via `qOut.xnam`/`qOut.nint`). By the model's own adding-up identity, a
shifter's coefficients sum to `0` across the `N` goods -- it reallocates
budget shares between goods for a given household characteristic, it does
not change total expenditure itself.

## 7. Sampling Weights, Clusters, Replicate Weights, and Strata

This library supports three, independent, opt-in mechanisms for survey/
microdata designs -- pick whichever matches what your dataset ships:

- **`weight`** (`quaidsFit`'s own optional keyword argument): a genuine
  sampling weight applied to the point estimate itself (the standard
  survey-WLS `sqrt(weight)` scaling), with a matching weighted sandwich SE
  via `quaidsRobustFit`'s own `weight` argument (a **different** scaling
  convention -- see that page). Use when your data has known sampling
  probabilities/design weights.
- **`clusterId`** (`quaidsRobustFit`'s keyword argument): a group-label
  vector for cluster-robust standard errors (CR1 small-sample
  correction). Requires **at least two clusters** -- `quaidsRobustFit`
  errors clearly ("cluster-robust SE require at least two clusters") on a
  single-cluster input, and cluster-robust asymptotics are only reliable
  with a reasonably large number of clusters (a few dozen or more is the
  usual applied guidance; this library does not enforce a minimum beyond
  two).
- **`replicateWeights`/`scaleFactor`** ([quaidsReplicateWeightFit](command-reference/quaidsReplicateWeightFit.md)):
  for survey extracts that ship **pre-computed replicate weight columns**
  (a common jackknife/BRR/Fay's-BRR design) rather than requiring you to
  implement your own resampling. Both arguments are always required --
  this library does not auto-detect or assume a specific replication
  design.

**Explicitly not supported**: formal strata as a concept distinct from
clustering, and finite-population correction. If your survey design
depends on either, this library's weighted/clustered/replicate-weight
support covers the point estimate and a reasonable variance estimate, but
is not a complete design-based survey estimator -- see the [Feature
Support Matrix](FEATURE_SUPPORT_MATRIX.md#support-tier-summary) for the
current, honest scope.

## 8. Minimum Sample/Design Size and Recommended Preflight Checks

There is no single hardcoded minimum sample size -- the right check is
whether your design has enough degrees of freedom and variation for a
stable fit, which is exactly what
[quaidsPreflight](command-reference/quaidsPreflight.md) screens for
automatically:

- `designCols`/`designDf`/`designInvOk` -- whether the shared regressor
  design matrix is even invertible given your sample size, number of
  goods, intercept shifters, and instruments. A hard failure here means
  your data cannot support this model at all (too few observations
  relative to parameters, or a rank-deficient design, e.g. a constant
  column duplicated by a shifter).
- `minPriceStd`/`totexpStd`/`minInstrStd` and their low-variation
  warnings -- a price, expenditure, or instrument column with too little
  variation across observations gives the estimator little to identify
  parameters from, even if the design is technically invertible.
- `nClusters`/`minClusterSize`/`singletonClusters` -- if using
  `clusterId`, these report whether you have enough clusters and whether
  any cluster is too small (a singleton cluster contributes no
  within-cluster variance information).
- `convergenceRisk` (`0`/`1`/`2`) -- a simple screen combining the above;
  higher risk does not guarantee a convergence failure, but is worth
  noting before spending time on an iterated (`aCtl.maxiter > 1`) fit --
  see [Model & Feature Support Tiers](../README.md#model--feature-support-tiers)
  for the measured convergence-failure rates for iterated AIDS/QUAIDS.

Always run `quaidsPreflight()` and resolve every reported error (not just
warnings) before fitting -- it is silent, estimator-free, and cheap to
run, and every issue above is far easier to diagnose from its structured
output than from a confusing downstream numerical failure.

## Final Input-Contract Checklist

Work through this list on your own data before your first
`quaidsFit()`/`quaidsFull()` call:

- [ ] `w` is `TxN`, every row sums to exactly `1` (row-normalize real/
      rounded data first -- see [Section 1](#1-budget-shares-and-total-expenditure)).
- [ ] `prices` is `TxN` **log** prices (not levels), consistent units/base
      across every column and row (see [Section 2](#2-price-and-expenditure-transformations-and-unit-consistency)).
- [ ] `totexp` is `Tx1` **log** total expenditure, computed as the log of
      the exact same total `w`'s shares were divided by.
- [ ] If using `quaidsFull()`, the dataframe's price/total-expenditure
      columns are already logged -- `quaidsFull()` does not transform
      them for you.
- [ ] `shareVars`/`priceVars` (or matrix columns) are ordered so good `i`
      matches good `i` in both -- double-checked explicitly, not assumed
      (see [Section 3](#3-goodcategory-ordering-across-shares-and-prices)).
- [ ] No missing/non-finite values remain in any of `w`/`intcpt`/
      `prices`/`totexp`/`instr`.
- [ ] No negative shares remain (a real data error if present).
- [ ] A decision has been made on zero shares: `quaidsFit()` if rare/
      absent, [quaidsZeroFit](command-reference/quaidsZeroFit.md) if
      common (see [Section 4](#4-missing-values-invalid-observations-zeros-and-corner-solutions)).
- [ ] `instr` has at least one column genuinely correlated with `totexp`,
      and its first-stage strength has been checked (not assumed) via
      `quaidsPreflight`'s `ivFstat`/`weakIV`.
- [ ] `intcpt` is set deliberately (`0` for none, or a `TxK` matrix of
      demographic shifters) -- not confused with `clusterId`/`weight`.
- [ ] A survey-design mechanism has been chosen deliberately if
      applicable: `weight`, `clusterId`, or `replicateWeights`/
      `scaleFactor` (see [Section 7](#7-sampling-weights-clusters-replicate-weights-and-strata)) --
      or explicitly none, for a simple random sample.
- [ ] `quaidsPreflight()` has been run and every reported error resolved
      (warnings reviewed, not necessarily all resolved).

## See Also

- [examples/00_real_data_quickstart.e](../examples/00_real_data_quickstart.e) --
  every checklist item above applied to a real published dataset.
- [Usage Guide](USAGE_GUIDE.md) -- choosing an API and a model.
- [Feature Support Matrix](FEATURE_SUPPORT_MATRIX.md) -- what is and is
  not supported, by feature.
- [quaidsPreflight command reference](command-reference/quaidsPreflight.md) --
  the full diagnostic field list.
