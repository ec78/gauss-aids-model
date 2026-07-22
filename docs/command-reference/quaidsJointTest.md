# quaidsJointTest

## Purpose

Wald test of the joint null hypothesis that homogeneity **and** symmetry
hold, against an unconstrained fit.

## Format

```gauss
{ stat, pval, df } = quaidsJointTest(qOut);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md) called
  with `aCtl.homogenous = 0`. **Errors clearly if `qOut.homogenous /= 0`.**

## Returns

- `stat` (*scalar*) - Wald chi-squared statistic.
- `pval` (*scalar*) - p-value, `cdfchic(stat, df)`.
- `df` (*scalar*) - degrees of freedom, `(n - 1) + (n - 1)(n - 2)/2`
  (homogeneity's `n - 1` restrictions plus symmetry restrictions
  `gamma_ij = gamma_ji` for `i < j`, `i, j = 1..n-1`).

## Remarks

A symmetric gamma matrix with adding-up already imposed automatically
satisfies homogeneity too (row `i` sum = column `i` sum by symmetry = `0`
by adding-up), so there is no separate "symmetry only, given adding-up"
test on a fully unconstrained fit -- use this joint test, or the existing
symmetry-given-homogeneity test (built into [quaidsFit](quaidsFit.md)'s
own `qOut.symStat`/`qOut.symPval`, valid when `aCtl.homogenous == 1`) if
homogeneity itself isn't in question.

Validated for both size and power in
`tests/quaids_hypothesis_tests_test.e`.

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.homogenous = 0;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ stat, pval, df } = quaidsJointTest(qOut);
```

## Source

`quaidstests.src`

## See Also

[quaidsHomogeneityTest](quaidsHomogeneityTest.md), [quaidsFit](quaidsFit.md)
