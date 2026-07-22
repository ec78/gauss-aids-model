# quaidsHomogeneityTest

## Purpose

Wald test of the null hypothesis that homogeneity holds
(`sum_j gamma_ij = 0` for every good `i`), against an unconstrained fit.

## Format

```gauss
{ stat, pval, df } = quaidsHomogeneityTest(qOut);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md) called
  with `aCtl.homogenous = 0`. **Errors clearly if `qOut.homogenous /= 0`.**

## Returns

- `stat` (*scalar*) - Wald chi-squared statistic.
- `pval` (*scalar*) - p-value, `cdfchic(stat, df)`.
- `df` (*scalar*) - degrees of freedom, `n - 1` (one row-sum restriction
  per independently-estimated equation; equation `n` is recovered via
  adding-up and contributes no new information to the test).

## Remarks

Reads `qOut.b` (the final, absolute-price-form, unconstrained gamma
matrix) and `qOut.v` (its covariance), and builds the restriction vector
via a selection matrix `L`: `stat = (L'*vec(b))' * inv(L'*V*L) *
(L'*vec(b))`. This is the standard delta-method/Wald construction for a
linear hypothesis. Asymptotically valid under the same regularity
conditions as the rest of this library's asymptotic inference -- large-T,
correctly specified instruments, no claim of finite-sample exactness.

Validated for both size (does it correctly fail to reject a true null?)
and power (does it correctly reject a false null?) in
`tests/quaids_hypothesis_tests_test.e` -- not just "it runs."

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.homogenous = 0;   // required: unconstrained fit

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ stat, pval, df } = quaidsHomogeneityTest(qOut);
if pval < 0.05;
    print "Reject homogeneity at the 5% level.";
endif;
```

## Source

`quaidstests.src`

## See Also

[quaidsJointTest](quaidsJointTest.md), [quaidsFit](quaidsFit.md)
