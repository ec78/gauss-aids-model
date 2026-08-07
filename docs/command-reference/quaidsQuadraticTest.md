# quaidsQuadraticTest

## Purpose

Wald test of the null hypothesis that the QUAIDS quadratic log-expenditure
term is not needed (`lambda_i = 0` for every good `i`) -- i.e. whether
plain AIDS would fit equally well -- against an unconstrained QUAIDS fit.

## Format

```gauss
{ stat, pval, df } = quaidsQuadraticTest(qOut);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md) called
  with `aCtl.homogenous = 0` **and** `aCtl.linear = 0`. **Errors clearly
  if `qOut.homogenous /= 0` (needs an unconstrained fit, same as
  [quaidsHomogeneityTest](quaidsHomogeneityTest.md)/
  [quaidsJointTest](quaidsJointTest.md)) or if `qOut.linear /= 0`** -- an
  AIDS fit never estimates `lambda` at all, so there is nothing to
  restrict.

## Returns

- `stat` (*scalar*) - Wald chi-squared statistic.
- `pval` (*scalar*) - p-value, `cdfchic(stat, df)`.
- `df` (*scalar*) - degrees of freedom, `n - 1` (`lambda_i = 0` for each
  independently-estimated equation `i`; equation `n`'s restriction is
  implied by adding-up once the other `n - 1` hold, same reasoning as
  [quaidsHomogeneityTest](quaidsHomogeneityTest.md)'s own `df`).

## Remarks

**Only ever callable on a QUAIDS fit.** When `aCtl.linear = 1` (AIDS),
`lambda` is not a nuisance parameter fixed at zero -- it is not estimated
at all, and `qOut.b` is one row shorter with no `lambda` row to select.
There is no way to "test whether QUAIDS is needed starting from an AIDS
fit"; this test only makes sense as a check on an already-fitted QUAIDS
model, telling you whether the simpler AIDS specification would have done
just as well.

A rejection (small `pval`) means the quadratic term is doing real work --
QUAIDS's added flexibility is statistically justified for this data. A
failure to reject means the data are consistent with plain AIDS, and the
simpler, more stable (see [quaidsFit](quaidsFit.md)'s convergence
Remarks) model may be preferable.

Validated for both size and power in
`tests/quaids_hypothesis_tests_test.e`, reusing the existing
`quadratic=1`/`quadratic=0` fixture pair from
`tests/quaidsfixtures.src`'s `_quaidsSyntheticDGP()`.

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.homogenous = 0;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ stat, pval, df } = quaidsQuadraticTest(qOut);
if pval > 0.05;
    print "Quadratic term not statistically justified -- consider AIDS (aCtl.linear=1).";
endif;
```

## Source

`quaidstests.src`

## See Also

[quaidsHomogeneityTest](quaidsHomogeneityTest.md),
[quaidsJointTest](quaidsJointTest.md), [quaidsFit](quaidsFit.md)
