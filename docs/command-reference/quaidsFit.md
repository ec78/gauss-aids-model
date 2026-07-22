# quaidsFit

## Purpose

Estimates a linearized AIDS, iterated AIDS, or QUAIDS demand system with
instrumental-variables treatment of log total expenditure, silently (no
console output), returning a `quaidsOut` structure.

## Format

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

## Parameters

- `w` (*TxN matrix*) - budget shares.
- `intcpt` (*TxK matrix, or scalar `0`*) - extra intercept-shifter
  variables (demographics etc.); `0` for none.
- `prices` (*TxN matrix*) - **absolute** log prices (not relative -- the
  proc converts internally).
- `totexp` (*Tx1 vector*) - log total expenditure, treated as endogenous.
- `instr` (*TxH matrix*) - instruments for log total expenditure.
- `aCtl` (*`quaidsControl` structure*) - see
  [quaidsControlCreate](quaidsControlCreate.md) for fields and defaults.

## Returns

`qOut` is a `quaidsOut` structure (~75 fields; see `src/quaids.sdf`),
grouped by phase:

- Metadata: `model` (`"LA-AIDS"`/`"AIDS"`/`"QUAIDS"`), `n`, `nint`, `ninst`,
  `nu`, `nobs`, name vectors (`xnam`, `wnam`, `znam`, `unam`), `intcptFull`.
- IV first-stage diagnostics: `ivB`, `ivSE`, `ivT`, `ivPvt`, `ivRsq`, `ivFstat`, ...
- Homogeneity-constrained (or unconstrained) estimates: `homogB`, `homogV`, ...
- Overidentification test (valid only when `ninst > nu`): `overidValid`,
  `overidGamma`, `overidFstat`, `overidPvf`, ...
- Symmetry test given homogeneity, and symmetry-constrained estimates
  (valid only when `aCtl.homogenous == 1`): `symValid`, `symStat`,
  `symPval`, `symcB`, `symcV`, ...
- Final absolute-price-form estimates: `b`, `v`, `bS`, `vS` -- the same
  four matrices the legacy [quaids](quaids.md) wrapper returns as
  `(b1, v1, b2, v2)`. If `aCtl.homogenous == 1`, `b`/`v` are
  homogeneity-constrained and `bS`/`vS` are homogeneity+symmetry
  constrained. If `aCtl.homogenous == 0`, `b`/`v` are the unconstrained
  (reparametrized) estimates and `bS`/`vS` are `0`.
- `bestB`/`bestV`: whichever is the most-constrained estimate actually fit
  (symmetric if homogeneity was imposed, else the recovered unconstrained
  fit) -- what elasticities ([quaidsElasFit](quaidsElasFit.md)) and the
  Slutzky diagnostic ([quaidsSlutzky](quaidsSlutzky.md)) should be
  evaluated against.

## Remarks

`quaidsFit()` prints nothing. Use [printQuaids](printQuaids.md) to
reproduce the estimation-stage console report from the returned struct, or
use the backward-compatible [quaids](quaids.md) wrapper to fit, print, and
reproduce the full legacy report in one call.

Log total expenditure is **always** treated as endogenous -- `instr` is a
required argument, not optional. There is no "no IV" estimation mode.

Set `aCtl.homogenous = 0` before calling
[quaidsHomogeneityTest](quaidsHomogeneityTest.md) or
[quaidsJointTest](quaidsJointTest.md), which both require an unconstrained
fit and error clearly otherwise.

Convergence is not guaranteed for the iterated estimator (`aCtl.maxiter >
1`): roughly half of random seeds in one synthetic-DGP family fail to
converge cleanly (see `GOLD_STANDARD_TODO.md`'s Milestone 3 findings).
Check `qOut.converged` and `qOut.iterations` after fitting.

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

call printQuaids(qOut);
print qOut.bestB;
```

Dataframe/column-name entry point instead of assembling matrices by hand:

```gauss
data = loadd("mydata.csv");
qOut = quaidsFull(data, shareVars, priceVars, "totexp", "instr", extraVars, aCtl);
```

## Source

`quaids.src`

## See Also

[printQuaids](printQuaids.md), [quaids](quaids.md),
[quaidsFull](quaidsFull.md), [quaidsControlCreate](quaidsControlCreate.md),
[quaidsHomogeneityTest](quaidsHomogeneityTest.md),
[quaidsJointTest](quaidsJointTest.md),
[quaidsElasFit](quaidsElasFit.md), [quaidsSlutzky](quaidsSlutzky.md)
