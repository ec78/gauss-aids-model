# quaidsControlCreate

## Purpose

Creates a `quaidsControl` struct populated with default estimation options.

## Format

```gauss
aCtl = quaidsControlCreate();
```

## Parameters

None.

## Returns

`aCtl` is a `quaidsControl` structure with fields:

| Field | Default | Meaning |
| --- | --- | --- |
| `linear` | `0` | `1` = LA-AIDS (linear); `0` = QUAIDS (quadratic log-expenditure term) |
| `maxiter` | `50` | `1` = one-step linearized AIDS with Stone price index; `>1` = iterate |
| `homogenous` | `1` | `1` = impose homogeneity (and test/report symmetry); `0` = unconstrained |
| `alpha0` | `0` | Fixed value of the translog price-index intercept `alpha_0` |
| `err` | `.0001` | Relative parameter-change convergence tolerance |
| `othnam` | `""` | Optional alternate variable names for printed output |
| `b0` | `0` | Optional user-supplied starting values; `0` = use built-in starting values. For `quaidsFit()`, a supplied matrix must match the reduced raw coefficient matrix shape used by the homogeneity stage (`qOut.homogB`). For `quaidsZeroFit()`, it must match the raw, pre-recovery coefficient shape (`zOut.bRaw`) -- **not** `zOut.b`, which is in recovered, absolute-price form and can differ in both shape and basis |
| `relax` | `1` | Under-relaxation factor for the iterated (`aCtl.maxiter>1`) fixed-point update, `(0,1]`; `1` = no damping. See Remarks |

Structure-inference return typing means callers do not need to pre-declare
`struct quaidsControl aCtl;` before assignment.

## Remarks

Set `aCtl.linear = 1` for LA-AIDS/iterated AIDS, `0` for QUAIDS. Set
`aCtl.maxiter = 1` for the one-step Stone-index LA-AIDS special case
(implies `aCtl.linear`'s value is irrelevant to the price index used, since
`maxiter == 1` always uses the Stone index regardless). Set
`aCtl.homogenous = 0` to fit unconstrained, e.g. before calling
`quaidsHomogeneityTest`/`quaidsJointTest`, which both require an
unconstrained fit.

**The defaults above (`linear=0`, `maxiter=50`) select QUAIDS with
iteration -- the highest-risk combination on the support-tier list**: a
200-seed sweep measured a 76% combined convergence-failure rate for QUAIDS
at these settings, versus 0% for `maxiter=1` LA-AIDS. These coded defaults
are kept unchanged for `0.x` compatibility; always check `qOut.converged`
after fitting, and see the [README's Model & Feature Support
Tiers](../../README.md#model--feature-support-tiers) or the [Feature
Support Matrix's Support Tier
Summary](../FEATURE_SUPPORT_MATRIX.md#support-tier-summary) before relying
on an unexamined default `quaidsControlCreate()` call.

`aCtl.homogenous` is a historical field-name misspelling retained for
source compatibility through the `0.x` series. New application code
should read/write it via [quaidsGetHomogeneity](quaidsGetHomogeneity.md)/
[quaidsSetHomogeneity](quaidsSetHomogeneity.md) instead of referencing
the misspelled field name directly -- see the
[compatibility policy](../../README.md#compatibility-policy) and
[`docs/public-api.json`](../public-api.json).

`aCtl.relax` trades convergence speed for stability on the
iterated estimator's fixed-point update: `b_new = relax*b_solved +
(1-relax)*b_old` each iteration. A 200-seed sweep
(`tests/quaids_convergence_sweep.e`) found `relax=.75` measurably reduced
the estimator's convergence-failure rate versus the default `relax=1`;
more aggressive damping (`.5`, `.3`) did not help further and often made
things worse. Not a convergence guarantee -- see
[Feature Support Matrix](../FEATURE_SUPPORT_MATRIX.md#notes).

## Examples

```gauss
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

## Source

`quaidsutil.src`

## See Also

[getDefaultQuaidsControl](getDefaultQuaidsControl.md), [quaidsFit](quaidsFit.md),
[quaidsSetHomogeneity](quaidsSetHomogeneity.md), [quaidsGetHomogeneity](quaidsGetHomogeneity.md)
