# quaidsControlCreate

## Purpose

Creates a `quaidsControl` struct populated with default estimation options.

## Format

```gauss
struct quaidsControl aCtl;
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
| `b0` | `0` | Optional user-supplied starting values; `0` = use linearized-AIDS starting values |
| `relax` | `1` | Under-relaxation factor for the iterated (`aCtl.maxiter>1`) fixed-point update, `(0,1]`; `1` = no damping (byte-identical to every release before Milestone 12). See Remarks |

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

`aCtl.relax` (Milestone 12) trades convergence speed for stability on the
iterated estimator's fixed-point update: `b_new = relax*b_solved +
(1-relax)*b_old` each iteration. A 200-seed sweep
(`tests/quaids_convergence_sweep.e`) found `relax=.75` measurably reduced
the estimator's convergence-failure rate versus the default `relax=1`;
more aggressive damping (`.5`, `.3`) did not help further and often made
things worse. Not a convergence guarantee -- see
[Feature Support Matrix](../FEATURE_SUPPORT_MATRIX.md#notes).

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

## Source

`quaidsutil.src`

## See Also

[getDefaultQuaidsControl](getDefaultQuaidsControl.md), [quaidsFit](quaidsFit.md)
