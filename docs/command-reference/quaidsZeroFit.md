# quaidsZeroFit

## Purpose

Estimates AIDS/QUAIDS budget-share equations corrected for zero budget
shares (corner solutions) via the Shonkwiler & Yen (1999) two-step
procedure. Unconstrained only in this first pass (no homogeneity/symmetry
imposition on the corrected model). Silent, no printing -- see
[printQuaidsZero](printQuaidsZero.md).

## Format

```gauss
struct quaidsZeroOut zOut;
zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
```

## Parameters

- `w` (*TxN matrix*) - budget shares, may contain zeros.
- `intcpt` (*TxK matrix, or `0`*) - extra intercept-shifter variables
  (demographics etc.), or `0` for none. Also reused, unchanged, as part of
  the first-stage probit's regressor set -- see Remarks.
- `prices` (*TxN matrix*) - absolute log prices.
- `totexp` (*Tx1 vector*) - log total expenditure (treated as endogenous,
  same IV handling as [quaidsFit](quaidsFit.md)).
- `instr` (*TxH matrix*) - instruments for log total expenditure.
- `aCtl` (*`quaidsControl` structure*) - same fields as
  [quaidsFit](quaidsFit.md) (`linear`, `maxiter`, `err`, `alpha0`,
  `relax`, `b0`). **`aCtl.homogenous` must be `0`** -- errors clearly
  otherwise (see Remarks).

## Returns

`zOut` is a `quaidsZeroOut` structure:

- `probitB`, `probitSE` - first-stage probit coefficients/standard errors,
  one column per good, regressed on `intcpt~prices~totexp`.
- `probitConverged` - `n x 1`, `1` if that good's probit converged.
- `shareZeroFrac` - `n x 1`, observed fraction of zero shares per good
  (diagnostic).
- `b`, `v`, `se` - the corrected coefficient matrix (intercept | gamma |
  beta | `[lambda]` | `u` | a trailing `n x n` hazard/`delta` block, see
  Remarks), its covariance, and standard errors.
- `converged`, `iterations`, `finalErr` - outer-iteration diagnostics.

## Remarks

**The problem**: real survey/microdata routinely has corner solutions
(some households report zero expenditure on some goods). Fitting
[quaidsFit](quaidsFit.md) directly on such data is a known source of
bias -- the linear/log-linear share equation has no mechanism for a
censored dependent variable.

**The reformulation**: a literal, textbook Shonkwiler-Yen correction
rescales *every* regressor in equation `i` by that equation's own
first-stage probability `F_i`, which would break the shared-design-matrix
Kronecker-product identity every stage of `quaidsFit()` depends on.
Instead, `quaidsZeroFit()` divides the whole equation by `F_i`:

```
w_i = F_i*(alpha_i + sum_j gamma_ij*ln(p_j) + beta_i*lx [+ lambda_i*lx2]) + f_i*delta_i + e_i
```
becomes
```
w_i/F_i = alpha_i + sum_j gamma_ij*ln(p_j) + beta_i*lx [+ lambda_i*lx2] + (f_i/F_i)*delta_i + e_i/F_i
```
This turns the problematic *regressor* rescaling into a dependent-variable
transform (`wTilde_i = w_i/F_i`, computed once from the first-stage
probit) plus one new shared regressor column per equation
(`h_i = f_i/F_i`) -- structurally the same kind of addition as the `u`
(IV-residual) columns [quaidsFit](quaidsFit.md) already appends. The
shared-design-matrix machinery, and everything built on it, survives
unchanged.

**The diagonal-delta restriction**: appending `n` hazard columns to a
shared design matrix means the one-shot GLS solve initially estimates a
full `n x n` cross-equation `delta` block (every equation's response to
*every* good's hazard term), when Shonkwiler-Yen only wants the diagonal
(good `i`'s own hazard term in good `i`'s own equation). `quaidsZeroFit()`
imposes this via the same `design()`-based minimum-distance restriction
[quaidsFit](quaidsFit.md)'s own symmetry stage uses, just with a diagonal
(not symmetric) restriction pattern -- so `b`'s trailing block is `n`
rows, not the single row every other block has, with off-diagonal entries
forced to exactly zero.

**Scope, deliberately** (matching this project's established honest-
scoping precedent):

- **Unconstrained only.** Every good is directly, independently estimated
  (errors clearly if `aCtl.homogenous = 1`). Imposing homogeneity/symmetry
  on top of the Shonkwiler-Yen correction is real, additional work
  (combining two different minimum-distance restrictions in one pass)
  left for a follow-up milestone.
- **Standard errors are a simplified formula**
  (`V(vec(b)) = S .*. inv(gg)`, the classic SUR-with-shared-regressors
  covariance) -- it does *not* account for the nonlinear translog-price-
  index feedback the way [quaidsFit](quaidsFit.md)'s own Jacobian-
  corrected variance does, nor for first-stage probit/IV
  generated-regressor uncertainty. Matches this project's established
  precedent for [quaidsCurvatureFit](quaidsCurvatureFit.md)'s own
  "simplified, not a full sandwich" delta-method SE.
- **Adding-up does not hold exactly** for the corrected coefficients --
  a real, known property of Shonkwiler-Yen itself (each equation is
  rescaled by its own good-specific `F_i`), not a bug. No post-hoc
  renormalization is applied.
- **The probit regressors are not a separate argument**: they reuse
  `intcpt`/`prices`/`totexp`, matching this library's "no separate
  exogenous-mode arguments" philosophy.
- **Starting values**: scalar `aCtl.b0 = 0` uses the built-in Stone-index
  starting values. A supplied matrix must match `zOut.b`'s coefficient
  shape for a compatible zero-corrected fit (including the trailing
  hazard/`delta` block).

**A known, unresolved limitation**: GAUSS's built-in `glm()` can hard-crash
on some degenerate inputs (e.g. `error: Intel MKL ERROR ... DGELS`), a
failure mode that is *not* caught by this library's usual
`trap 1,1;`/`scalmiss()` guard idiom -- the same class of non-trappable
failure already documented for `eighv()` in
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md)'s own
Remarks. Not hardened against in this pass; if a fit crashes this way,
try a different starting sample, a coarser probit regressor set, or
report it.

**Uses GAUSS's built-in `glm()`** (`c:\gauss26\src\glm.src`, base GAUSS
runtime) for the first-stage probits -- no new package dependency.

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;          // 1 for AIDS, 0 for QUAIDS
aCtl.maxiter = 100;
aCtl.homogenous = 0;      // required

struct quaidsZeroOut zOut;
zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOut);
print "fraction of zero shares per good:" zOut.shareZeroFrac';
```

## Source

`quaidszerocorrect.src`

## See Also

[printQuaidsZero](printQuaidsZero.md), [quaidsFit](quaidsFit.md),
[quaidsCurvatureFit](quaidsCurvatureFit.md) (a similarly-scoped
"simplified SE, honestly documented limitation" precedent)
