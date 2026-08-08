# quaidsZeroFit

## Purpose

Estimates AIDS/QUAIDS budget-share equations corrected for zero budget
shares (corner solutions) via the Shonkwiler & Yen (1999) two-step
procedure. Supports both unconstrained estimation and homogeneity/symmetry
imposition on top of the correction (`aCtl.homogenous`). Silent, no
printing -- see [printQuaidsZero](printQuaidsZero.md).

## Format

```gauss
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
  `relax`, `b0`). `aCtl.homogenous = 0` (default) leaves every good
  directly, independently estimated; `aCtl.homogenous = 1` additionally
  imposes homogeneity and symmetry on the corrected model -- see Remarks.
  If supplied, `aCtl.b0` must match `zOut.bRaw`'s shape/basis, not
  `zOut.b`'s (see Remarks).

## Returns

`zOut` is a `quaidsZeroOut` structure:

- `n`, `n1` - number of goods, and `n - aCtl.homogenous` (the number of
  independently-estimated price/gamma rows; equation columns are always
  all `n`).
- `homogenous` - echoes `aCtl.homogenous`.
- `probitB`, `probitSE` - first-stage probit coefficients/standard errors,
  one column per good, regressed on `intcpt~prices~totexp`.
- `probitConverged` - `n x 1`, `1` if that good's probit converged.
- `shareZeroFrac` - `n x 1`, observed fraction of zero shares per good
  (diagnostic).
- `bRaw` - the raw, pre-recovery coefficient matrix, in the internal
  relative-price/reduced-row basis. This, not `b`, is the only valid
  shape/basis for a supplied `aCtl.b0` -- mirrors
  [quaidsFit](quaidsFit.md)'s own `qOut.homogB`/`aCtl.b0` pairing.
- `b`, `v`, `se` - the diagonal-delta-constrained coefficient matrix, in
  genuine **absolute-price** form (intercept | gamma | beta | `[lambda]`
  | `u` | a trailing `n x n` hazard/`delta` block, see Remarks), its
  covariance, and standard errors. Shape is `1+nint+nendog+nu+2n` rows,
  independent of `aCtl.homogenous` -- only the values differ.
- `symValid` - `== aCtl.homogenous`.
- `symStat`, `symPval`, `symDf` - the symmetry-given-homogeneity Wald
  test (`symDf = n1*(n1+1)/2` -- see Remarks for why this is *not*
  `n1*(n1-1)/2`).
- `bS`, `vS`, `seS` - the homogeneity+symmetry+diagonal-delta-constrained
  coefficient matrix (same shape/layout as `b`), covariance, and standard
  errors. `0` if `aCtl.homogenous = 0`.
- `bestB`, `bestV` - `= bS`/`vS` if `aCtl.homogenous`, else `= b`/`v` --
  what downstream post-estimation should be evaluated against.
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

**Homogeneity** (`aCtl.homogenous = 1`) is imposed the same way
[quaidsFit](quaidsFit.md) imposes it: only `n1 = n-1` (relative) price
columns are used as regressors, instead of `n`. Unlike `quaidsFit()`,
no equation *column* needs recovering via adding-up -- `quaidsZeroFit()`'s
`wTilde = w/F` residuals are not singular the way `quaidsFit()`'s raw
share residuals are (each equation is independently rescaled by its own
`F_i`), so every one of the `n` equations stays directly, independently
estimated throughout; only gamma's *rows* (which prices are used as
regressors) shrink from `n` to `n1`.

**Symmetry** (also `aCtl.homogenous = 1`) is imposed simultaneously with
the diagonal-delta restriction in one combined minimum-distance GLS
projection. This is *not* simply "symmetrize the `n1 x n1` gamma
sub-block, leave the redundant price row's column free" -- that was an
earlier draft's design, found by direct testing against a known-symmetric
DGP to leave the *recovered* absolute-price gamma measurably asymmetric.
The correct restriction additionally constrains the redundant row's
entries to be a deterministic function of that same symmetric sub-block
(via homogeneity's own row-sum-zero identity), which is why
`symDf = n1*(n1+1)/2` (the within-sub-block restrictions, `n1*(n1-1)/2`,
*plus* `n1` more from this cross-constraint) rather than the smaller
`n1*(n1-1)/2` an incomplete restriction would suggest -- confirmed by an
independent free-parameter count (`n1*n` free gammaRel parameters before
symmetry, `n1*(n1+1)/2` after, matching the standard AIDS
homogeneity+symmetry count `n(n-1)/2`). See `quaidszerocorrect.src`'s own
`_quaidsZeroSymDiagRestrict()` header for the full derivation.

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

**A real, previously-shipped bug, found and fixed**: prior to this, `b`'s
gamma columns were left in a mixed relative/absolute price basis (columns
`1..n-1` held `gamma_ij - gamma_in`, not genuine `gamma_ij`) --
`quaidsZeroFit()` never had an analog of [quaidsFit](quaidsFit.md)'s own
"RECOVERS ABSOLUTE PRICE EFFECTS FROM RELATIVE" conversion. Fixed in both
modes -- `b`/`bS` are now genuine absolute-price coefficients, matching
[quaidsFit](quaidsFit.md)'s own output convention. This changes `b`'s
*values* (not shape) relative to earlier releases.

**Scope, deliberately** (matching this project's established honest-
scoping precedent):

- **Standard errors are a simplified formula**
  (`V(vec(b)) = S .*. inv(gg)`, the classic SUR-with-shared-regressors
  covariance) -- it does *not* account for the nonlinear translog-price-
  index feedback the way [quaidsFit](quaidsFit.md)'s own Jacobian-
  corrected variance does, nor for first-stage probit/IV
  generated-regressor uncertainty. Matches this project's established
  precedent for [quaidsCurvatureFit](quaidsCurvatureFit.md)'s own
  "simplified, not a full sandwich" delta-method SE. The symmetry test
  statistic reuses this same simplified covariance for consistency.
- **Adding-up does not hold exactly** for the corrected coefficients --
  a real, known property of Shonkwiler-Yen itself (each equation is
  rescaled by its own good-specific `F_i`), not a bug. No post-hoc
  renormalization is applied.
- **The probit regressors are not a separate argument**: they reuse
  `intcpt`/`prices`/`totexp`, matching this library's "no separate
  exogenous-mode arguments" philosophy.
- **Starting values**: scalar `aCtl.b0 = 0` uses the built-in Stone-index
  starting values. A supplied matrix must match `zOut.bRaw`'s shape (the
  raw, pre-recovery internal basis) -- not `zOut.b`'s recovered shape,
  which can differ (both in row count when `aCtl.homogenous=1`, and in
  basis regardless of mode).

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
aCtl = quaidsControlCreate();
aCtl.linear = 0;          // 1 for AIDS, 0 for QUAIDS
aCtl.maxiter = 100;
aCtl.homogenous = 0;      // 1 to also impose symmetry

zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOut);
print "fraction of zero shares per good:" zOut.shareZeroFrac';
```

With homogeneity/symmetry imposed:

```gauss
aCtl.homogenous = 1;
zOutH = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
print "symmetry test p-value:" zOutH.symPval;
print zOutH.bestB;   // == zOutH.bS
```

## Source

`quaidszerocorrect.src`

## See Also

[printQuaidsZero](printQuaidsZero.md), [quaidsFit](quaidsFit.md),
[quaidsCurvatureFit](quaidsCurvatureFit.md) (a similarly-scoped
"simplified SE, honestly documented limitation" precedent)
