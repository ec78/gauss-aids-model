# quaidsCurvatureFit

## Purpose

Re-estimates a homogeneity+symmetry-constrained LA-AIDS/AIDS fit under an
additional local curvature (Slutzky negative semidefiniteness) restriction
at the sample mean, via the Diewert-Wales (1987) Cholesky
reparametrization. Silent, no printing -- see
[printQuaidsCurvature](printQuaidsCurvature.md).

## Format

```gauss
struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md), called
  with `aCtl.homogenous = 1` and `aCtl.linear = 1`. **Errors clearly if
  either requirement is not met** -- QUAIDS's quadratic term is not yet
  supported (see Remarks).
- `w` (*TxN matrix*) - budget shares, the same sample used to fit `qOut`.
- `prices` (*TxN matrix*) - absolute log prices, same sample.
- `totexp` (*Tx1 vector*) - log total expenditure, same sample.
- `aCtl` (*`quaidsControl` structure*) - `aCtl.err` controls the outer
  iteration's convergence tolerance, same as [quaidsFit](quaidsFit.md).

## Returns

`cOut` is a `quaidsCurvOut` structure:

- `cholA` - the estimated `(n-1) x (n-1)` lower-triangular Cholesky factor.
- `b`, `v`, `se` - the curvature-constrained coefficient matrix (same row
  layout as `qOut.bestB`: intercept | gamma | beta), its covariance, and
  standard errors.
- `gama` - the `n x n` curvature-constrained price-effect matrix (also
  `b`'s gamma block).
- `converged`, `iterations`, `finalErr` - outer-iteration diagnostics.
- `intcptPt`, `pricesPt`, `totexpPt` - the reference point curvature was
  imposed at (the sample means).
- `eigenvalues` - the Slutzky matrix's eigenvalues at the reference point;
  all should be `<= 0` to the outer iteration's convergence tolerance.

## Remarks

**Scope**: LA-AIDS/AIDS only (`aCtl.linear = 1`) in this release. QUAIDS's
quadratic log-expenditure term adds cross-terms to the Slutzky matrix that
entangle three nonlinear parameter blocks instead of two, a bounded but
not-yet-implemented follow-on -- see
`GOLD_STANDARD_TODO.md`'s Milestone 10 section.

**Why the sample mean, and not a caller-supplied point**: concavity of a
flexible functional form cannot be imposed globally (a standard
demand-theory result, not a limitation of this implementation) -- Diewert
& Wales impose it locally, at one reference point, and evaluating at the
observed sample mean is their own standard practice. A caller-supplied
arbitrary point (matching [quaidsElasFit](quaidsElasFit.md)'s flexibility)
would additionally require solving an implicit fixed-point equation for
the model-implied share at that point; deferred, not attempted in this
release.

**Estimation**: within each outer iteration (mirroring
[quaidsFit](quaidsFit.md)'s own translog-price-index iteration), the free
elements of the Cholesky factor (`vech(A)`, a small parameter vector) are
found by minimizing the IV residual sum of squares via GAUSS's `optmt`
package. For any candidate `A`, the remaining coefficients are exactly
identified by ordinary least squares once the reparametrized `gamma` is
substituted in as a fixed offset -- a profiled/concentrated nonlinear
least squares problem, so `optmt` only ever searches over `vech(A)`.

**Standard errors**: a simplified, homoskedastic nonlinear-least-squares
delta-method approximation (not a full SUR/GMM sandwich). **Known
limitation**: the estimated Cholesky factor frequently has boundary
(exactly zero) entries -- the constrained optimum lies on the edge of the
negative-semidefinite cone rather than its interior -- where classical
delta-method inference is known to be unreliable (the same boundary-
inference complication that arises for non-negativity-constrained
variance components elsewhere in econometrics). Point estimates and the
exact curvature property at the reference point are unaffected; treat
standard errors for near-zero Cholesky-factor entries with caution.

## Examples

```gauss
library optmt;

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);

call printQuaidsCurvature(cOut);
print "Slutzky eigenvalues at the sample mean:" cOut.eigenvalues';
```

## Source

`quaidscurvature.src`

## See Also

[printQuaidsCurvature](printQuaidsCurvature.md), [quaidsFit](quaidsFit.md),
[quaidsSlutzky](quaidsSlutzky.md)
