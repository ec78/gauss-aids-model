# quaidsCurvatureFit

## Purpose

Re-estimates a homogeneity+symmetry-constrained AIDS or QUAIDS fit under
an additional local curvature (Slutzky negative semidefiniteness)
restriction at the sample mean, via the Diewert-Wales (1987) Cholesky
reparametrization. Silent, no printing -- see
[printQuaidsCurvature](printQuaidsCurvature.md).

## Format

```gauss
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md), called
  with `aCtl.homogenous = 1` (either `aCtl.linear = 1` for AIDS or
  `aCtl.linear = 0` for QUAIDS -- both supported since Milestone 13).
  **Errors clearly if `aCtl.homogenous = 1` was not used.**
- `w` (*TxN matrix*) - budget shares, the same sample used to fit `qOut`.
- `prices` (*TxN matrix*) - absolute log prices, same sample.
- `totexp` (*Tx1 vector*) - log total expenditure, same sample.
- `aCtl` (*`quaidsControl` structure*) - `aCtl.err` controls the outer
  iteration's convergence tolerance, same as [quaidsFit](quaidsFit.md).
  `aCtl.relax` (see Remarks) is effectively required for QUAIDS.

## Returns

`cOut` is a `quaidsCurvOut` structure:

- `cholA` - the estimated `(n-1) x (n-1)` lower-triangular Cholesky factor.
- `b`, `v`, `se` - the curvature-constrained coefficient matrix (intercept
  | gamma | beta, plus a trailing `lambda` row for QUAIDS fits -- no `u`
  row either way, unlike `qOut.bestB`), its covariance, and standard
  errors.
- `gama` - the `n x n` curvature-constrained price-effect matrix (also
  `b`'s gamma block).
- `converged`, `iterations`, `finalErr` - outer-iteration diagnostics.
- `intcptPt`, `pricesPt`, `totexpPt` - the reference point curvature was
  imposed at (the sample means).
- `eigenvalues` - the Slutzky matrix's eigenvalues at the reference point;
  all should be `<= 0` to the outer iteration's convergence tolerance.

## Remarks

**Scope**: LA-AIDS/AIDS and QUAIDS, both supported. QUAIDS was initially
deferred at Milestone 10 -- its quadratic log-expenditure term adds a
`lambda`-dependent cross-term to the Slutzky matrix, entangling three
nonlinear parameter blocks instead of two -- but this resolved at
Milestone 13 using the same lag-then-solve trick
[quaidsFit](quaidsFit.md)'s own iteration already uses for its `beta`/
`lambda` coefficients: they're profiled out by OLS every outer round
rather than joining `optmt`'s search, which stays `vech(A)`-only,
unchanged in size. See `GOLD_STANDARD_TODO.md`'s Milestone 10 and 13
sections.

**QUAIDS needs `aCtl.relax`**: QUAIDS's curvature outer loop is
measurably less stable than AIDS's own two-block version (found
empirically, not assumed) -- undamped (`aCtl.relax=1`, the default) runs
can diverge to `NaN` within a handful of iterations. `aCtl.relax=.25`
(the same opt-in damping field [quaidsFit](quaidsFit.md) uses, Milestone
12) was found to converge cleanly on this library's own validation
fixture; not needed for AIDS, whose outer loop converges reliably
undamped.

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
standard errors for near-zero Cholesky-factor entries with caution. See
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) (Milestone
15) for a bootstrap alternative that does not share this boundary-
inference weakness.

**A real bug, found and fixed (Milestone 18)**: `se`'s individual cells
were silently scrambled relative to `b` from this proc's original
Milestone 10 release until it was caught while building
[quaidsCurvatureBootstrapCI](quaidsCurvatureBootstrapCI.md)'s ground-
truth cross-check (GAUSS's `reshape()` fills row-major, not column-major
like `vec()` -- a subtle, easily-missed distinction). Invisible to
shape/sign/finiteness checks, since those are permutation-invariant. `v`
(the full covariance) was never affected, only the reshaped `se` display
convenience. If you used `se` from a version before this fix, re-run
with the current release.

## Examples

```gauss
library optmt;

aCtl = quaidsControlCreate();
aCtl.linear = 0;         // 0 for QUAIDS, 1 for AIDS
aCtl.maxiter = 100;
aCtl.homogenous = 1;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

aCtl.relax = .25;    // recommended for QUAIDS; not needed for AIDS

cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);

call printQuaidsCurvature(cOut);
print "Slutzky eigenvalues at the sample mean:" cOut.eigenvalues';
```

## Source

`quaidscurvature.src`

## See Also

[printQuaidsCurvature](printQuaidsCurvature.md), [quaidsFit](quaidsFit.md),
[quaidsSlutzky](quaidsSlutzky.md),
[quaidsCurvatureBootstrapFit](quaidsCurvatureBootstrapFit.md) (a bootstrap
alternative to this proc's delta-method standard errors),
[quaidsCurvatureBootstrapCI](quaidsCurvatureBootstrapCI.md)
