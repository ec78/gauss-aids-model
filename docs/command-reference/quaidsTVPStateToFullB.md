# quaidsTVPStateToFullB

## Purpose

Unpacks one period's TVP-AIDS reduced state vector (from
[quaidsTVPFit](quaidsTVPFit.md)'s `filteredState`/`smoothedState`, or any
column of either) into the full `n`-good coefficient matrix, in the same
`intcpt | prices | lx` row layout [quaidsElasFit](quaidsElasFit.md)/
`_quaidsElas()` expect. Closes two gaps the underlying reduced state
representation deliberately leaves open: equation `n`'s own
adding-up-implied coefficients (never separately estimated -- the reduced
state only carries `n1 = n-1` independently-estimated equations), and
converting the state's relative-price gamma sub-block to absolute-price
form.

## Format

```gauss
bFull = quaidsTVPStateToFullB(state, n1);
```

## Parameters

- `state` (*k_states x 1 vector*, `k_states = 2*n1 + n1*(n1+1)/2`) - one
  period's filtered or smoothed state (a column of
  [quaidsTVPFit](quaidsTVPFit.md)'s `filteredState`/`smoothedState`, or
  the lower-level `_quaidsTVPKalmanFit()`/`_quaidsTVPSmoothFit()` outputs
  directly).
- `n1` (*scalar*) - number of independently-estimated equations
  (`n = n1+1` goods total).

## Returns

`bFull` -- `(n+2) x n` matrix: intercept alpha (1 row) | absolute-price
gamma (`n` rows, symmetric, homogeneity-respecting) | beta (1 row), one
column per good (all `n`, including good `n`). Directly usable as the `b`
argument to `_quaidsElas()`/[quaidsElasFit](quaidsElasFit.md) with
`intcpt=1` and `aCtl.linear=1`.

## Remarks

**Equation `n`'s coefficients** are recovered via the same adding-up
identities [quaidsFit](quaidsFit.md)/[quaidsZeroFit](quaidsZeroFit.md)
already use: `alpha_n = 1 - sum(alpha_1..n1)`,
`beta_n = -sum(beta_1..n1)`.

**Relative-to-absolute-price gamma conversion** uses homogeneity's own
row-sum-zero identity (`sum_j gamma_abs[i,j] = 0` for every row `i`,
including row `n`): `gamma_abs[i,n] = -sum_{j=1}^{n1} gamma_rel[i,j]` for
`i=1..n1`, then row `n` itself via symmetry
(`gamma_abs[n,j] = gamma_abs[j,n]`) and homogeneity applied to row `n`.
Adding-up on gamma (every column also summing to zero) is **not**
separately imposed -- it falls out automatically once both symmetry and
homogeneity hold.

**No `sslib` dependency** -- `state` is accepted as a plain vector, so
this proc works whether `state` came from a filtered or smoothed path,
and requires no `library cmlmt, tsmt, sslib;` of its own (only
`quaidstvp.src` itself).

## Examples

```gauss
n1 = 3;

// At the final period's smoothed state:
bFull = quaidsTVPStateToFullB(tvOut.smoothedState[., tvOut.nobs], n1);

// Homogeneity/symmetry/adding-up all hold exactly on the result:
gammaAbs = bFull[2:n1+2, .];
"row sums (homogeneity):"; gammaAbs*ones(n1+1, 1);
"symmetry:"; maxc(maxc(abs(gammaAbs - gammaAbs')));
```

## Source

`quaidstvp.src`

## See Also

[quaidsTVPFit](quaidsTVPFit.md), [quaidsTVPElasFit](quaidsTVPElasFit.md)
(the next step -- elasticities from this same `bFull`),
[quaidsElasFit](quaidsElasFit.md)
