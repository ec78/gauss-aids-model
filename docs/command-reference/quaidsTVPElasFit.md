# quaidsTVPElasFit

## Purpose

Income and price elasticities at one period's TVP-AIDS state, via
[quaidsTVPStateToFullB](quaidsTVPStateToFullB.md) (state -> full `n`-good
coefficient matrix) then `_quaidsElas()` (the same already-correct,
already-tested elasticity math [quaidsElasFit](quaidsElasFit.md) uses) --
no new elasticity formula is introduced.

## Format

```gauss
{ er, ep, epc } = quaidsTVPElasFit(state, n1, prices, totexp, aCtl);
```

## Parameters

- `state` (*k_states x 1 vector*) - one period's filtered or smoothed
  state (`k_states = 2*n1 + n1*(n1+1)/2` -- see
  [quaidsTVPStateToFullB](quaidsTVPStateToFullB.md)).
- `n1` (*scalar*) - number of independently-estimated equations
  (`n = n1+1` goods total).
- `prices` (*n x 1 vector*) - ABSOLUTE log prices (not the relative
  prices the state space itself was built from) at the evaluation point.
- `totexp` (*scalar*) - log total expenditure at the evaluation point.
- `aCtl` (*`quaidsControl` structure*) - only `aCtl.alpha0` and
  `aCtl.linear` are read. **`aCtl.linear` must be `1`** -- this
  initiative's own Stage 1 scope (Stone index, no quadratic term) means
  the recovered state has no lambda row to read otherwise; errors clearly
  if not. A caller with no other reason to build a `quaidsControl` can
  just call `quaidsControlCreate()` then set `aCtl.linear = 1`.

## Returns

- `er` (*n x 1 vector*) - income elasticities.
- `ep` (*n x n matrix*) - uncompensated (Marshallian) price elasticities.
- `epc` (*n x n matrix*) - compensated (Hicksian) price elasticities.

## Remarks

**Point elasticities only** -- no standard errors, unlike
[quaidsElasFit](quaidsElasFit.md)'s own delta-method SEs. Propagating the
reduced state's own filtered/smoothed covariance through
[quaidsTVPStateToFullB](quaidsTVPStateToFullB.md)'s recovery step is real,
tractable future work (that recovery is itself linear in the state), not
attempted yet.

**No `sslib` dependency** -- `state` is a plain vector, so "filtered vs.
smoothed" is entirely the caller's own choice of which column of
[quaidsTVPFit](quaidsTVPFit.md)'s `filteredState`/`smoothedState` to pass.

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;   // required

n1 = tvOut.n1;
finalPrices = prices[tvOut.nobs, .]';
finalTotexp = totexp[tvOut.nobs];

{ er, ep, epc } = quaidsTVPElasFit(tvOut.smoothedState[., tvOut.nobs], n1, finalPrices, finalTotexp, aCtl);

"income elasticities:"; er;
```

## Source

`quaidstvpelas.src`

## See Also

[quaidsTVPFit](quaidsTVPFit.md), [quaidsTVPStateToFullB](quaidsTVPStateToFullB.md),
[quaidsElasFit](quaidsElasFit.md), [quaidsElas](quaidsElas.md)
