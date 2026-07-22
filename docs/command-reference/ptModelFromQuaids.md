# ptModelFromQuaids

## Purpose

Builds a `pubtable` `ptModel` (one coefficient column) for a single
good/equation from a `quaidsOut` fit. Optional -- requires the `pubtable`
package.

## Format

```gauss
struct ptModel mdl;
mdl = ptModelFromQuaids(name, qOut, eqIdx);
```

## Parameters

- `name` (*string*) - column label (e.g. the good's name).
- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md).
- `eqIdx` (*scalar*) - 1-based index of the good/equation to extract.

## Returns

`mdl` is a `pubtable` `ptModel` structure: `estimates`/`stdErrors` from
`qOut.bestB[., eqIdx]`/the corresponding block of `qOut.bestV`, p-values,
and term names. Row order: intercept block (`qOut.xnam`) | price/gamma
block (`GAMMA_` + good name, one row per good) | `BETA_LX` |
`LAMBDA_LX2` (QUAIDS only) | one row per IV control-function residual term
(`qOut.unam`, `qOut.nu` rows -- always present, since
[quaidsFit](quaidsFit.md) always treats log total expenditure as
endogenous).

## Remarks

Uses `qOut.bestB`/`qOut.bestV` -- the same "most-constrained estimate
actually fit" that elasticities and the Slutzky diagnostic are evaluated
against. Requires `pubtable.sdf`'s `ptModel` struct type declared before
this file compiles (`#include pubtable.sdf`/`pubtable.src` or `library
pubtable;`), and errors with a clear stub message if `quaids.sdf` was not
included first (guarded by `#ifDef QUAIDS_SDF_INCLUDED`).

Not part of `package.json`'s `src` array -- see
[ptFromQuaids](ptFromQuaids.md) for why, and how to include this file.

## Examples

```gauss
library pubtable, quaids;
#include pubtable_quaids.src

struct ptModel mdl;
mdl = ptModelFromQuaids("Good 1", qOut, 1);
```

## Source

`pubtable_quaids.src`

## See Also

[ptFromQuaids](ptFromQuaids.md), [ptFromQuaidsFamily](ptFromQuaidsFamily.md)
