# ptFromQuaids

## Purpose

Builds a `pubtable` `ptTable` coefficient comparison table -- one column
per good/equation -- from a `quaidsOut` fit, ready to export. Optional --
requires the `pubtable` package.

## Format

```gauss
struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);
```

## Parameters

- `qOut` (*`quaidsOut` structure*) - from [quaidsFit](quaidsFit.md).

## Returns

`coefTbl` is a `pubtable` `ptTable`: one comparison column per good, built
via `ptModelCompare` over [ptModelFromQuaids](ptModelFromQuaids.md) called
for each good. Title reflects the fitted model (`"LA-AIDS results"` /
`"AIDS results"` / `"QUAIDS results"`).

## Remarks

Mirrors `pubtable`'s own bundled `pubtable_qardl.src` adapter pattern
(`ptFromQardl`'s per-quantile comparison table), but lives inside this
repo (`src/pubtable_quaids.src`) rather than physically inside the
installed `pubtable` package -- see `GOLD_STANDARD_TODO.md`'s Milestone 6
section for why. **Not** listed in `package.json`'s `src` array: its
return-type annotation (`proc (struct ptTable) = ...`) needs
`pubtable.sdf`'s struct types declared unconditionally, which would make
`pubtable` a hard compile-time dependency for the whole `quaids` package
to compile. A caller `#include`s it directly, after `quaids.sdf`,
`quaids.src` and `pubtable.sdf`/`pubtable.src` (or `library quaids,
pubtable;`).

## Examples

```gauss
library pubtable, quaids;
#include pubtable_quaids.src

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);
call ptExport(coefTbl, "results.tex");
```

## Source

`pubtable_quaids.src`

## See Also

[ptModelFromQuaids](ptModelFromQuaids.md), [ptFromQuaidsFamily](ptFromQuaidsFamily.md)
