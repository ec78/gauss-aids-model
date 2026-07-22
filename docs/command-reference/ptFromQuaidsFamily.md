# ptFromQuaidsFamily

## Purpose

Convenience dispatcher: builds the appropriate `pubtable` `ptTable` for
either a `quaidsOut` or a `quaidsElasOut`, without the caller needing to
know which adapter proc to call. Optional -- requires the `pubtable`
package.

## Format

```gauss
struct ptTable tbl;
tbl = ptFromQuaidsFamily(x);
```

## Parameters

- `x` (*`quaidsOut` or `quaidsElasOut` structure*) - dispatches via
  `isStructType`.

## Returns

`tbl` is a `pubtable` `ptTable`:

- `quaidsOut` -> [ptFromQuaids](ptFromQuaids.md)
- `quaidsElasOut` -> [ptFromQuaidsElas](ptFromQuaidsElas.md)

Errors (`errorlog` + `end`) for any other struct type.

## Remarks

Mirrors `pubtable_qardl.src`'s `ptFromArdlFamily` dispatcher. Defined
unconditionally (no `#ifDef QUAIDS_SDF_INCLUDED` guard needed) since it
just delegates to the `(...)`-variadic adapters above, which already
handle the QUAIDS-missing case themselves.

## Examples

```gauss
struct ptTable coefTbl;
coefTbl = ptFromQuaidsFamily(qOut);        // dispatches to ptFromQuaids

struct ptTable elasTbl;
elasTbl = ptFromQuaidsFamily(elasOut);     // dispatches to ptFromQuaidsElas
```

## Source

`pubtable_quaids.src`

## See Also

[ptFromQuaids](ptFromQuaids.md), [ptFromQuaidsElas](ptFromQuaidsElas.md)
