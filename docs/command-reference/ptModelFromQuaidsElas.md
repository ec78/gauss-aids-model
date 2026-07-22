# ptModelFromQuaidsElas

## Purpose

Builds a `pubtable` `ptModel` for income elasticities (one column) from a
`quaidsElasOut`. Optional -- requires the `pubtable` package.

## Format

```gauss
struct ptModel mdl;
mdl = ptModelFromQuaidsElas(name, elasOut);
```

## Parameters

- `name` (*string*) - column label.
- `elasOut` (*`quaidsElasOut` structure*) - from
  [quaidsElasFit](quaidsElasFit.md).

## Returns

`mdl` is a `pubtable` `ptModel`: `estimates`/`stdErrors` from
`elasOut.er`/`elasOut.ser`, one row per good, named from `elasOut.wnam`.

## Remarks

Covers income elasticities only -- the `n x n` uncompensated/compensated
price elasticity matrices don't fit `ptModel`'s single-coefficient-vector
shape; see [ptTablesFromQuaidsElas](ptTablesFromQuaidsElas.md) for those.

## Examples

```gauss
struct ptModel mdl;
mdl = ptModelFromQuaidsElas("Income elasticities", elasOut);
```

## Source

`pubtable_quaids.src`

## See Also

[ptFromQuaidsElas](ptFromQuaidsElas.md), [ptTablesFromQuaidsElas](ptTablesFromQuaidsElas.md)
