# ptFromQuaidsElas

## Purpose

Builds the single canonical `pubtable` `ptTable` for a `quaidsElasOut` --
the income elasticities table. Optional -- requires the `pubtable`
package.

## Format

```gauss
struct ptTable elasIncomeTbl;
elasIncomeTbl = ptFromQuaidsElas(elasOut);
```

## Parameters

- `elasOut` (*`quaidsElasOut` structure*) - from
  [quaidsElasFit](quaidsElasFit.md).

## Returns

`elasIncomeTbl` is a `pubtable` `ptTable`, titled `"Income elasticities"`.

## Remarks

For the full report including the `n x n` uncompensated/compensated price
elasticity matrices, use
[ptTablesFromQuaidsElas](ptTablesFromQuaidsElas.md) instead, which returns
a 3-table bundle.

## Examples

```gauss
struct ptTable elasIncomeTbl;
elasIncomeTbl = ptFromQuaidsElas(elasOut);
call ptExport(elasIncomeTbl, "income_elasticities.md");
```

## Source

`pubtable_quaids.src`

## See Also

[ptModelFromQuaidsElas](ptModelFromQuaidsElas.md), [ptTablesFromQuaidsElas](ptTablesFromQuaidsElas.md),
[ptFromQuaidsFamily](ptFromQuaidsFamily.md)
