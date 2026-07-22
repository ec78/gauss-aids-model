# ptTablesFromQuaidsElas

## Purpose

Builds a 3-table `pubtable` bundle -- income elasticities, uncompensated
(Marshallian) price elasticities, compensated (Hicksian) price elasticities
-- from a `quaidsElasOut`, suitable for `ptExportAll()`. Optional --
requires the `pubtable` package.

## Format

```gauss
struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);
```

## Parameters

- `elasOut` (*`quaidsElasOut` structure*) - from
  [quaidsElasFit](quaidsElasFit.md).

## Returns

`elasTbls` is a 3x1 `struct ptTable` array:

1. Income elasticities (same as [ptFromQuaidsElas](ptFromQuaidsElas.md)).
2. Uncompensated (Marshallian) price elasticities: `n x n` matrix, one
   value row and one `(se)` row per good.
3. Compensated (Hicksian) price elasticities: same shape as (2).

The two price-elasticity tables are built directly as `ptTable`s (not via
`ptModel`, since an `n x n` matrix with paired value/SE rows doesn't fit
`ptModel`'s single-coefficient-vector design) via the internal
`_ptQuaidsElasMatrixTable()` helper.

## Remarks

Mirrors `pubtable_qardl.src`'s `*Full`-workflow table bundles (e.g.
`ptTablesFromQardlFull`).

## Examples

```gauss
struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);
call ptExport(elasTbls[1], "income_elasticities.md");
call ptExport(elasTbls[2], "uncompensated_elasticities.tex");
call ptExport(elasTbls[3], "compensated_elasticities.csv");
```

## Source

`pubtable_quaids.src`

## See Also

[ptFromQuaidsElas](ptFromQuaidsElas.md), [ptModelFromQuaidsElas](ptModelFromQuaidsElas.md)
