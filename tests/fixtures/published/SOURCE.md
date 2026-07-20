# `blanciforti86_food32.csv` — source and license

**Data**: annual U.S. food consumption data, 1947–1978 (32 years), 4 food
groups. Columns: `pFood1`–`pFood4` (price indices), `wFood1`–`wFood4`
(budget shares), `xFood` (total food expenditure), `xAgg` (total aggregate
consumption expenditure across all 11 commodity groups in the original
source — used here as an instrument for `xFood`).

**Original source**: Blanciforti, L., Green, R., King, G. (1986). *U.S.
Consumer Behavior Over the Postwar Period: An Almost Ideal Demand System
Analysis*. Giannini Foundation Monograph No. 40, University of California.

**Immediate source**: the `Blanciforti86` dataset bundled with the R
package [`micEconAids`](https://cran.r-project.org/package=micEconAids)
(Arne Henningsen, GPL ≥ 2), rows 1–32, columns `pFood1`–`pFood4`,
`wFood1`–`wFood4`, `xFood`, `xAgg`. Extracted 2026-07-20 via:

```r
library(micEconAids); data("Blanciforti86")
B86 <- Blanciforti86[1:32, ]
write.csv(B86[, c("pFood1","pFood2","pFood3","pFood4",
                   "wFood1","wFood2","wFood3","wFood4","xFood","xAgg")], ...)
```

**License note**: `micEconAids` itself is GPL ≥ 2. The underlying data is
aggregate published U.S. consumption statistics from a 1986 academic
monograph; it is used here for non-commercial validation/replication
purposes, consistent with the package's own inclusion and documentation of
it as a standard teaching/replication dataset for AIDS demand-system
estimation (used in the package's own vignette). Committed to this
MIT-licensed repository with explicit approval from the repo owner
(2026-07-20) after the licensing question was raised — see
`GOLD_STANDARD_TODO.md` Milestone 3.

**Used by**: `tests/quaids_published_validation_test.e`.

**Reference scripts** (reproduce the R/Python numbers embedded in that
test): `generate_r_reference.R` (requires R + `micEconAids`) and
`python_reference_check.py` (requires `numpy`/`pandas`; documented as
supplementary, not the pass/fail source — see its header comment).
