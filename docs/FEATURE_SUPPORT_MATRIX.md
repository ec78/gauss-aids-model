# Feature Support Matrix

This matrix summarizes the current public support surface across the three
model choices `quaidsFit()` (and its wrappers `quaids()`/`quaidsFull()`)
selects via `aCtl.linear`/`aCtl.maxiter` -- see the
[usage guide](USAGE_GUIDE.md#choosing-a-model-la-aids-vs-iterated-aids-vs-quaids)
for the exact switch values.

| Feature | LA-AIDS | Iterated AIDS | QUAIDS |
| --- | --- | --- | --- |
| Estimator | `quaidsFit`/`quaids`/`quaidsFull` | `quaidsFit`/`quaids`/`quaidsFull` | `quaidsFit`/`quaids`/`quaidsFull` |
| Price index | Stone (one-step) | Nonlinear translog, iterated | Nonlinear translog, iterated |
| Quadratic log-expenditure term | No | No | Yes |
| IV treatment of total expenditure | Always (control-function) | Always (control-function) | Always (control-function) |
| Overidentification test | Yes, if `ninst > nu` | Yes, if `ninst > nu` | Yes, if `ninst > nu` |
| Homogeneity imposition | Yes (`aCtl.homogenous=1`) | Yes (`aCtl.homogenous=1`) | Yes (`aCtl.homogenous=1`) |
| Symmetry imposition (given homogeneity) | Yes | Yes | Yes |
| Symmetry-given-homogeneity test | Yes (built into `quaidsFit`) | Yes (built into `quaidsFit`) | Yes (built into `quaidsFit`) |
| Standalone homogeneity test | Yes (`quaidsHomogeneityTest`, needs `aCtl.homogenous=0` fit) | Yes | Yes |
| Standalone joint homogeneity+symmetry test | Yes (`quaidsJointTest`) | Yes | Yes |
| Elasticities at arbitrary points | Yes (`quaidsElasFit`) | Yes | Yes |
| Delta-method elasticity standard errors | Yes | Yes | Yes |
| Exact algebraic identity validation (Engel/Cournot/homogeneity) | Yes | Yes | Yes |
| Slutzky negativity diagnostic | Yes (`quaidsSlutzky`) | Yes | Yes |
| Curvature imposition | No (diagnosis only) | No (diagnosis only) | No (diagnosis only) |
| Dataframe/column-name entry point | Yes (`quaidsFull`) | Yes | Yes |
| Formula-string (`"y ~ x"`) API | Not applicable (multi-equation system) | Not applicable | Not applicable |
| `pubtable` export (LaTeX/Markdown/CSV/...) | Yes (`src/pubtable_quaids.src`, optional) | Yes | Yes |
| Synthetic deterministic validation | Yes (`tests/quaids_synthetic_validation_test.e`) | Yes | Yes |
| Published-data cross-validation vs. R | Yes (`Blanciforti86` vs. 3SLS, `tests/quaids_published_validation_test.e`) | Yes (`Blanciforti86` vs. `method="IL"`, wider tolerance -- see Notes) | No independent reference implementation exists (see Notes) |
| Iteration convergence guarantee | Not applicable (one-step) | No -- roughly half of random seeds in one synthetic-DGP family failed to converge cleanly; check `qOut.converged` | No -- same caveat |
| Installed-package (`library quaids;`) support | Yes | Yes | Yes |

## Notes

- "Always (control-function)" means `instr` is a required argument to
  every estimator entry point -- there is no exogenous-total-expenditure
  estimation mode in this library.
- The published-data cross-validation (`Blanciforti86` vs. R's
  `micEconAids`) covers both LA-AIDS (`aCtl.linear=1, aCtl.maxiter=1`, vs.
  `aidsEst(..., instNames=...)`, 3SLS -- max abs difference ~0.021) and
  iterated AIDS (`aCtl.linear=1, aCtl.maxiter>1`, vs.
  `aidsEst(method="IL", ...)`, the Iterated Linear Least Squares Estimator
  -- max abs difference ~0.11, tolerance `0.15`). The iterated-AIDS
  comparison has a wider gap for a real, understood reason, not
  approximation slop: `micEconAids`'s `method="IL"` does not support
  instrumental variables (combining it with `instNames` segfaults R's
  `aidsEst` rather than erroring cleanly -- confirmed by direct testing),
  so that reference is SUR-estimated, while GAUSS's iterated fit always
  instruments log total expenditure. The comparison therefore spans both a
  different estimation algorithm *and* an IV-vs-no-IV difference. **QUAIDS
  has no independent reference implementation available**: `micEconAids`
  does not implement a quadratic log-expenditure term at all, and no other
  comparably-established QUAIDS implementation was found (see
  `GOLD_STANDARD_TODO.md`'s Milestone 3 section on the Python from-scratch
  replica, kept as supplementary evidence only). QUAIDS's validation is
  therefore the known-true synthetic-DGP recovery in
  `tests/quaids_synthetic_validation_test.e` -- a real, non-circular check
  (independently-generated data with known-true parameters, not just
  re-running the estimator on its own output), but a different, weaker
  tier of evidence than cross-implementation agreement on real published
  data. Documented here rather than silently claimed as equivalent.
- Curvature imposition (Diewert-Wales Cholesky reparametrization or
  similar) is scoped and explicitly deferred, not silently absent -- see
  [Methodology Notes](METHODOLOGY_NOTES.md#slutzky-negativity) and
  `GOLD_STANDARD_TODO.md`'s Milestone 5 section. Even the R `micEconAids`
  reference implementation used for this library's own cross-validation
  only offers curvature diagnosis, not imposition.

Related documentation:

- [Usage guide](USAGE_GUIDE.md)
- [Methodology notes](METHODOLOGY_NOTES.md)
- [Command reference](COMMAND_REFERENCE.md)
