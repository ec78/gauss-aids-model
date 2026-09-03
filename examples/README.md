# QUAIDS Examples

A numbered, read-in-order suite of runnable GAUSS programs, one per major
feature area. Each is a focused, standalone script -- read it top to
bottom, then run it and compare the printed output. All 13 share one
small, well-commented synthetic dataset generator
([`example_data.src`](example_data.src)): a simulated household budget
survey with 5 spending categories (Food, Housing, Transportation,
Recreation, Other), a household-size demographic shifter, and log total
expenditure instrumented by a log wage-income variable. See that file's
own header comment for a real, honest caveat about this synthetic data
(individual budget shares can fall outside [0,1], even though they
always sum to exactly 1 -- explained there, and in
[`01_basic_estimation.e`](01_basic_estimation.e)).

Run each example from this directory:

```powershell
cd examples
tgauss -b -x 01_basic_estimation.e
```

(`10_curvature_imposition.e` and `13_pubtable_reporting.e` need the
optional `optmt`/`pubtable` packages installed and loaded first, noted
in their own header comments.)

## Reading order

| # | File | Demonstrates | Requires |
| --- | --- | --- | --- |
| 01 | [`01_basic_estimation.e`](01_basic_estimation.e) | `quaidsControlCreate`, `quaidsFit`/`printQuaids`, `quaids()`, the LA-AIDS/iterated-AIDS/QUAIDS model switch | -- |
| 02 | [`02_dataframe_input.e`](02_dataframe_input.e) | `quaidsFull()` -- selecting columns from a named-column dataframe instead of assembling matrices by hand | -- |
| 03 | [`03_preflight_diagnostics.e`](03_preflight_diagnostics.e) | `quaidsPreflight`/`printQuaidsPreflight` -- a warning, a second warning, and a hard failure | -- |
| 04 | [`04_hypothesis_tests.e`](04_hypothesis_tests.e) | `quaidsHomogeneityTest`, `quaidsJointTest`, `quaidsQuadraticTest` | -- |
| 05 | [`05_elasticities_shares_slutzky.e`](05_elasticities_shares_slutzky.e) | `quaidsElasFit`, `quaidsSharesFit`, `quaidsSlutzky` at the mean and a counterfactual price scenario | -- |
| 06 | [`06_welfare_analysis.e`](06_welfare_analysis.e) | `quaidsWelfareFit`/`printQuaidsWelfare` -- compensating/equivalent variation for a price change | -- |
| 07 | [`07_zero_share_correction.e`](07_zero_share_correction.e) | `quaidsZeroFit`/`printQuaidsZero` -- Shonkwiler-Yen correction for corner solutions, unconstrained and homogeneity+symmetry-constrained | -- |
| 08 | [`08_robust_standard_errors.e`](08_robust_standard_errors.e) | `quaidsRobustFit`/`quaidsRobustCovariance`, `quaidsRobustBootstrapFit` -- heteroskedasticity- and cluster-robust SE | -- |
| 09 | [`09_replicate_weights.e`](09_replicate_weights.e) | `quaidsReplicateWeightFit`/`printQuaidsReplicateWeight` -- a hand-built JK1 delete-one-cluster design | -- |
| 10 | [`10_curvature_imposition.e`](10_curvature_imposition.e) | `quaidsCurvatureFit` (AIDS and QUAIDS), `quaidsCurvatureBootstrapFit`, `quaidsCurvatureBootstrapCI` | `optmt` |
| 11 | [`11_survey_weighted_estimation.e`](11_survey_weighted_estimation.e) | `quaidsFit`'s `weight=` argument, `quaidsSurveyWorkflowFit` -- naive vs. weighted estimation on an informatively-sampled dataset | -- |
| 12 | [`12_applied_workflow.e`](12_applied_workflow.e) | `quaidsWorkflowFit`/`quaidsWorkflowScenarioFit` -- preflight, fit, shares/elasticities, robust SE, and a CV/EV scenario in one call | -- |
| 13 | [`13_pubtable_reporting.e`](13_pubtable_reporting.e) | `ptFromQuaids`, `ptTablesFromQuaidsElas`, `ptTablesFromQuaidsWorkflow` -- LaTeX/Markdown/CSV export | `pubtable` |

See [docs/USAGE_GUIDE.md](../docs/USAGE_GUIDE.md) for the prose reference
each example pairs with, and
[docs/COMMAND_REFERENCE.md](../docs/COMMAND_REFERENCE.md) for the full
per-procedure documentation.
