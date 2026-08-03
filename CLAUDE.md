# CLAUDE.md — GAUSS QUAIDS Library

Context file for Claude Code sessions working on this repository.

## What this library does

Estimates **Almost Ideal Demand System** models: linearized AIDS (Stone price
index), iterated AIDS (nonlinear translog price index), and **QUAIDS**
(Banks, Blundell & Lewbel 1997 quadratic-log-expenditure extension), with
instrumental-variables treatment of (endogenous) log total expenditure.
Estimation imposes homogeneity and/or Slutzky symmetry via
iterated FGLS with cross-equation restrictions applied through a
minimum-distance reparametrization. Use cases: consumer demand estimation,
welfare analysis, elasticity calculation, testing demand-theory restrictions.

The library is **pre-alpha** (package version `0.21.0`) and is packaged as an
installable GAUSS application package (`library quaids;`). See
`GOLD_STANDARD_TODO.md` for the full roadmap — this file is the
quick-orientation companion to it, and should be kept synchronized with it.

**Naming**: the package and its public procs use a `quaids` prefix (decided
at Milestone 0), even though the estimator also covers plain linear AIDS —
QUAIDS is the more general model actually implemented, and `aids` was judged
too likely to collide/confuse as a bare identifier. "AIDS"/"Almost Ideal
Demand System" remains the correct term for the model family in docs, papers,
and comments; only the GAUSS identifier prefix changed.

## Repository layout (post-Milestone-26)

```
src/
  quaids.sdf        # Struct definitions: quaidsControl, quaidsOut,
                    #   quaidsElasOut. Guarded by #ifndef/#define
                    #   QUAIDS_SDF_INCLUDED / #endif (added Milestone 6) so
                    #   pubtable_quaids.src has a symbol to detect.
  quaidsutil.src    # quaidsControlCreate() / getDefaultQuaidsControl().
  quaidsiv.src      # _quaidsIVFirstStage() -- private helper, the
                    #   instrumental-variables first-stage regression of log
                    #   total expenditure. The one internal phase of
                    #   quaidsFit() that was cleanly separable; see the
                    #   Milestone 2 scoping note below for why the rest of
                    #   quaidsFit() (starting values, iteration, variance,
                    #   overidentification test, symmetry test) was not
                    #   further split. Milestone 26 adds a required `weight`
                    #   argument (all three callers -- quaidsFit(),
                    #   quaidsPreflight(), quaidsZeroFit() -- updated).
  quaidszerocorrect.src # Milestone 19: quaidsZeroFit()/printQuaidsZero()
                    #   -- AIDS/QUAIDS estimation corrected for zero
                    #   budget shares via the Shonkwiler-Yen (1999)
                    #   two-step procedure. Divides the whole share
                    #   equation by the first-stage probit's fitted
                    #   probability F_i (rather than rescaling every
                    #   regressor, which would break the shared-design-
                    #   matrix Kronecker-product identity every other
                    #   estimation stage relies on), then imposes a
                    #   diagonal-delta minimum-distance restriction
                    #   mirroring quaids.src's own symmetry-restriction
                    #   stage. Unconstrained only (errors if
                    #   aCtl.homogenous=1). Uses GAUSS's built-in glm()
                    #   for the per-good first-stage probit -- no new
                    #   package dependency. See "Milestone 19: zero
                    #   budget share correction (Shonkwiler-Yen)" below.
  quaidselas.src    # quaidsElas_() (silent, low-level), quaidsElasFit()
                    #   (silent, struct-returning: point estimates +
                    #   standard errors), printQuaidsElas() (the separated
                    #   printer), quaidsElas() (backward-compatible
                    #   wrapper: fit then print) -- elasticities at a
                    #   point. See "Milestone 5: elasticities
                    #   generalization" below.
  quaidsshares.src  # Milestone 16: quaidsSharesFit()/printQuaidsShares()
                    #   -- the model-implied predicted budget share vector
                    #   (with full delta-method covariance) at an
                    #   arbitrary point. A deliberately independent third
                    #   copy of the share formula quaidsElas_() already
                    #   computes internally (not a refactor of it -- see
                    #   "Milestone 16" below). No new package dependency.
  quaidsslutzky.src # quaidsSlutzky() -- Slutzky negativity diagnostic.
  quaids.src        # quaidsFit() (silent, struct-returning estimation core;
                    #   calls _quaidsIVFirstStage()), printQuaids() (the
                    #   separated estimation-report printer), quaids()
                    #   (backward-compatible wrapper: fits, prints,
                    #   reproduces the legacy elasticities/descriptive
                    #   stats/Slutzky report via quaidsElas()/
                    #   quaidsSlutzky(), returns the 4 legacy matrices).
                    #   Milestone 26 adds an OPTIONAL trailing `weight`
                    #   argument (GAUSS dynargs) -- a genuine sampling-
                    #   weighted point estimate via the standard
                    #   sqrt(weight)-scaled cross-product WLS trick, an
                    #   exact no-op when omitted or uniform.
  quaidsformula.src # quaidsFull() -- dataframe/column-name entry point;
                    #   selects w/intcpt/prices/totexp/instr from a
                    #   dataframe by column name and calls quaidsFit().
  quaidstests.src   # quaidsHomogeneityTest(), quaidsJointTest() -- standalone
                    #   Wald tests, operating on an already-computed
                    #   unconstrained (aCtl.homogenous=0) quaidsOut. See
                    #   "Milestone 4: new hypothesis tests" below.
                    #   Milestone 17 adds quaidsQuadraticTest() -- a Wald
                    #   test of whether the QUAIDS quadratic term is
                    #   needed, additionally requiring qOut.linear==0
                    #   (an AIDS fit never estimates lambda, so there is
                    #   nothing to test). See "Milestone 17" below.
  quaidscurvature.src # Milestone 10: quaidsCurvatureFit()/
                    #   printQuaidsCurvature() -- local curvature (Slutzky
                    #   negative semidefiniteness) imposition for LA-AIDS/
                    #   AIDS via the Diewert-Wales Cholesky
                    #   reparametrization, at the sample mean. Hard
                    #   compile-time dependency on the installed `optmt`
                    #   package (struct PV/optmtControl/optmtResults,
                    #   pvPacki/pvUnpack/optmt) -- unlike
                    #   pubtable_quaids.src, this IS listed in
                    #   package.json's src array (required public API),
                    #   so `optmt` is now a real package dependency
                    #   (package.json's deps array). See "Milestone 10:
                    #   curvature imposition" below. Milestone 15 adds
                    #   quaidsCurvatureBootstrapFit()/
                    #   printQuaidsCurvatureBootstrap() -- an i.i.d. row
                    #   bootstrap standard error for quaidsCurvatureFit(),
                    #   reported alongside (not replacing) its delta-method
                    #   SE. No new package dependency. See "Milestone 15:
                    #   bootstrap standard errors" below. Milestone 18 adds
                    #   quaidsCurvatureBootstrapCI() -- percentile CIs from
                    #   the bootstrap's raw draws, no new resampling -- and
                    #   fixes a real row-major-vs-column-major reshape()
                    #   bug (present since Milestone 10/15) in both
                    #   quaidsCurvatureFit()'s se and
                    #   quaidsCurvatureBootstrapFit()'s seBoot. See
                    #   "Milestone 18: percentile bootstrap confidence
                    #   intervals" below.
  quaidswelfare.src # Milestone 11: quaidsWelfareFit()/printQuaidsWelfare()
                    #   -- exact compensating/equivalent variation for a
                    #   price change, holding nominal expenditure fixed.
                    #   Unified across LA-AIDS/iterated-AIDS/QUAIDS (no
                    #   scoping-out needed, unlike Milestone 10), pure
                    #   closed-form algebra, no new package dependency.
                    #   See "Milestone 11: welfare measures" below.
  quaidsrobust.src  # Milestone 20: quaidsRobustFit()/printQuaidsRobust()
                    #   -- heteroskedasticity-robust and cluster-robust
                    #   standard errors for an already-fitted quaidsFit()
                    #   result, generalizing the pooled, homoskedastic
                    #   S.*.inv(gg) sandwich every other covariance in
                    #   this library uses to a per-observation or
                    #   per-cluster score aggregation, unified through
                    #   one clusterId argument (0 = robust, a Tx1 group
                    #   vector = cluster-robust). Also
                    #   quaidsRobustBootstrapFit()/
                    #   printQuaidsRobustBootstrap() -- a cluster-aware
                    #   nonparametric bootstrap alternative. Uses a
                    #   SIMPLIFIED bread (inv(gg)-based, not quaidsFit()'s
                    #   own nonlinear-feedback-corrected Jacobian), found
                    #   empirically to make its se dramatically more
                    #   conservative than qOut's own classical SE -- an
                    #   expected consequence of the simplification, not a
                    #   bug. Milestone 22 adds quaidsRobustCovariance()
                    #   and quaidsRobustBootstrapCovariance() to expand
                    #   reduced robust/bootstrap covariance into
                    #   qOut.bestB's full basis for shares/elasticities/
                    #   welfare delta-method SE. See "Milestone 20" and
                    #   "Milestone 22" below. Milestone 26 adds the same
                    #   OPTIONAL trailing `weight` argument to both procs
                    #   -- a DIFFERENT scaling convention from quaidsFit()'s
                    #   own sqrt(weight): the bread keeps sqrt(weight), but
                    #   the per-observation score contribution uses PLAIN
                    #   weight (Horvitz-Thompson pweight-robust
                    #   convention). See "Milestone 26" below.
  quaidsdiagnostics.src # Milestone 23: quaidsPreflight()/
                    #   printQuaidsPreflight() -- silent, estimator-free
                    #   preflight diagnostics for dimensions, finite data,
                    #   share adding-up, zero/negative shares, variation,
                    #   first-stage IV strength, design invertibility,
                    #   cluster counts, and convergence-risk screening.
                    #   Milestone 26 adds a required `weight` argument
                    #   (mirroring clusterId's own required-positional
                    #   convention here) plus weightValid/weightSum/effN
                    #   fields.
  quaidsworkflow.src # Milestone 21 seed: quaidsWorkflowFit(), a thin
                    #   applied workflow layer composing quaidsFit(),
                    #   quaidsSharesFit(), quaidsElasFit(), and
                    #   quaidsRobustFit()/quaidsRobustCovariance() into
                    #   one silent struct-returning call. Milestone 24
                    #   also calls quaidsPreflight() first and echoes a
                    #   compact non-gating diagnostic summary in the flat
                    #   workflow struct. Plus
                    #   quaidsWorkflowScenarioFit(), which fills the same
                    #   struct's welfare fields and robust welfare SE for
                    #   one explicit CV/EV price-change scenario. This is
                    #   not a new estimator; it is the first public
                    #   one-call workflow bundle for applied scripts.
                    #   Milestone 26 adds an OPTIONAL trailing `weight`
                    #   argument to both procs, threaded into their own
                    #   quaidsFit()/quaidsPreflight()/quaidsRobustFit()
                    #   calls.
  quaidssurvey.src  # Milestone 25 seed: quaidsSurveyWorkflowFit(), an
                    #   opt-in survey/microdata workflow wrapper that
                    #   recomputes the workflow evaluation point as a
                    #   sampling-weighted mean and then recomputes
                    #   shares/elasticities there. Milestone 26: its own
                    #   `weight` argument now ALSO fits the estimator
                    #   (forwarded into quaidsWorkflowFit()'s new
                    #   argument) -- a deliberate behavior change from the
                    #   Milestone 25 release, which left the estimator
                    #   unweighted. See "Milestone 26" below.
  pubtable_quaids.src # Optional pubtable adapter -- ptModelFromQuaids()/
                    #   ptFromQuaids() (coefficient tables),
                    #   ptModelFromQuaidsElas()/ptFromQuaidsElas()/
                    #   ptTablesFromQuaidsElas() (elasticity tables),
                    #   ptTablesFromQuaidsWorkflow() (Milestone 21
                    #   applied workflow table bundle),
                    #   ptFromQuaidsFamily() dispatcher. NOT in
                    #   package.json's src array or self-included by any
                    #   other src/ file -- has a hard compile-time
                    #   dependency on pubtable.sdf's ptModel/ptTable struct
                    #   types, so a caller must #include pubtable.sdf/
                    #   pubtable.src (or `library pubtable;`) before this
                    #   file, same as quaids.sdf/quaids.src. See
                    #   "Milestone 6: reporting via pubtable" below.
examples/
  quaids_example.e  # One synthetic 5-good dataset (homogeneity/symmetry true
                    #   by construction), run through quaids() with
                    #   eyeballed comparison of printed estimates to true
                    #   parameters. Not an automated test — no assertions.
                    #   Uses `library quaids;` against the installed
                    #   package (Milestone 9 -- previously #included the
                    #   source tree directly via ../src/..., from before
                    #   Milestone 7 made the package installable; switched
                    #   once README.md/USAGE_GUIDE.md started documenting
                    #   `library quaids;` as the primary usage pattern, so
                    #   the examples actually demonstrate what the docs
                    #   promise, matching gauss-qardl's own examples/,
                    #   which are all `library qardl;`-based). Source-tree
                    #   #include-based testing still lives on in tests/
                    #   (except tests/package_public_api.e, which is also
                    #   library-based, by design -- see "Testing status"
                    #   below).
  workflow_example.e # Milestone 21: installed-package one-call applied
                    #   workflow example using quaidsWorkflowFit() and
                    #   quaidsWorkflowScenarioFit() for mean-point shares/
                    #   elasticities/robust SE and a CV/EV scenario.
  pubtable_export_example.e # Milestone 6: same style, but exports a
                    #   quaidsFit() coefficient table and a
                    #   quaidsElasFit() elasticity report to
                    #   LaTeX/Markdown/CSV via pubtable_quaids.src.
                    #   Requires the pubtable package installed. Uses
                    #   `library quaids, pubtable;` (Milestone 9) plus a
                    #   bare `#include quaids.sdf` -- required so that
                    #   pubtable_quaids.src's #ifDef QUAIDS_SDF_INCLUDED
                    #   guard is active; `library quaids;` alone lazily
                    #   loads procs on demand and does not run quaids.sdf's
                    #   #define (confirmed empirically), matching how
                    #   pubtable's own bundled pubtable_qardl.src documents
                    #   the identical requirement for qardl.sdf. Still
                    #   #includes ../src/pubtable_quaids.src by relative
                    #   path, since that adapter is not part of the
                    #   installed quaids package (see "Milestone 6" below).
tests/
  quaids_schema_test.e         # Milestone 1: asserts quaidsOut field
                    #   values/shapes, that quaidsFit() prints nothing, and
                    #   that the legacy quaids() wrapper's returned matrices
                    #   are byte-identical to the struct fields they're
                    #   drawn from.
  quaids_formula_parity_test.e # Milestone 2: builds the same synthetic
                    #   dataset as both plain matrices and a named-column
                    #   dataframe, and asserts quaidsFull(dataframe, ...)
                    #   matches quaidsFit(matrices...) exactly, including
                    #   the extraVars == 0 ("no extra intercept shifters")
                    #   path.
  quaidsfixtures.src            # Milestone 3: _quaidsSyntheticDGP(), a
                    #   shared 5-good synthetic-data generator (not part of
                    #   the public src/ API) parameterized by whether the
                    #   true model has a quadratic term and whether total
                    #   expenditure is genuinely endogenous. Returns the
                    #   true parameters pre-stacked in qOut.bS's row order.
                    #   Milestone 10 adds _quaidsCurvatureSyntheticDGP(): a
                    #   dedicated 5-good LA-AIDS dataset whose true gamma
                    #   is curvature-consistent at its own actual sample
                    #   mean by construction (built via a self-consistent
                    #   fixed-point iteration, seed=500 found by direct
                    #   screening -- see "Milestone 10" below). Milestone
                    #   19 adds _quaidsZeroSyntheticDGP(): a 5-good QUAIDS
                    #   dataset with structural coefficients scaled to
                    #   0.25x _quaidsSyntheticDGP()'s own magnitudes
                    #   (screened, not guessed -- the deterministic
                    #   price*gamma/expenditure*beta swing, not noise or
                    #   price variance, is the dominant driver of negative
                    #   latent shares at the original scale), censored the
                    #   economically correct way (an accounting identity,
                    #   not an ad hoc redistribution) to produce a genuine,
                    #   uneven-but-non-degenerate zero-share pattern at
                    #   seed=1 -- see "Milestone 19" below. Milestone 20
                    #   adds _quaidsClusterSyntheticDGP(): a 5-good AIDS
                    #   dataset with a genuine within-cluster-correlated
                    #   noise component (idiosyncratic noise plus a shock
                    #   shared by every row in the same cluster) -- needed
                    #   because every other fixture in this file has
                    #   plain i.i.d. noise, which cannot distinguish a
                    #   correct cluster-robust SE from an incorrect one
                    #   that ignores clustering entirely. Milestone 26
                    #   adds _quaidsSurveyWeightedDGP(): a 5-good AIDS
                    #   population subjected to informative two-stratum
                    #   sampling on an UNOBSERVED error term (not an
                    #   included covariate -- selecting on an included
                    #   regressor was tried first and found not to bias
                    #   naive estimation, the standard Heckman result),
                    #   returning the sampled data, each row's Horvitz-
                    #   Thompson weight, and the population's true
                    #   parameters -- see "Milestone 26" below.
  quaids_synthetic_validation_test.e # Milestone 3: 22 checks recovering
                    #   true DGP parameters within a documented tolerance
                    #   across all 6 (LA-AIDS/AIDS/QUAIDS) x (with/without
                    #   IV) combinations. See the "Seed sensitivity" note
                    #   below for why this uses one specific, documented
                    #   seed rather than an arbitrary one.
  quaids_published_validation_test.e # Milestone 3: fits quaidsFit() on
                    #   real published data (Blanciforti86 food-demand
                    #   data) and checks against an independent R
                    #   reference. This is the test that caught the
                    #   Stone-index starting-value bug -- see "Milestone 3:
                    #   real bug found and fixed" below. Milestone 9:
                    #   extended to also validate iterated AIDS against R
                    #   aidsEst(method="IL", ...) -- 19 checks total (was
                    #   11) -- see "Milestone 9" below.
  fixtures/published/
    blanciforti86_food32.csv  # The data itself; see SOURCE.md for
                    #   attribution/license.
    SOURCE.md                 # Attribution, license note, repo-owner
                    #   approval record.
    generate_r_reference.R    # Reproduces the R numbers hardcoded into
                    #   quaids_published_validation_test.e. Requires R +
                    #   the CRAN package micEconAids.
    python_reference_check.py # Independent from-scratch Python replica;
                    #   supplementary evidence, not the assertion source
                    #   (see its header for why). Requires numpy/pandas.
  quaids_hypothesis_tests_test.e # Milestone 4: 19 checks -- size AND power
                    #   for quaidsHomogeneityTest()/quaidsJointTest(), a
                    #   power check for the existing symmetry-given-
                    #   homogeneity test, and the first-ever exercise of
                    #   the overidentification test (every prior fixture
                    #   was exactly identified, ninst==nu). Milestone 17
                    #   adds 3 checks (22 total) for quaidsQuadraticTest()
                    #   -- size and power, reusing the file's existing
                    #   quadratic=1 fixture for power and one fresh
                    #   quadratic=0 fixture for size.
  quaids_elasticities_test.e   # Milestone 5: 17 checks -- parity between
                    #   quaidsElasFit()/printQuaidsElas() and the
                    #   pre-split quaidsElas_(), plus three EXACT
                    #   algebraic identities (Engel aggregation, Cournot
                    #   aggregation, elasticity homogeneity) checked at a
                    #   real out-of-sample observation and a synthetic
                    #   counterfactual price scenario -- not just the four
                    #   points quaids() has always used. Milestone 16
                    #   replaced this file's own private modelShareAt()
                    #   helper with a direct call to the new
                    #   quaidsSharesFit(), removing that duplication.
  quaids_shares_test.e         # Milestone 16: 21 checks -- quaidsSharesFit()'s
                    #   point estimate matches an independent, freshly
                    #   hand-evaluated share formula (both AIDS and QUAIDS
                    #   fixtures); exact adding-up (sum(w)==1); shape/
                    #   finiteness/non-negativity of se/v and se==sqrt(diag(v));
                    #   a shifted evaluation point gives a genuinely
                    #   different share (non-vacuous).
  quaids_pubtable_test.e       # Milestone 6/Milestone 21: pubtable adapter
                    #   checks -- exact numeric
                    #   parity between pubtable ptModel.estimates/
                    #   stdErrors and the qOut.bestB/qOut.bestV/
                    #   elasOut.er/ser values they're built from, shape/
                    #   title checks, the ptFromQuaidsFamily dispatcher,
                    #   ptTablesFromQuaidsWorkflow() workflow bundle
                    #   checks,
                    #   and an end-to-end export smoke test that writes
                    #   real .tex/.md/.csv files and reads them back.
                    #   Requires the pubtable package installed.
  quaids_curvature_test.e      # Milestone 10: 17 checks (AIDS block) --
                    #   recovery against a known-curvature-consistent true
                    #   gamma, exact adding-up/homogeneity/symmetry,
                    #   near-exact negative-semidefiniteness at the
                    #   reference point, and a non-vacuousness check (the
                    #   unconstrained fit genuinely violates curvature on
                    #   this fixture). Milestone 13 added a 14-check QUAIDS
                    #   block (31 total) -- convergence/exact NSD/non-
                    #   vacuousness/shape against the general QUAIDS
                    #   fixture, deliberately NOT true-gamma recovery (no
                    #   curvature-consistent QUAIDS fixture with plausible
                    #   shares was found despite a broad screen). Requires
                    #   the optmt package installed. See "Milestone 10:
                    #   curvature imposition" and "Milestone 13: QUAIDS
                    #   curvature imposition" below. Milestone 18 added 2
                    #   checks (33 total) confirming cOut.se's individual
                    #   cells are correctly positioned relative to cOut.b
                    #   (a regression guard for the reshape() bug found and
                    #   fixed that milestone).
  quaids_curvature_bootstrap_test.e  # Milestone 15: 26 checks -- bootstrap
                    #   run bookkeeping (requested/completed/failed/
                    #   attempts), shape/finiteness of the bootstrap SE,
                    #   exact echo of the base point estimate/delta-method
                    #   SE, and a plausibility check that the bootstrap SE
                    #   stays well-behaved where the delta-method SE does
                    #   not -- on both an AIDS (B=15) and a QUAIDS (B=5)
                    #   fixture, small B chosen to bound this file's own
                    #   runtime (~45-50s) given the ~0.9s/~7.3s per-
                    #   replication cost. NOT run by run_source_tests.ps1's
                    #   default (unflagged) invocation -- see
                    #   -SkipBootstrap below. Requires the optmt package
                    #   installed. See "Milestone 15: bootstrap standard
                    #   errors" below. Milestone 18 added 11 checks (37
                    #   total) -- quaidsCurvatureBootstrapCI() shape/
                    #   ordering/containment checks, a direct quantile()
                    #   cross-check at a specific flattened index, and a
                    #   seBoot position-correctness regression guard. See
                    #   "Milestone 18: percentile bootstrap confidence
                    #   intervals" below.
  quaids_welfare_test.e        # Milestone 11: 20 checks -- exact zero-
                    #   price-change identity, exact round-trip inverse-
                    #   function identity (feeding e(p1,u0) back into
                    #   V(.,p1) returns u0 exactly), first-order
                    #   (Marshallian) limiting-case consistency, and CV/EV
                    #   sign agreement -- on both a QUAIDS and an AIDS fit.
                    #   No extra package required. See "Milestone 11:
                    #   welfare measures" below.
                    #   All test/fixture files run from tests/ (or
                    #   tests/fixtures/published/ for the R/Python
                    #   scripts) as the working directory.
  quaids_convergence_sweep.e   # Milestone 12: a real, committed 200-seed
                    #   x 2-model convergence-reliability diagnostic for
                    #   the iterated estimator, replacing an informal
                    #   8-seed probe referenced since Milestone 3 that
                    #   never survived as a repo artifact (confirmed by a
                    #   full git-history search). Prints a per-seed row
                    #   plus a three-bucket summary (never-converged /
                    #   converged-but-wrong / converged-correctly) --
                    #   deliberately a diagnostic report generator, not a
                    #   pass/fail gate (no "ALL N CHECKS PASSED" line, not
                    #   in run_source_tests.ps1). Run via
                    #   run_convergence_sweep.ps1. See "Milestone 12:
                    #   numerical reliability" below.
  quaids_reliability_regression_test.e  # Milestone 12: 8 checks --
                    #   regression guard for the three changes the sweep
                    #   above led to: aCtl.relax=1 (default) is byte-
                    #   identical to leaving it unset; a previously-
                    #   crashing seed (an unguarded invpd() in the
                    #   symmetry-test block) now degrades gracefully
                    #   instead of aborting the call; aCtl.relax=.75 is
                    #   pinned against a real seed where it measurably
                    #   changes a never-converged fit into a correct one.
  quaids_zero_test.e  # Milestone 19: 17 checks -- the fixture's own
                    #   adding-up identity (exact, by construction); the
                    #   diagonal-delta restriction holds exactly
                    #   (off-diagonal hazard-coefficient entries exactly
                    #   0, on-diagonal entries genuinely estimated);
                    #   shape/finiteness of probitB/se/b; all n first-
                    #   stage probits converged; and the core validation
                    #   -- quaidsZeroFit()'s corrected coefficients
                    #   recover the true latent (uncensored) DGP
                    #   parameters better than naively fitting
                    #   quaidsFit() on the same censored data, on both a
                    #   max- and mean-absolute-difference basis. See
                    #   "Milestone 19: zero budget share correction
                    #   (Shonkwiler-Yen)" below.
  quaids_robust_test.e  # Milestone 20/22: 26 checks -- the point estimate
                    #   matches a fresh, independent hand-evaluation of
                    #   the sandwich formula; an exact-identity regression
                    #   guard (clusterId=0 vs. an explicit
                    #   seqa(1,1,nobs) per-row label give byte-identical
                    #   output); the reshape/cell-position regression
                    #   guard (written from day one, not found the hard
                    #   way a third time); shape/finiteness/non-
                    #   negativity; and the core non-vacuous check --
                    #   cluster-robust se is measurably larger than the
                    #   naive se on _quaidsClusterSyntheticDGP's genuinely
                    #   clustered data; and Milestone 22 full-basis
                    #   covariance expansion drives shares/elasticities/
                    #   welfare robust delta-method SE without changing
                    #   point estimates. See "Milestone 20" and
                    #   "Milestone 22" below.
  quaids_robust_bootstrap_test.e  # Milestone 20/22: 17 checks -- bootstrap
                    #   run bookkeeping, shape/finiteness, the reshape
                    #   regression guard for seBoot, exact echo of the
                    #   base point estimate/seRobust, full-basis
                    #   covariance expansion for downstream post-
                    #   estimation, and a plausibility
                    #   check that a cluster-aware bootstrap's seBoot
                    #   exceeds a plain-row bootstrap's on the same
                    #   genuinely clustered data. Not run by
                    #   run_source_tests.ps1's default invocation -- added
                    #   to the same -SkipBootstrap-gated group as
                    #   quaids_curvature_bootstrap_test.e.
  quaids_preflight_test.e  # Milestone 23: 16 checks (13 + 3 at Milestone
                    #   26) -- clean preflight, zero-share warning,
                    #   negative-share hard failure, adding-up failure,
                    #   cluster summaries, low price variation warning,
                    #   dimension-mismatch return, and (Milestone 26) a
                    #   genuinely unequal valid weight's effN/weightSum and
                    #   a non-finite weight's hard-error flag.
  quaids_workflow_test.e  # Milestone 21 seed: 32 checks (27 + 5 at
                    #   Milestone 26) -- parity checks proving
                    #   quaidsWorkflowFit() returns the same fit,
                    #   mean-point shares/elasticities, robust coefficient
                    #   SE, and robust propagated shares/elasticity SE as
                    #   explicit calls to the underlying public APIs; also
                    #   checks restriction/model-choice summaries,
                    #   Milestone 24's compact preflight summary,
                    #   quaidsWorkflowScenarioFit() classical/robust
                    #   welfare parity, and (Milestone 26)
                    #   quaidsWorkflowFit()'s own optional estimator
                    #   `weight` argument against direct quaidsFit()/
                    #   quaidsPreflight()/quaidsRobustFit() calls.
  quaids_survey_test.e  # Milestone 26: 19 checks -- exact-identity
                    #   regression guard (weight omitted / an explicit
                    #   uniform weight reproduce the pre-Milestone-26
                    #   unweighted behavior byte-for-byte, across
                    #   quaidsFit()/quaidsPreflight()/quaidsRobustFit()/
                    #   quaidsRobustBootstrapFit()); invalid-weight
                    #   preflight error; the sqrt(weight)-bread vs.
                    #   plain-weight-meat sandwich convention, checked
                    #   directly against a deliberately-wrong alternative
                    #   (this milestone's single easiest detail to get
                    #   backwards); and the core non-vacuous check --
                    #   weighted recovery beats naive on
                    #   _quaidsSurveyWeightedDGP's informatively-sampled
                    #   fixture.
  package_public_api.e   # Milestone 7: installed-package release gate --
                    #   `library quaids;` (not #include) against a real
                    #   install, exercising quaidsControlCreate/
                    #   getDefaultQuaidsControl, quaidsFit/printQuaids/
                    #   quaids, quaidsFull, quaidsElasFit/quaidsElas/
                    #   printQuaidsElas, quaidsSlutzky,
                    #   quaidsPreflight/printQuaidsPreflight,
                    #   quaidsWorkflowFit preflight summary fields,
                    #   quaidsHomogeneityTest/quaidsJointTest. Builds its
                    #   own small inline synthetic dataset rather than
                    #   reusing quaidsfixtures.src (tests/-only, not part
                    #   of the installed package). Run after
                    #   scripts/run_release_verification.ps1
                    #   -InstallArtifact. See "Milestone 7: package build
                    #   and release tooling" below. Milestone 9: found and
                    #   fixed a real gap here -- printQuaids() and
                    #   quaidsElas() were never actually called (only
                    #   their split components were), so a load-order bug
                    #   or stale .lcg entry specific to either one could
                    #   have passed this gate undetected. See "Milestone 9"
                    #   below.
  run_source_tests.ps1    # Milestone 7: runs verify_package_manifest.ps1
                    #   then all tgauss test files above (except
                    #   quaids_convergence_sweep.e, deliberately not a
                    #   pass/fail gate -- see "Testing status" below),
                    #   checking this repo's own PASS/FAIL-line convention
                    #   (not just tgauss's exit code). Milestone 15 adds
                    #   -SkipBootstrap (default: NOT skipped when this
                    #   script is run directly, matching the existing
                    #   -SkipPubtable/-SkipCurvature naming pattern) --
                    #   but .github/workflows/tests.yml's automatic
                    #   push-triggered CI run passes -SkipBootstrap
                    #   explicitly, since quaids_curvature_bootstrap_test.e
                    #   adds ~45-50s, close to doubling that run's
                    #   baseline. Run without the flag locally or during
                    #   release verification to exercise it.
  run_convergence_sweep.ps1  # Milestone 12: runs quaids_convergence_sweep.e
                    #   and captures its output to
                    #   tests/convergence_sweep_report.txt (gitignored --
                    #   regenerate on demand). Not run by
                    #   run_source_tests.ps1.
  verify_package_manifest.ps1  # Milestone 7: package.json src array vs.
                    #   actual src/ directory consistency (no dupes,
                    #   nothing missing, nothing unlisted except the
                    #   documented pubtable_quaids.src allowlist entry).
scripts/
  build_lcg.ps1     # Milestone 7: writes lib/<name>.lcg, the plain-text
                    #   symbol-location catalog GAUSS's `library`
                    #   mechanism reads (verified against the real,
                    #   installed qardl.lcg/pubtable.lcg format -- not a
                    #   stub). Run against an INSTALLED copy of the
                    #   package (e.g. the staging dir build_package.ps1 +
                    #   Expand-Archive produce), not this repo directly.
  build_package.ps1 # Milestone 7: stages package.json plus whichever of
                    #   README.md/CHANGELOG.md/CITATION.cff/CITATION.md/
                    #   LICENSE/llms.txt exist, plus src/docs/examples/
                    #   scripts/tests, strips generated run artifacts, and
                    #   zips "<name> <version>.zip" in the repo root (already
                    #   gitignored via the pre-existing `quaids *.zip`
                    #   pattern from Milestone 0).
  verify_release_artifact.ps1  # Milestone 7: checks a built .zip's name/
                    #   CHANGELOG.md entry match package.json's version, no
                    #   stale artifacts sit in the repo root, and the
                    #   archive contains every file package.json's src
                    #   array promises plus this repo's current root
                    #   files/dirs -- README.md and all four docs/*.md
                    #   files added to the required-entries list at
                    #   Milestone 8, once they existed.
  run_release_verification.ps1 # Milestone 7: orchestrator -- source tests,
                    #   build/verify the release artifact, optionally
                    #   install it into a GAUSS package directory
                    #   (defaults to <GaussHome>/pkgs, i.e. for real,
                    #   alongside every other package on the machine, only
                    #   with -InstallArtifact), then the installed-package
                    #   public API gate.
docs/
  COMMAND_REFERENCE.md  # Milestone 8: index of every public proc, grouped
                    #   (control struct, estimation, hypothesis tests,
                    #   elasticities/diagnostics, optional pubtable
                    #   reporting), linking to one command-reference/*.md
                    #   page each. Cross-checked against the actual source
                    #   by tests/verify_package_manifest.ps1 (every
                    #   documented proc must exist in src/; every linked
                    #   page must exist).
  USAGE_GUIDE.md    # Milestone 8: choosing an API (quaidsFit vs. quaids
                    #   vs. quaidsFull), the LA-AIDS/iterated-AIDS/QUAIDS
                    #   switch table, why IV is always required, symmetry/
                    #   homogeneity/overidentification workflow,
                    #   elasticities at arbitrary points, pubtable
                    #   reporting, current limitations.
  METHODOLOGY_NOTES.md  # Milestone 8: the estimator itself -- Stone vs.
                    #   translog price index, the QUAIDS quadratic term,
                    #   the IV control-function approach (and how it
                    #   differs from R micEconAids's 3SLS), the full
                    #   estimation algorithm phase-by-phase, elasticity/
                    #   Slutzky formulas, citing Deaton & Muellbauer (1980)
                    #   and Banks, Blundell & Lewbel (1997).
  FEATURE_SUPPORT_MATRIX.md  # Milestone 8: LA-AIDS x iterated-AIDS x
                    #   QUAIDS support for IV, hypothesis tests,
                    #   elasticities, Slutzky, curvature, zero correction,
                    #   robust/cluster inference, dataframe API,
                    #   pubtable export, synthetic/published validation.
  command-reference/  # Milestone 8: one *.md page per public proc (18
                    #   pages) -- Purpose/Format/Parameters/Returns/
                    #   Remarks/Examples/Source/See Also, matching
                    #   gauss-qardl's page template. Covers every proc in
                    #   quaids.sdf's load-bearing src/ files plus the
                    #   optional pubtable_quaids.src adapter's 6 procs
                    #   (documented despite being outside package.json's
                    #   src array, since they're real public API surface).
.github/
  workflows/
    tests.yml       # Milestone 14: CI on a self-hosted GitHub Actions
                  #   runner (this repo's tests need licensed GAUSS,
                  #   unavailable on GitHub-hosted runners). Triggers on
                  #   push to master only (plus workflow_dispatch), never
                  #   pull_request -- a deliberate public-repo security
                  #   choice, since self-hosted runners on public repos
                  #   are a real fork/PR code-execution risk otherwise.
                  #   Runs tests/run_source_tests.ps1 -SkipBootstrap (the
                  #   -SkipBootstrap flag added at Milestone 15 to keep
                  #   the routine per-push run fast -- see "Milestone 15"
                  #   below). See "Milestone 14: continuous integration"
                  #   below for the shell-invocation subtlety this
                  #   required (shell: cmd, not the default
                  #   shell: powershell).
package.json      # GAUSS package manifest (name: quaids, version: 0.21.0,
                  #   license: MIT). pubtable_quaids.src deliberately not
                  #   listed in its src array -- see "Milestone 6" below.
                  #   quaidscurvature.src IS listed (required public API),
                  #   so deps now lists "optmt" -- this library's first
                  #   real external package dependency (Milestone 10).
                  #   quaidswelfare.src (Milestone 11) is also listed but
                  #   adds no new dependency. Milestone 15 (bootstrap
                  #   standard errors) adds new public API inside
                  #   quaidscurvature.src but no new src array entry or
                  #   new dependency. quaidsworkflow.src (Milestone 21
                  #   seed) is listed last because it composes existing
                  #   public fit/post-estimation APIs. Milestone 26
                  #   (sampling-weighted estimation) likewise adds no new
                  #   src array entry or dependency -- every change is an
                  #   in-place edit to six already-listed files.
LICENSE           # MIT, copyright Eric Clower.
CITATION.cff      # Citation metadata; cites Deaton & Muellbauer (1980) and
                  #   Banks, Blundell & Lewbel (1997).
CHANGELOG.md      # Keep a Changelog style version history, reconstructed
                  #   from this file's own milestone records and maintained
                  #   through the current package version. This repo now has
                  #   an automated commit/push process (discovered at
                  #   Milestone 11 via `git
                  #   log`, not set up by any tool call in this
                  #   conversation) plus explicit commits at milestone
                  #   breakpoints per the repo owner's request -- "nothing
                  #   in this repo has been committed" is no longer
                  #   accurate as of Milestone 11 (see "Milestone 11:
                  #   welfare measures" below and the repo's own commit
                  #   history for current status).
README.md         # Milestone 8: front door -- install (Tools > Install
                  #   Application, or the scripts/ build+install
                  #   one-liner), quick start, model-choice summary,
                  #   feature list, links into docs/, testing/release
                  #   commands.
.gitignore        # Compiled .gcg artifacts, tmp/, .claude/, packaged zips,
                  #   generated `output file=...` run artifacts.
GOLD_STANDARD_TODO.md  # Living roadmap: release blockers, milestones,
                  #   definition of done. Read this before any nontrivial
                  #   change and update it as milestones close.
```

The original ten-milestone roadmap is complete, plus Milestones 11-26:
0 (repo hygiene), 1 (API/output-schema baseline), 2 (modular source split +
dataframe entry point), 3 (validation fixtures, including published-data
cross-implementation validation), 4 (hypothesis testing completeness), 5
(elasticities/diagnostics generalization), 6 (reporting via `pubtable`),
7 (package build and release tooling), 8 (documentation), 9 (final gold
standard integration gate), 10 (curvature imposition via Diewert-Wales,
requested by the repo owner after Milestone 9 closed), 11 (exact
welfare measures, requested by the repo owner after Milestone 10 closed),
12 (numerical reliability of the iterated estimator, requested by the
repo owner after Milestone 11 closed), 13 (QUAIDS curvature
imposition -- extending Milestone 10's AIDS-only support, requested by
the repo owner after Milestone 12 closed), 14 (continuous integration via
a self-hosted GitHub Actions runner, requested by the repo owner after
Milestone 13 closed), 15 (bootstrap standard errors for
`quaidsCurvatureFit()`, requested by the repo owner after Milestone 14
closed), 16 (predicted budget shares at an arbitrary point,
`quaidsSharesFit`, the first item of a full-demand-system-workflow
outline the repo owner asked for after Milestone 15 closed, being worked
through in order), 17 (a standalone AIDS-vs-QUAIDS specification
test, `quaidsQuadraticTest` -- a Wald test of whether the quadratic
log-expenditure term is needed, the second item of that same outline),
18 (percentile bootstrap confidence intervals,
`quaidsCurvatureBootstrapCI`, the third item of that outline -- building
it found and fixed a real, silent bug present since Milestone 10/15: a
row-major-vs-column-major `reshape()` mismatch had scrambled
`quaidsCurvatureFit()`'s `se` and `quaidsCurvatureBootstrapFit()`'s
`seBoot`, invisible to permutation-invariant shape/sign/finiteness
checks), 19 (zero-budget-share correction via Shonkwiler-Yen,
`quaidsZeroFit`, the fourth and largest item of that outline -- required
a real reformulation of the textbook method to fit onto this codebase's
shared-design-matrix Kronecker-product estimation core, and is
deliberately unconstrained-only in this first pass), 20 (robust and
cluster-robust standard errors, `quaidsRobustFit`/
`quaidsRobustBootstrapFit`, the fifth and final item of that outline --
genuinely new math generalizing every other covariance in this library's
pooled, homoskedastic sandwich to a per-observation or per-cluster score
aggregation, unified through one `clusterId` argument), 21 (the applied
workflow driver, `quaidsWorkflowFit`/`quaidsWorkflowScenarioFit`,
bundling the fit and its most common post-estimation outputs into one
call -- opening a second, post-outline phase focused on applied workflow
support), 22 (robust inference propagation, `quaidsRobustCovariance`/
`quaidsRobustBootstrapCovariance`, closing Milestone 20's own "does not
propagate automatically" gap on request), 23 (preflight data/design
diagnostics, `quaidsPreflight`, a non-gating pre-estimation check layer),
24 (a compact preflight summary echoed into the workflow output), 25
(an opt-in sampling-weighted workflow evaluation point,
`quaidsSurveyWorkflowFit`, explicitly scoped short of full survey-design
estimation), and 26 (sampling-weighted estimation itself,
`quaidsFit`'s optional `weight` argument with a matching weighted/
clustered sandwich SE via `quaidsRobustFit`/`quaidsRobustBootstrapFit` --
closing Milestone 25's own biggest deferred item, and this project's
biggest departure yet from the "new sibling proc, don't touch shipped
code" convention, since weighting touches essentially the whole
estimation core -- see each milestone's own section below for the real
bugs found and fixed along the way, including several in already-shipped
Milestone 20 code).

**The package is now actually installed** at `c:\gauss26\pkgs\quaids`
(Milestone 7), alongside `qardl` and `pubtable` on this machine --
`library quaids;` works. Rebuild/reinstall after any `src/` change with
`powershell scripts\run_release_verification.ps1 -BuildArtifact
-ForceArtifact -InstallArtifact` (run from the repo root).

**pubtable is installed in this environment** at `c:\gauss26\pkgs\pubtable`
(package version `1.0.0`) and is what `src/pubtable_quaids.src` (Milestone
6) targets. Like R/Python below, it is not a repo dependency — nothing in
`src/` or the core test suite requires it; only
`examples/pubtable_export_example.e` and `tests/quaids_pubtable_test.e` do.

**R and Python are installed in this environment** (as of 2026-07-20, with
explicit repo-owner approval) for cross-implementation validation: R 4.5.0
at `C:\Program Files\R\R-4.5.0\bin\` (not on `PATH` — invoke `Rscript.exe`
by full path) with the CRAN package `micEconAids` (and its dependencies
`micEcon`, `systemfit`) installed to the user library; Python 3.12 with
`numpy`/`pandas`/`scipy` installed via `pip`. Neither is a repo dependency
(nothing in `src/` or the test suite requires them to run) — they exist
only to regenerate/extend the published-validation reference numbers in
`tests/fixtures/published/`.

**Milestone 2 scoping note**: `quaidsFit()`'s starting-value construction,
iteration loop, variance computation, overidentification test, and symmetry
test/symmetry-constrained stage were deliberately *not* further split into
separate procs/files. They share heavily mutated intermediate state across
phases (`m`, `gg`, `gw`, `ng`, and iteration-loop-local `_beta`/`lambda`/
`lx`/`b_p`/`lx2` all flow into the variance and overidentification-test
formulas) — splitting them now would mean either passing a long, brittle
parameter list between new procs or bundling "in-progress fit state" into a
new struct, for a maintainability win that's real but not urgent on
pre-alpha code with no other consumers yet. Revisit after Milestone 3's
validation/fixture harness exists to safety-net a deeper refactor of this
kind — see `GOLD_STANDARD_TODO.md`'s "Roadmap Rules" for why that ordering
matters. The IV first stage (`quaidsiv.src`) was extracted because it
*is* cleanly separable: it only reads `intcpt`/`prices`/`instr`/`totexp`/`n`
and its outputs (`prices` converted to relative, `u`, `zzi`, `m1`, and the
IV diagnostics) are either consumed once downstream or not at all after this
phase.

## The `quaidsFit()` / `printQuaids()` / `quaids()` split

**`quaidsFit()` is the primary, silent, struct-returning entry point.** It
does 100% of the estimation with zero printing and returns a `quaidsOut`
struct (defined in `src/quaids.sdf`, ~75 fields, grouped by phase: metadata,
IV first-stage diagnostics, homogeneity-constrained stage, overidentification
test, symmetry test + symmetry-constrained stage, and the final
absolute-price-form `b`/`v`/`bS`/`vS`).

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
```

**`printQuaids(qOut)`** reproduces the estimation-stage console report (IV
first-stage table, iteration summary, homogeneity-constrained coefficient
table, overidentification test, symmetry test + symmetry-constrained table)
from a `quaidsOut` struct alone. It does *not* print elasticities,
descriptive statistics, or the Slutzky diagnostic — those are separate,
explicitly-callable reports (`quaidsElas()`, `quaidsSlutzky()`), since
Milestone 5 plans to generalize elasticities to arbitrary evaluation points
rather than bake a fixed set into the struct.

**`quaids()` is the original, backward-compatible call** — unchanged
signature and unchanged printed behavior:

```gauss
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, struct quaidsControl aCtl);
```

It calls `quaidsFit()`, calls `printQuaids(qOut)`, then reproduces the
legacy elasticities-at-four-points / descriptive-statistics / Slutzky report
exactly as the pre-Milestone-1 `aids()` proc did, then returns
`(qOut.b, qOut.v, qOut.bS, qOut.vS)`. **Verified byte-for-byte at the time**:
running `quaids()` post-refactor against the same fixed-seed synthetic
dataset used at Milestone 0 produced output identical to the pre-Milestone-1
code, including the per-iteration convergence log (reproduced from a stored
`qOut.iterHistory` matrix rather than printed live during iteration, since
`quaidsFit()` cannot print, but reproduced with matching iteration
numbers/error values so the printed table was unchanged). **This byte-parity
claim is historical, not current**: Milestone 3 found and fixed a real bug
in the Stone-index starting-value computation (see "Milestone 3: real bug
found and fixed" below), which intentionally changes numerical output for
`aCtl.maxiter==1` calls, and shifts the iteration path for `aCtl.maxiter>1`
calls, relative to this old baseline. The *refactor* (Milestone 1's
estimation/printing split) was and remains behavior-preserving with respect
to the code as it existed at that time; a later, separate, intentional
correctness fix is not a violation of that.

- `w` — `TxN` budget shares.
- `intcpt` — `TxK` extra intercept-shifter variables (demographics etc.), or
  `0` for none.
- `prices` — `TxN` **log** prices (absolute, not relative — the proc converts
  internally).
- `totexp` — `Tx1` log total expenditure (treated as endogenous).
- `instr` — `TxH` instruments for log total expenditure.
- `aCtl` — `quaidsControl` struct (see below).

`quaids()`'s (`b1`/`v1`/`b2`/`v2`) return semantics are unchanged: if
`aCtl.homogenous == 1`, `b1`/`v1` are homogeneity-constrained
estimates/covariance and `b2`/`v2` are homogeneity+symmetry-constrained. If
`aCtl.homogenous == 0`, `b1`/`v1` are the unconstrained (reparametrized)
estimates and `b2`/`v2` are `0`. Same values are available struct-side as
`qOut.b`/`qOut.v`/`qOut.bS`/`qOut.vS`; `qOut.bestB`/`qOut.bestV` always hold
"whichever is the most-constrained estimate actually fit" (symmetric if
homogeneity was imposed, else the recovered unconstrained fit) — this is
what elasticities/Slutzky are evaluated against, matching what the original
code's reused `b_s`/`v_s` locals held at that point regardless of branch.

**A pre-existing anomaly, preserved not fixed**: in the symmetry-constrained
table's "Residuals of instrumental regressions" row, the original code pairs
point estimates from the *homogeneity-stage* `b` with standard errors/t/p
from the *symmetry-stage* fit (both stages leave the IV-residual coefficient
block numerically identical, but the printed SE/t/p come from the
symmetry-constrained covariance). This looks like it could be an original
authoring inconsistency. It was carried over unchanged (`qOut.homogB` paired
with `qOut.symcSE`/`qOut.symcT`/`qOut.symcPvt` in `printQuaids()`) since
Milestone 1 is a structure-preserving refactor, not a correctness pass —
flag for review during Milestone 3/4 validation work.

## `quaidsFull()` — dataframe entry point

```gauss
struct quaidsOut qOut;
qOut = quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars, aCtl);
```

`data` is an already-loaded dataframe (`data = loadd(fname);` or built with
`asDF()`/`dfaddcol()` — `quaidsFull()` does not load files itself).
`shareVars`/`priceVars` are string arrays of column names, **matched by
position** (`shareVars[i]` and `priceVars[i]` must be the same good — there
is no name-matching magic). `totexpVar` is a single column-name string.
`instrVars` is a string array. `extraVars` is a string array, or the scalar
`0` for "no extra intercept shifters" (matching `quaidsFit()`'s `intcpt ==
0` convention — checked via `type(extraVars)`, not value equality, since
comparing a string array to `0` with `==` is itself a type-mismatch trap).
Internally it just selects the named columns and calls `quaidsFit()` — no
duplicated estimation logic. Verified numerically identical to the
equivalent `quaidsFit(matrices...)` call, including the `extraVars == 0`
path, by `tests/quaids_formula_parity_test.e` (17 checks).

### GAUSS type-system notes learned while building this

- `type()` codes actually seen on GAUSS 26.1.1: `6` = plain matrix (and, it
  turns out, **also** a dataframe and a dataframe column selection —
  `data[., "col"]` and `data[., stringArray]` both return plain type-`6`
  matrices, not some distinct "dataframe" type; dataframe-ness doesn't leak
  into downstream numeric code). `13` = native `string`. `15` = native
  `string array`. Legacy `0$+"X"$+ftocv(...)`-built character matrices are
  *also* type `6` (matrix) — this is why `quaidsOut`'s name-vector fields
  are declared `matrix`, not `string array` (see below).
- Struct-field assignment is strict about matching one of these types
  exactly; ordinary expressions are more forgiving. `"CONSTANT"|xnam`
  (vertical-concatenating a native string literal onto a legacy char-matrix)
  silently coerces and works fine as an *expression*. Directly assigning a
  bare native string literal to a `matrix`-typed struct field
  (`qOut.xnam = "CONSTANT";`) does **not** coerce — it throws `G0071 Type
  mismatch`. Use the `0$+` prefix idiom (`0$+"CONSTANT"`) to get the
  legacy-char-matrix form explicitly before assigning into a `matrix`
  struct field.

### Two pre-existing bugs found (and fixed) by the formula parity test

`tests/quaids_formula_parity_test.e` exercises `intcpt == 0` (no extra
intercept-shifter columns) — a code path that existed unchanged from the
original `aids_rev.src` but had **never actually been run** by anything in
this repo before, because the synthetic example always passes a non-zero
`intcpt`. Writing a real caller for that path immediately surfaced two bugs,
both in `quaidsFit()`'s name-setup block (`src/quaids.src`, just after the
`intcpt == 0` branch):

1. `xnam` was never assigned inside the `if intcpt == 0;` branch (only in
   `else;`), so the following unconditional `xnam = "CONSTANT"|xnam;` read
   an uninitialized variable — `G0152 Variable not initialized`.
2. After fixing that by assigning `xnam = "CONSTANT";` in the `if` branch,
   assigning it into `qOut.xnam` (`matrix` type) threw `G0071 Type
   mismatch`, per the type-system note above — a native string literal
   isn't the same type as the char-matrix `xnam` holds in the `else`
   branch. Fixed with `xnam = 0$+"CONSTANT";`.

Both fixes are scoped to the previously-dead branch only; the `else` branch
(exercised by every prior test) was not touched, and re-running the
Milestone 0/1 parity and schema tests after the fix confirmed zero change to
already-verified behavior. Lesson: a code path with no caller is a code path
with no evidence it works, regardless of how long it's sat there unchanged.

## Milestone 3: synthetic validation findings

`tests/quaids_synthetic_validation_test.e` fits `quaidsFit()` against
`tests/quaidsfixtures.src`'s `_quaidsSyntheticDGP()` (a 5-good dataset with
homogeneity/adding-up true by construction, parameterized by whether the
true model has a quadratic term and whether total expenditure is genuinely
endogenous) and checks recovered parameters against the true ones, not just
that the code ran. Full tolerance rationale is in `GOLD_STANDARD_TODO.md`'s
Milestone 3 "Tolerance Policy" section — summary: structural coefficients
within `0.10`, the IV-residual coefficient row within `0.50` (it's
consistently the noisiest row by roughly 10x — a control-function term on
an estimated regressor, not a "deep" structural parameter), LA-AIDS within
`1.20` throughout (Stone-index approximation bias is a real property of
that method).

### Seed sensitivity — a real numerical-reliability finding

Calibrating those tolerances required an 8-seed probe (`tobs=3000`,
`aCtl.err=.0001`, `aCtl.maxiter=100`) that turned up something worth
flagging prominently: **the iterative estimator (QUAIDS and iterated linear
AIDS) fails to converge, or "converges" to wildly wrong estimates (errors
of magnitude 200–2500 against true parameters of magnitude ~0.1–2), for
roughly half of random seeds** in this DGP family. The failure pattern
tracked with the *price* draw for a given seed, not with whether the model
had a quadratic term or genuine endogeneity — the same seeds failed (or
succeeded) whether `quadratic`/`endogenous` were on or off, holding
`prices` fixed. This points at the iteration's conditioning being sensitive
to the specific price data realization, not at a bug specific to QUAIDS or
IV handling.

The synthetic fixtures use `seed=204`, one of the seeds confirmed (by that
probe) to converge cleanly across all six model/endogeneity combinations —
documented as such in the test file's comments, not silently cherry-picked.
**Partially explained, not fully fixed**: the Stone-index starting-value
bug (below) likely explains *some* of this non-convergence (a materially
wrong starting point makes convergence failure more likely), and was fixed.
But re-running the same 8-seed probe after the fix would be needed to know
how much of the non-convergence rate that actually accounts for — not done,
since it wasn't necessary to validate the fix itself (the published-data
comparison already did that directly). `quaidsFit()`'s iteration still has
no globally-convergent guarantee, no damping, and no fallback for a bad
starting point. Worth a dedicated numerical-reliability pass (analogous to
`gauss-qardl`'s Milestone 13) once more of the roadmap exists to build on.

### Milestone 3: real bug found and fixed — Stone-index starting value

The published-data cross-check below initially disagreed with R by roughly
5x on `beta` — far beyond the few-percent gap expected between two
different-but-valid IV algorithms. Root cause, in `quaidsFit()`'s "STARTING
VALUE" block (`src/quaids.src`): `stone = prices*meanc(w)` was applied to
`prices` *after* it had already been converted to relative form in columns
`1:n-1` (each minus the reference good's price) while column `n` stayed
absolute. Weighting that mixed matrix by mean shares does not compute the
Stone price index — algebraically (verified by direct derivation, then
confirmed empirically by patching the formula in an isolated Python check
and watching the gap with R close) it computes
`StandardStoneIndex − ln(p_n)·(1 − meanShare_n)`, a distortion tracking the
reference good's own price trend.

**Impact**: `aCtl.maxiter==1` (LA-AIDS) never iterates past this starting
value, so the distorted deflator *was* the final answer for every LA-AIDS
call using default starting values (`aCtl.b0==0`) — not an occasional
glitch, every time. For `aCtl.maxiter>1`, this was only a bad *starting
point*; the iteration loop's own `a_p` formula is correct and unaffected.

**Fix** (`src/quaids.src`, "STARTING VALUE" block): reconstruct absolute
prices before computing `stone`:
```gauss
stone = (prices[., 1:n-1] + prices[., n])~prices[., n];
stone = stone*meanc(w);
```
rather than changing the relative-price convention used elsewhere in the
proc. Confirmed this fix (not something else) closed the gap: after
patching, GAUSS's estimates matched R's independent `micEconAids`
3SLS-with-instrument reference to within ~0.021 max absolute difference
(see below) instead of ~5x off.

**This changes numerical output** for `aCtl.maxiter==1` calls, and shifts
(not necessarily worsens) the iteration path for `aCtl.maxiter>1` calls,
relative to every prior milestone's frozen baseline — expected and correct,
since those baselines were captured from the original, buggy code.
`examples/quaids_example.e` (seed 11) was already one of the non-converging
seeds identified above even pre-fix; its output changed too and was not
re-tuned to a better-behaved seed as part of this fix.

### Published-data validation against R and Python

`tests/quaids_published_validation_test.e` fits `quaidsFit()`
(`aCtl.linear=1, aCtl.maxiter=1` — LA-AIDS) on `Blanciforti86` (annual U.S.
food-consumption data, 1947–1978, 4 food groups — see
`tests/fixtures/published/SOURCE.md`), instrumenting `log(xFood)` (total
food expenditure) with `log(xAgg)` (total aggregate expenditure, `corr ≈
0.97` in logs — a strong, genuinely informative instrument, not a
near-degenerate one: first-stage `R² ≈ 0.998`).

- **R** (`tests/fixtures/published/generate_r_reference.R`,
  `micEconAids::aidsEst(..., instNames=...)`, dispatches to 3SLS via
  `systemfit`): **max abs difference from GAUSS ≈ 0.021** across
  alpha/beta/gamma. This is the hard assertion target in the test
  (tolerance `0.05`, real headroom above the observed gap). The residual
  ~0.02 is attributable to R's 3SLS and GAUSS's control-function
  (residual-inclusion) approach being different, both valid, IV algorithms
  for the same model — not expected to match bit-for-bit.
- **Python** (`tests/fixtures/published/python_reference_check.py`,
  hand-coded from the Deaton-Muellbauer equations — no comparably
  established Python AIDS package exists the way R has `micEconAids`):
  broadly consistent (one equation matches GAUSS almost exactly) but with
  larger residual differences on another equation than either the R
  comparison or an earlier, simpler Python no-IV replica showed. Since R —
  an independently authored, widely-used implementation — agrees closely
  with GAUSS on exactly the coefficients where this Python script diverges
  most, that divergence is attributed to the from-scratch replica (most
  likely in how it forms the IV-residual-augmented GLS weighting), not to
  GAUSS. Kept in the repo for transparency and as a documented starting
  point; **not used as a pass/fail assertion source**, only R is.
- Both scripts are runnable standalone from `tests/fixtures/published/`
  (`Rscript generate_r_reference.R`, `python python_reference_check.py`)
  to regenerate or extend these numbers.

## Milestone 4: new hypothesis tests

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);  // aCtl.homogenous = 0
{ stat, pval, df } = quaidsHomogeneityTest(qOut);
{ statJ, pvalJ, dfJ } = quaidsJointTest(qOut);
```

Both (`src/quaidstests.src`) are Wald chi2 tests and **require an
unconstrained fit** (`qOut.homogenous == 0` — they error clearly if not).
They read `qOut.b`/`qOut.v` (the final, absolute-price-form unconstrained
gamma matrix and its covariance) and build a restriction vector/covariance
via a selection matrix `L`: `stat = (L'*vec(b))' * inv(L'*V*L) * (L'*vec(b))`.

- **`quaidsHomogeneityTest`**, `df = n-1`: tests `sum_j gamma_ij = 0` for
  each independently-estimated equation (equation `n` is recovered via
  adding-up and adds no information to a Wald test).
- **`quaidsJointTest`**, `df = (n-1) + (n-1)(n-2)/2`: homogeneity's
  restrictions plus symmetry (`gamma_ij = gamma_ji`, `i<j`, `i,j=1..n-1`).
  A symmetric gamma with adding-up already implies homogeneity (row `i`
  sum = column `i` sum by symmetry = 0 by adding-up) — there's no separate
  "symmetry given adding-up alone" test on an unconstrained fit; use this
  joint test, or the existing symmetry-given-homogeneity test if
  homogeneity itself isn't in question.

Both are validated for size *and* power in `tests/quaids_hypothesis_tests_test.e`
(19 checks) — full derivation, both dead ends hit while getting there, and
why size+power (not just "it runs") is the standard applied, are in
`GOLD_STANDARD_TODO.md`'s Milestone 4 section. Short version: the first
implementation read `qOut.b`'s row layout wrong (mistook an internal
pre-recovery reparametrization, described in `quaidsFit()`'s own
docstring, for the final post-recovery form) and rejected a true null with
`p≈0` — caught immediately by running the size check before trusting the
formula, not by code review.

Same milestone also **exercised the overidentification test for the first
time in this repo's history**: every fixture through Milestone 3 used
exactly-identified instruments (`ninst==nu`), so `qOut.overidValid` was
always `0` and that branch of `quaidsFit()` had never actually run. A new
2-instrument fixture in the same test file confirms it runs, has the
right shape (`overidGamma` is `ninst x n`) and df (`ninst-nu`), and
doesn't spuriously reject when both instruments are valid by construction.

**Why the existing symmetry/overID tests weren't extracted into their own
procs** (unlike the two new ones above): they depend on heavily mutated
intermediate iteration/variance state (`Ji`, `Sgma`, `S`, `O`, `D`, `zzi`,
`m1`, iteration-final `_beta`/`lambda`/`lx`/`b_p`/`lx2`) — the same coupling
Milestone 2 already declined to untangle for the estimation core itself,
unchanged here since nothing about the estimation core moved. The two new
tests avoid this entirely: they only need the *finished* `qOut.b`/`qOut.v`,
not any of that intermediate state.

## Milestone 5: elasticities generalization

```gauss
struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPoint, pricesPoint, totexpPoint, aCtl);
call printQuaidsElas(elasOut);   // optional -- omit for a silent call
```

**The roadmap's framing turned out to be slightly off, worth recording**:
it described `quaidsElas()`/`quaidsElas_()` as "hardcoded to mean/Q1/
median/Q3," but `quaidsElas_()` already took `intcpt`/`prices`/`totexp` as
*point arguments*, not sample statistics — it always could evaluate
anywhere. The actual gap was that `quaidsElas()` mixed computation with
printing (no way to get a result object silently) and that `quaids()`
(the legacy wrapper) only ever *called* it at four fixed points. Fixed by
splitting `quaidsElas()` the same way Milestone 1 split `quaidsFit()`:

- **`quaidsElasFit(b, v, intcpt, prices, totexp, aCtl)`** — silent, returns
  `struct quaidsElasOut` (`er`/`ep`/`epc` point estimates, `ser`/`sep`/`sepc`
  raw numeric delta-method standard errors — not the pre-formatted
  `"(0.123)"` display strings the old code produced inline, which stay in
  the printer where they belong).
- **`printQuaidsElas(elasOut)`** — the separated printer.
- **`quaidsElas(...)`** — unchanged signature, now a thin wrapper
  (`quaidsElasFit()` then `printQuaidsElas()`). **Verified byte-for-byte**
  identical printed output against the pre-split version on
  `examples/quaids_example.e`.
- `quaidsElas_()` itself is untouched — it was already the right shape.

**Validating "arbitrary point" actually works, not just "compiles"**:
`tests/quaids_elasticities_test.e` checks three **exact** algebraic
identities that any valid AIDS/QUAIDS elasticity set must satisfy given
adding-up/homogeneity (which this estimator always imposes by
construction) — these are consequences of the functional form, not
separately-estimated statistical quantities, so they should hold to
floating-point precision, not just "within tolerance":

- Engel aggregation: `sum_i(w_i * er_i) = 1`
- Cournot aggregation: `sum_i(w_i * ep_ij) + w_j = 0`, for each price `j`
- Elasticity homogeneity: `sum_j(ep_ij) + er_i = 0`, for each good `i`

Checked at a real out-of-sample observation and at a fully synthetic
counterfactual (a hypothetical 20% price increase on one good, evaluated
at the sample-mean point otherwise) — both held to ~1e-16. A test proc
bug surfaced immediately while writing this (`struct quaidsElasOut
elasOut` must be declared in the parameter list for field access to work,
same as any GAUSS struct-typed parameter — omitting it gives a `G0008
Syntax error` on the first `.field` access, not a silent failure).

**`quaidsSlutzky()` needed no change**: it already accepts an arbitrary
`intcpt`/`prices`/`totexp` sample (any number of rows) as input, so "keep
the Slutzky diagnostic general" was already true going into this
milestone.

**Curvature imposition (Diewert-Wales Cholesky reparametrization) — scoped,
not implemented**: the roadmap listed this as P2/"if justified." Even
`micEconAids` (the reference implementation used for Milestone 3's
cross-check) only offers curvature *diagnosis* (`aidsMono()`,
`aidsConcav()` — post-hoc checks, matching what `quaidsSlutzky()` already
does here), not *imposition* as a constrained-estimation mode. That's real
evidence against this library needing to leap ahead of the reference
implementation on a deliberately-optional item — revisit only if a
concrete use case shows up.

### `quaidsControl` fields (`src/quaids.sdf` / `quaidsControlCreate()`)

| Field | Default | Meaning |
|---|---|---|
| `linear` | `0` | `1` = LA-AIDS (linear); `0` = QUAIDS (quadratic log-expenditure term) |
| `maxiter` | `50` | `1` = one-step linearized AIDS with Stone price index; `>1` = iterate |
| `homogenous` | `1` | `1` = impose homogeneity (and test/report symmetry); `0` = unconstrained |
| `alpha0` | `0` | Fixed value of the translog price-index intercept `alpha_0` |
| `err` | `.0001` | Relative parameter-change convergence tolerance |
| `othnam` | `""` | Optional alternate variable names for printed output |
| `b0` | `0` | Optional user-supplied starting values; `0` = use built-in starting values. For `quaidsFit()`, a supplied matrix must match the reduced raw coefficient matrix shape used by the homogeneity stage (`qOut.homogB`). For `quaidsZeroFit()`, it must match the zero-corrected coefficient shape (`zOut.b`) |
| `relax` | `1` | Milestone 12: under-relaxation factor for the iterated fixed-point update, `(0,1]`; `1` = no damping |

The `stone`, `aids`, and `varname` fields flagged as dead at Milestone 0 were
removed from `quaidsControl` at Milestone 1 (also dropped from
`quaidsControlCreate()`), since grep confirmed no read anywhere in
`quaids.src` and there are no external consumers of this pre-alpha struct yet.
Also added: `getDefaultQuaidsControl()` (alias for `quaidsControlCreate()`,
naming-convention parity with `gauss-qardl`'s `getDefault...Control` procs),
and structure-inference return typing (`proc (struct quaidsControl) =
quaidsControlCreate();`) so callers no longer need to pre-declare
`struct quaidsControl aCtl;` before assignment.

## Milestone 6: reporting via pubtable

```gauss
library pubtable;
#include quaids.sdf
#include quaidsutil.src
#include quaidsiv.src
#include quaidselas.src
#include quaidsslutzky.src
#include quaids.src
#include pubtable_quaids.src

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);            // one comparison column per good
call ptExport(coefTbl, "results.tex");   // .tex/.md/.csv/.rtf/.html/.xlsx by extension

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);

struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);  // 3x1: income, uncompensated, compensated
call ptExport(elasTbls[1], "income_elasticities.md");
```

`src/pubtable_quaids.src` is an optional adapter onto the `pubtable`
package (`c:\gauss26\pkgs\pubtable`), following the same pattern as
`pubtable`'s own bundled `pubtable_qardl.src`
(`c:\gauss26\pkgs\pubtable\src\pubtable_qardl.src`): every proc is guarded
by `#ifDef QUAIDS_SDF_INCLUDED` with a clear `_library_missing_error()` stub
otherwise, using the `(...)`/`dynargsGet()` variadic-argument idiom for the
struct-typed inputs.

**Location deliberately diverges from the `pubtable_qardl.src` precedent**:
that file lives *inside the installed pubtable package itself*, not inside
`gauss-qardl`'s own repo — pubtable ships adapters for both Aptech's own
estimators (`cmlmt`/`maxlikmt`/`optmt`/`tsmt`) and third-party ones
(`qardl`) bundled together in its own `src/`. Matching that physically for
`pubtable_quaids.src` would mean writing into a shared, installed package
outside this repo's git history, affecting every other project on this
machine that loads `pubtable` — asked of the repo owner explicitly rather
than assumed (see `GOLD_STANDARD_TODO.md`'s Milestone 6 section for the
full reasoning). Decision: `src/pubtable_quaids.src` stays inside this
repo, self-contained and git-tracked, matching every other file in `src/`
(none of which `#include` their own dependencies — the caller does). It is
**not** listed in `package.json`'s `src` array, since (unlike every other
file there) its `proc (struct ptModel) = ...`/`proc (struct ptTable) = ...`
return-type annotations sit outside the `#ifDef` guard and need
`pubtable.sdf`'s struct types declared unconditionally — adding it to
`src` would make `pubtable` a hard compile-time dependency for the whole
package, contradicting `package.json`'s empty `deps` array.

**Public procs**:
- `ptModelFromQuaids(name, qOut, eqIdx)` / `ptFromQuaids(qOut)` — one
  equation's (good's) coefficient column from `qOut.bestB`/`qOut.bestV`
  (the same "most-constrained estimate actually fit" elasticities/Slutzky
  are evaluated against), or an `n`-good comparison table via
  `ptModelCompare` (mirrors `pubtable_qardl.src`'s per-quantile
  `ptFromQardl`). Row order: intercept block (`qOut.xnam`) | price/gamma
  block (`GAMMA_`+good name, one row per good) | `BETA_LX` |
  `LAMBDA_LX2` (QUAIDS only) | one row per IV control-function residual
  term (`qOut.unam`, `qOut.nu` rows — always present, since `quaidsFit()`
  always treats log total expenditure as endogenous).
- `ptModelFromQuaidsElas(name, elasOut)` / `ptFromQuaidsElas(elasOut)` —
  income elasticities (`elasOut.er`/`elasOut.ser`) as a single-column
  `ptModel`/`ptTable`.
- `ptTablesFromQuaidsElas(elasOut)` — 3x1 `struct ptTable` array (income
  elasticities, uncompensated/Marshallian price elasticities, compensated/
  Hicksian price elasticities), suitable for `ptExportAll()`. The two
  price-elasticity tables are `n x n` matrices with a value row and an
  `(se)` row per good — built directly via the internal
  `_ptQuaidsElasMatrixTable()` helper rather than through `ptModel`, since
  that shape doesn't fit `ptModel`'s single-coefficient-vector design.
- `ptFromQuaidsFamily(x)` — `isStructType`-based dispatcher (`quaidsOut` ->
  `ptFromQuaids`, `quaidsElasOut` -> `ptFromQuaidsElas`), mirroring
  `pubtable_qardl.src`'s `ptFromArdlFamily`.

**A real type-system finding, empirically verified (not assumed from the
Milestone 1 lesson)**: `quaidsOut`/`quaidsElasOut` name vectors (`xnam`,
`wnam`, `unam`, ...) are legacy character matrices (GAUSS type 6); `pubtable`
struct fields for names/titles are natively typed `string array` (type 15)
or scalar `string` (type 13). Direct assignment throws `G0071 Type
mismatch` in *this* direction too (matrix into string-typed field), not
just the matrix-field direction CLAUDE.md already documented. Tested
several candidate conversions directly with `tgauss` rather than guessing:
`strtrim()` looked like the obvious builtin and instead errors ("Invalid
argument type") on a legacy char-matrix input — **there is no dedicated
conversion builtin**. The idiom that does work: concatenating the char-matrix
with a native string via `$|` forces element-wise conversion to a native
string array (`cm $| ""`, slicing off the trailing blank row for a full
array; a *scalar*-indexed single row of the result, `sa[1]` and not
`sa[1:1]`, further demotes to a true scalar `string`). Encapsulated as
`_ptQuaidsToStrArray()`/`_ptQuaidsToStr()` in `pubtable_quaids.src`.

**A real row-count bug, found by running against a real fit**: the first
draft of `ptModelFromQuaids()` built row names only from the blocks
`quaidsElas_()` reads (`qOut.xnam` | `GAMMA_`+`qOut.wnam` | `BETA_LX` |
`LAMBDA_LX2`) — but `quaidsElas_()` only *reads* the first
`1+nint+n+nendog` rows of `b` and silently ignores the rest; `qOut.bestB`
itself always has `qOut.nu` additional trailing rows for the IV
control-function residual coefficient(s). Running the adapter against a
real `quaidsFit()` QUAIDS/IV result immediately threw
`ptModelSetNames: termNames must contain 10 labels` (9 built vs. 10
actual) — caught by execution, not by re-deriving the formula on paper.
Fixed by appending `qOut.unam`-derived names for the trailing block, which
also makes the table more complete: `printQuaids()`'s own console report
already shows these as the "Residuals of instrumental regressions" row, so
a pubtable coefficient report should too. Re-verified against LA-AIDS
(`linear=1`) and QUAIDS (`linear=0`), both `homogenous=1` and
`homogenous=0`, to confirm the row-count formula holds across all four
combinations.

## Milestone 7: package build and release tooling

```powershell
# From the repo root, after any src/ change:
powershell scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact -InstallArtifact

# Or step by step:
powershell tests\run_source_tests.ps1                       # manifest check + all 7 tgauss tests
powershell scripts\build_package.ps1 -Force                 # writes "quaids <version>.zip", self-verifies
powershell scripts\run_release_verification.ps1 -InstallArtifact -SkipInstalledPackageTest  # installs to <GaussHome>/pkgs/quaids
tgauss -b -x package_public_api.e                            # from tests/, exercises `library quaids;`
```

Adapted from `gauss-qardl`'s `scripts/`/`tests/` release tooling, scaled
down for this repo's current scope — see `GOLD_STANDARD_TODO.md`'s
Milestone 7 section for the full list of what was intentionally omitted
(separate benchmark/examples-smoke scripts; this repo's equivalent
validation already lives in `run_source_tests.ps1`'s 7 tgauss files) and
what's deferred to Milestone 8 (`verify_release_artifact.ps1`'s
`requiredEntries` list doesn't require `README.md`/
`docs/COMMAND_REFERENCE.md` yet, since neither exists).

**The package is now really installed** at `c:\gauss26\pkgs\quaids`,
verified via `library quaids;` — not just staged/zipped. This was an
explicit repo-owner decision (asked via the same "writing outside this git
repo into shared machine state" reasoning as the Milestone 6 adapter-
location question), since fully validating the release pipeline requires
`library quaids;` to actually resolve, which only happens against GAUSS's
real, shared package directory.

**Three real bugs found by actually building, installing, and running the
package — not by re-reading the scripts**, the same "never trust a
derived formula/script without running it" standard this whole project has
followed since Milestone 3:

1. **`build_package.ps1`'s cleanup step deleted the entire staged
   `examples/`/`tests/` directories.** `Get-ChildItem -LiteralPath <dir>
   -Include <patterns> -Recurse` silently ignores `-Include` when the base
   path isn't itself a wildcard (a known PowerShell footgun) — `-Recurse`
   returned every file with no effective filter, and `Remove-Item -Force`
   deleted all of them. `verify_release_artifact.ps1` caught it
   immediately: "missing required entry: examples/quaids_example.e" on the
   first real run. Fixed by removing named generated files by explicit
   literal path, plus a `-Filter` (not `-Include`) pass for `*.log`, which
   doesn't have this bug.
2. **`build_lcg.ps1`'s proc-detection regex only matched one of the three
   GAUSS proc-declaration forms this codebase actually uses.** It matched
   `proc (struct X) = name(...)`, but not the bare-digit form `proc N =
   name(...)` (e.g. `proc 3 = quaidsElas_(...)`) or the no-return-spec
   form `proc name(...)` (e.g. `proc quantile(x, s);`). This silently
   dropped `quaidsSlutzky`, `_quaidsIVFirstStage`, `quaids()`,
   `printQuaids`, `quaidsElas_`, `printQuaidsElas`, `quaidsElas`, and a
   private `quantile` helper from the generated `.lcg` catalog — invisible
   to every `#include`-based test (which never goes through the catalog),
   surfaced only as `Undefined symbol` errors from `library quaids;` when
   `tests/package_public_api.e` called them for real. Fixed the regex to
   match all three forms; verified the regenerated catalog lists every
   proc in every `src/` file.
3. **A genuinely pre-existing bug in `src/quaids.src`, dating from before
   Milestone 0**: a private `quantile(x, s)` helper duplicating GAUSS's
   builtin `quantile()`. The original author clearly intended to delete
   it — it sat inside a comment block — but GAUSS comments do not nest,
   and the wrapping comment's own doc-header used an inner
   `/**...**/`-style block whose closing marker closed the *outer*
   "delete this" comment early, leaving the entire proc live and
   uncommented. Silently duplicated the builtin's behavior at its 3 call
   sites (`quaids()`'s legacy elasticities-at-four-points block) for as
   long as this repo has existed, invisible via `#include`-based
   compilation (which just locally shadows the builtin name) — but a real
   GAUSS builtin name cannot be redefined via `library`-based lazy
   loading, so it surfaced as `Undefined symbol: 'quantile'` resolving
   `quaids.src`'s *own* proc definition, once bug #2's fix let `library
   quaids;` actually try to load it. Fixed by deleting the dead code
   outright (the original author's evident intent), not just repairing
   the comment nesting — the 3 call sites now resolve to the GAUSS
   builtin, same as `quaidsslutzky.src`'s identical calls always have.
   Full 7-test regression suite re-run afterward; no assertion anywhere
   depended on the old local implementation's exact output (only the
   legacy, unasserted quartile-elasticities console block uses it), so no
   test needed updating.

Also hit and fixed a self-inflicted issue while writing bug #3's own
explanatory comment: GAUSS's lexer does not tolerate literal `"`
characters inside a `/* ... */` block comment — an odd count throws
`error G0097 String not closed`, even though the characters are inside a
comment, not a string. Caught immediately by re-running the test suite;
fixed by rephrasing without quote characters. Worth remembering when
writing GAUSS comments generally, not just here.

**No version bump for this milestone**: build/release tooling doesn't
change GAUSS public API surface, so per this repo's established policy
(version bumps track public API changes, not every milestone — see
Milestone 3's R/Python scripts) the package stays at `0.5.0`. Not tagging
a git release either: nothing in this repo has been committed yet as of
Milestone 7 (every milestone's changes are staged but uncommitted, per the
"never commit unless asked" policy), and a tag needs a commit to point at
— `CHANGELOG.md`/`package.json` versioning infrastructure is in place for
whenever the repo owner chooses to commit and tag.

## Milestone 8: documentation

Full doc set added: `README.md`, `docs/COMMAND_REFERENCE.md` plus 18
`docs/command-reference/*.md` pages (one per public proc, Purpose/Format/
Parameters/Returns/Remarks/Examples/Source/See Also template, matching
`gauss-qardl`'s page style), `docs/USAGE_GUIDE.md`, `docs/METHODOLOGY_NOTES.md`,
`docs/FEATURE_SUPPORT_MATRIX.md`. `CLAUDE.md` itself was already complete
and kept synchronized with every milestone since Milestone 0 — nothing
further needed there.

**Coverage includes the optional `pubtable` adapter**: all 12 procs across
`quaids.sdf`'s load-bearing `src/` files, plus all 6 procs in
`src/pubtable_quaids.src` (not in `package.json`'s `src` array, but real,
working, documented public API surface — the same reasoning `README.md`
and `USAGE_GUIDE.md` already apply to the optional `pubtable` export
path).

**Followed through on two Milestone-7-deferred TODOs rather than leaving
them stale**, both explicitly flagged in their own header comments at the
time:

- `tests/verify_package_manifest.ps1` now cross-checks
  `docs/COMMAND_REFERENCE.md` against the actual source the same way
  `gauss-qardl`'s version does: every documented proc must be defined
  somewhere in `src/` (using the same intentionally-unlisted allowlist the
  `package.json`-`src`-array check already has, so `pubtable_quaids.src`'s
  documented procs don't spuriously fail), and every linked
  command-reference page must exist. **Verified this actually catches a
  real problem, not just that it runs** — the same "never trust a check
  without testing it" standard as everywhere else in this project:
  temporarily renamed one documented link to a nonexistent proc name,
  confirmed the script threw a clear, correctly-identifying error, then
  reverted.
- `scripts/verify_release_artifact.ps1`'s `requiredEntries` list now
  includes `README.md` and all four top-level `docs/*.md` files (previously
  deliberately omitted, per that script's own Milestone 7 comment, since
  they didn't exist yet).

**Verification**: rebuilt the release artifact (`scripts/build_package.ps1
-Force`), confirmed `docs/` and its 18 command-reference pages are
actually present inside the `.zip` (not just referenced by the build
script's directory list), then re-ran the full pipeline
(`scripts/run_release_verification.ps1 -InstallArtifact`) end to end — all
7 source-tree test files (150 checks), the extended manifest check, a real
reinstall to `c:\gauss26\pkgs\quaids`, and `tests/package_public_api.e`
against that install all still pass.

**No version bump**: documentation doesn't change GAUSS public API
surface, so per this repo's established policy the package stays at
`0.5.0`.

## Milestone 9: final gold standard integration gate

A genuine gate, not a rubber stamp: running it top-to-bottom surfaced
three real, previously-undetected gaps, each found by actually exercising
the system (not by re-reading it) — the same standard applied at every
prior milestone.

**Examples didn't match what the docs promised.** README.md,
`docs/USAGE_GUIDE.md`, and every command-reference page (all written at
Milestone 8) document `library quaids;` as the primary usage pattern,
since the package became genuinely installable at Milestone 7. But
`examples/quaids_example.e` and `examples/pubtable_export_example.e`
still `#include`d the source tree directly — a leftover from before
Milestone 7. Fixed: both now use `library quaids;`.
`pubtable_export_example.e` additionally needs a bare `#include
quaids.sdf` — confirmed empirically (not assumed) that `library quaids;`
alone does not activate `pubtable_quaids.src`'s `#ifDef
QUAIDS_SDF_INCLUDED` guard, since `library` lazily compiles individual
procs on demand rather than eagerly running `quaids.sdf`'s `#define`;
explicitly including the `.sdf` (resolved via the installed package's
search path) does activate it — matching exactly what `pubtable`'s own
bundled `pubtable_qardl.src` documents as required for `qardl.sdf` ("just
include qardl.sdf first"). Both examples re-verified against the freshly
rebuilt/reinstalled package.

**The installed-package gate didn't exercise two of its own required
procs.** Auditing "package exports match public docs" line by line found
`tests/package_public_api.e` never actually called `printQuaids()` or
`quaidsElas()` — both real, required, `package.json`-listed public procs,
previously only exercised indirectly via their split components
(`quaidsFit`, `quaidsElasFit`+`printQuaidsElas`). A load-order bug or
stale `.lcg` entry specific to either one could have passed this gate
undetected. Fixed by adding direct calls to both; re-ran the full pipeline
to confirm.

**Published-data validation only covered one of three model families.**
Extended `tests/quaids_published_validation_test.e` (11 -> 19 checks) to
also validate iterated AIDS (`aCtl.linear=1, aCtl.maxiter>1`) against R
`micEconAids::aidsEst(method="IL", ...)` — the Iterated Linear Least
Squares Estimator (Blundell & Robin 1999), which uses the same
starting-value/iteration structure `quaidsFit()` does (LA-AIDS starting
point, then iterate with the translog price index). Observed max absolute
difference ~0.11 (vs. ~0.021 for the existing LA-AIDS/3SLS check),
tolerance set to `0.15`. The wider gap has a real, understood cause, not
approximation slop: `micEconAids`'s `method="IL"` does not support
instrumental variables — confirmed by direct testing, not assumed:
combining `method="IL"` with `instNames` **segfaults** R's `aidsEst`
rather than erroring cleanly. So the R reference here is SUR-only, while
GAUSS's iterated fit always instruments log total expenditure — the
comparison spans both a different algorithm and an IV-vs-no-IV
difference. **QUAIDS has no independent reference implementation
available**: `micEconAids` does not implement a quadratic log-expenditure
term at all, and no other comparably-established QUAIDS implementation
exists (consistent with the Milestone 3 search that produced only a
from-scratch Python replica, itself kept as supplementary evidence, not a
QUAIDS reference). QUAIDS's validation remains the known-true
synthetic-DGP recovery in `tests/quaids_synthetic_validation_test.e` — a
real, non-circular check, but a different, weaker tier of evidence than
cross-implementation agreement on real published data. Documented
explicitly in `docs/FEATURE_SUPPORT_MATRIX.md` rather than silently
equated with the other two model families.

**Two stale Definition-of-Done aspirations reconciled against the actual,
deliberate, already-shipped design**, in `GOLD_STANDARD_TODO.md`: (1) "LA-
AIDS, iterated AIDS, and QUAIDS... with and without IV" assumed an
exogenous-total-expenditure mode that was deliberately never built —
`instr` has been a required argument since Milestone 2/3 by design, not
oversight. (2) "Homogeneity, symmetry, and overidentification each have
standalone... procs" is literally false for two of the three: only
homogeneity and joint homogeneity+symmetry are separately-callable procs
(`quaidsHomogeneityTest`/`quaidsJointTest`); symmetry-given-homogeneity
and overidentification are automatically computed and reported as part of
`quaidsFit()` itself (`qOut.symStat`/`qOut.symPval`,
`qOut.overidValid`/...), a deliberate Milestone 4 scoping decision. Both
aspirations were corrected to match what was actually built and tested,
rather than left as checkboxes nothing could ever satisfy.

**Version consistency verified directly**: `package.json` (`0.5.0`), the
rebuilt artifact filename (`quaids 0.5.0.zip`), the installed copy's
`package.json` (`0.5.0`), and `CHANGELOG.md`'s `## 0.5.0` entry all agree
— read each one directly rather than trusting
`verify_release_artifact.ps1`'s own pass alone, since this is the final
gate. `CHANGELOG.md`'s `0.5.0` entry gained `### Changed`/`### Fixed`
sections documenting Milestones 7-9's tooling/docs/test work (no version
bump -- no `src/` public API changed across any of it).

**Full pipeline re-run end to end after every fix above**: all 7
source-tree test files (158 checks total — schema 34, formula parity 17,
synthetic validation 22, published validation 19, hypothesis tests 19,
elasticities 17, pubtable adapter 30), the extended manifest check, a
real rebuild/reinstall to `c:\gauss26\pkgs\quaids`,
`tests/package_public_api.e` against that install, and both examples run
directly against the installed package.

## Milestone 10: curvature imposition

```gauss
library optmt, quaids;

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;          // Milestone 10 example; QUAIDS support added at Milestone 13
aCtl.maxiter = 100;
aCtl.homogenous = 1;      // required -- quaidsCurvatureFit needs a
                          // homogeneity+symmetry-constrained starting fit

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
call printQuaidsCurvature(cOut);
```

Requested directly by the repo owner after Milestone 9 closed ("are there
any next-level extensions?" → curvature imposition, already flagged as
P2/deferred since Milestone 0's original scoping). Planned via
`EnterPlanMode`/`ExitPlanMode` before writing any code, with two scope
questions put to the repo owner first: AIDS-only vs. also QUAIDS in the
same pass (chose AIDS-only), and full Diewert-Wales nonlinear
reparametrization vs. a cheaper eigenvalue-clipping heuristic (chose full
Diewert-Wales, for real delta-method standard errors rather than informal
point estimates).

**The math**: `quaidsSlutzky()` already computes, per observation, the
matrix that must be negative semidefinite for concavity:
`wepc = -diag(w) + w*w' + gama + (beta'beta)*lx` (AIDS: no quadratic
term). Concavity cannot be imposed globally for a flexible functional
form without over-restricting it (a standard demand-theory result, not a
gap here) -- Diewert & Wales (1987) impose it locally, at the **observed
sample mean** (matching their own practice and `quaidsElasFit()`'s
"evaluate at a given point" convention). Reparametrization: gamma's
upper-left `(n-1)x(n-1)` block becomes `-A*A' - K0` (`A` lower triangular,
`K0` the matrix's non-gamma part at the reference point) -- negative
semidefinite by construction, for any `A`.

**A real design refinement found during implementation** (within the
approved plan's intent, not verbatim): rather than searching the *full*
parameter vector via `optmt`, for *any* candidate `A` the remaining
coefficients are *exactly identified* by OLS once `gama(A)` is substituted
in as a fixed offset (profiled/concentrated nonlinear least squares) --
so `optmt` only searches `vech(A)` (6-15 free parameters typically),
reusing `quaidsFit()`'s own `moment()`/`solpd()` primitives for the "given
`A`" regression, within an outer iteration mirroring `quaidsFit()`'s own
translog-price-index iteration.

**Two real numerical-methods findings from building and testing this, not
from the algebra alone**:

1. **The synthetic fixture's reference point couldn't be an idealized
   population value.** Building the true, curvature-consistent gamma
   against an idealized reference point (uniform shares, or the DGP's
   analytically-derived population mean) left a persistent gap of several
   tenths to units between the constructed gamma's curvature property and
   the *actual* simulated sample's behavior -- confirmed empirically (not
   assumed), and unaffected by the Cholesky factor's scale. Root cause:
   the translog price index's own nonlinearity in gamma means the
   reference point isn't fixed independent of gamma. Fixed with a short
   self-consistent fixed-point iteration (simulate → recompute the
   reference point from the realized sample → rebuild gamma → repeat),
   converging to ~1e-16 agreement within 5-8 rounds for a well-behaved
   seed and diverging to astronomical magnitudes within a handful of
   rounds for most seeds -- the same iteration-sensitivity already
   documented for the main estimator (Milestone 3). `seed=500` was found
   by direct screening, not guessed, and separately confirmed to be a
   non-vacuous test case: the existing unconstrained `quaidsFit()`
   recovers this DGP reasonably (max abs diff ~0.16) but its own Slutzky
   matrix at the sample mean has a positive eigenvalue (~+0.17) -- a real
   violation the constrained fit then has to fix.
2. **Standard errors expose a genuine boundary-inference complication,
   not a bug.** The estimated Cholesky factor frequently converges with
   some entries at exactly zero -- the constrained optimum sits on the
   *edge* of the negative-semidefinite cone, not its interior -- where
   classical delta-method inference is known to be unreliable (the same
   complication that arises for non-negativity-constrained variance
   components elsewhere in econometrics). Confirmed directly (this run's
   own estimated factor has this property), not assumed, and documented
   everywhere standard errors are discussed rather than silently
   presented with false precision. Point estimates and the exact
   curvature property at the reference point are unaffected.

**What GAUSS already provides, reused rather than duplicated**: the
`optmt` package (`c:\gauss26\pkgs\optmt`), exactly as flagged in this
repo's very first roadmap draft. No solver was hand-rolled.
`package.json`'s `deps` array is no longer empty -- this library's first
real external package dependency.

**Version bump to `0.6.0`**: unlike Milestones 7-9 (pure tooling/docs/
testing), this milestone adds real, required new public API
(`quaidsCurvatureFit`/`printQuaidsCurvature`) and a new dependency,
matching this project's policy of bumping the version when public API
surface changes.

**Update, Milestone 13 (2026-07-25)**: the "QUAIDS not yet supported"
scoping above was resolved -- `quaidsCurvatureFit()` now accepts QUAIDS
fits too. See "Milestone 13: QUAIDS curvature imposition" below.

## Milestone 11: welfare measures

```gauss
n = qOut.n;
nint = qOut.nint;
intcptPt = meanc(qOut.intcptFull);
pricesPt0 = meanc(prices);
totexpPt0 = meanc(totexp);

pricesPt1 = pricesPt0;
pricesPt1[1] = pricesPt1[1] + ln(1.05);   // a hypothetical 5% price increase on good 1

struct quaidsWelfareOut wOut;
wOut = quaidsWelfareFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt0, pricesPt1, totexpPt0, aCtl);
call printQuaidsWelfare(wOut);
```

Requested by the repo owner after Milestone 10 closed ("what's next on
the expansion dev path?"). Unlike curvature imposition, this needed no
scoping question up front: computing CV/EV requires no new estimation,
just a closed-form evaluation of the already-fitted expenditure function
at two points, so the same formula either works correctly for both AIDS
and QUAIDS or it doesn't -- no "AIDS first, QUAIDS deferred" tradeoff to
make.

**A real correctness risk was found and resolved *before* writing any
code, not after** -- the same standard as every "verify before trusting a
derived formula" moment in this project's history (Milestone 3's Stone-
index bug, Milestone 5's homogeneity-test formula), just caught one step
earlier this time: an initial from-memory attempt to invert Banks-
Blundell-Lewbel (1997)'s QUAIDS indirect utility function into a closed-
form expenditure function was checked via Shephard's lemma
(`w_i = d ln(e)/d ln(p_i)`, holding utility fixed) against the *already
validated, already tested* QUAIDS share equation used everywhere else in
this codebase -- and did **not** reproduce it. The bug was a misreading of
which term was inverted (`b(p)/L` vs `L/b(p)`) inside the indirect utility
formula itself; the base formula was independently confirmed via a web
search before re-deriving. The corrected formula --

```
ln V(x,p) = 1 / ( b(p)/L(p,x) + lambda(p) )      // indirect utility
ln e(u,p) = a(p) + b(p) / ( 1/u - lambda(p) )     // this milestone's verified inversion
```

-- reproduces the QUAIDS share equation exactly term-by-term via
Shephard's lemma, and collapses correctly to the simpler, independently-
verified AIDS expenditure function `a(p) + u*b(p)` when `lambda=0`. No
external reference data was needed to validate the final implementation:
`tests/quaids_welfare_test.e` checks three *exact* algebraic identities
(zero price change gives `CV=EV=0`; feeding `e(p1,u0)` back into the
*original* indirect utility formula at `p1` returns `u0` exactly, the
defining inverse-function property; `CV`/`EV` for a small price change
converge to the standard Marshallian first-order approximation) rather
than comparing against a published or cross-implementation number,
mirroring Milestone 5's elasticities-identity approach.

Delta-method standard errors reuse the exact numerical-Jacobian-
perturbation technique already validated twice
(`quaidsElasFit`/`quaidsCurvatureFit`) -- no new SE machinery, and none of
Milestone 10's boundary-inference complications apply (no constrained/
reparametrized optimization here, just closed-form evaluation).

**No new package dependency**: pure closed-form algebra. `package.json`'s
`deps` stays at `["optmt"]` (from Milestone 10 only). Version bumped to
`0.7.0` (real new required public API), matching established policy.

**Version control**: this session discovered, via `git log` (not via any
tool call in this conversation), that this repo already has an automated
commit/push process that had captured Milestones 8-10 (`git log` showed
commits with content-accurate auto-generated messages; `local master` and
`origin/master` were confirmed in sync before this milestone started).
The repo owner also asked, starting this milestone, for explicit
commits/pushes at milestone breakpoints going forward -- done for this
milestone (see the repo's commit history), superseding the "nothing in
this repo has been committed" language in every prior milestone's
write-up, which was accurate at the time it was written but is no longer
current status.

## Milestone 12: numerical reliability of the iterated estimator

```gauss
library optmt, quaids;   // optmt only needed if also using quaidsCurvatureFit

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.relax = .75;    // optional -- default 1 (no damping); .75 measurably
                     // reduced the convergence-failure rate in testing

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
print "converged:" qOut.converged " iterations:" qOut.iterations;
```

Requested by the repo owner after Milestone 11 closed ("what's next on
the expansion dev path?" -> recommended closing the loop on a
reliability problem documented but never root-caused since Milestone 3:
the iterated estimator fails to converge, or converges to a wrong
answer, for a large fraction of random seeds, previously described only
as "roughly half" with no surviving evidence for that number -- a full
git-history search confirmed the original 8-seed probe never existed as
a committed artifact anywhere in this repo).

**A real, committed diagnostic, not a one-off probe**:
`tests/quaids_convergence_sweep.e` (run via
`tests/run_convergence_sweep.ps1`) sweeps 200 seeds x 2 models (iterated
AIDS, QUAIDS) at the original informal probe's exact settings
(`tobs=3000`, `aCtl.err=.0001`, `aCtl.maxiter=100`), classifying every
fit into one of three buckets the old "roughly half" phrasing conflated:
**never-converged** (hit `aCtl.maxiter` without meeting tolerance),
**converged-but-wrong** (`qOut.converged==1` but the recovered
coefficients are >10x the normal structural tolerance from the truth --
a self-consistent but wrong fixed point, a materially different failure
mode from simply running out of iterations), and **converged-correctly**.
Deliberately a diagnostic report generator, not a pass/fail gate -- it
prints no "ALL N CHECKS PASSED" line and is not run by
`run_source_tests.ps1`, since there is no convergence guarantee to gate
on and a metric with no such guarantee makes a flaky CI check.

**Measured baseline** (default settings, `aCtl.relax=1`): iterated AIDS
never-converges 39% of the time and converges wrong another 19% (58%
combined failure); QUAIDS never-converges 54.5% and converges wrong
another 21.5% (76% combined) -- both noticeably worse than the old
"roughly half" estimate, especially for QUAIDS.

**A real crash, found by running the sweep at scale, not a
hypothetical**: partway through the very first 200-seed run,
`quaidsFit()` itself threw `error G0121: Matrix not positive definite`
and aborted the whole batch. Root cause: an unguarded `invpd()` in the
symmetry-test block of `src/quaids.src`
(`vi = invpd(v[1:n1*ng, 1:n1*ng])`) crashes the entire call -- not just
returns a bad answer -- whenever a badly-diverged iterated fit produces a
non-positive-definite variance block. Verified the fix in isolation
first (a small reproduction script against the exact failing seed)
before touching production code, using GAUSS's own `trap`/`scalmiss`
idiom -- the identical pattern GAUSS's shipped `gdaols.src` uses for this
exact situation (`oldtrp = trapchk(1); trap 1,1; x = invpd(m); trap
oldtrp,1; if scalmiss(x); ... endif;`). Applied *inside* `quaidsFit()`
itself (not left to callers to work around): on failure, `qOut.symValid`
is forced to `0` and the symmetry-constrained refit falls back to the
homogeneity-constrained estimate as "best available" -- the same
fallback the code already uses when `aCtl.homogenous==0`. Confirmed
`printQuaids()` already gates the entire symmetry-constrained report on
`qOut.symValid` (read directly, not assumed), so the degenerate fallback
values are never surfaced to a human reader. Treated as an unconditional
bugfix (no version bump on its own, no opt-in flag) -- a crash on
unlucky-but-valid input is unambiguously a bug, same precedent as the
Milestone 3 Stone-index fix.

**Fix candidate (a): near-zero-denominator guard**
(`src/quaids.src`'s convergence check,
`err = maxc(maxc(abs((b-b0)./(b0 + (b0 .== 0)))))`, matching the
identical guard `src/quaidscurvature.src`'s own analogous outer-loop
check already uses). Real motivation: this codebase's own synthetic
fixtures build true coefficients via `round(rndns(...)*10)/10`, so exact
zeros are a frequent, reproducible occurrence, not an edge case.
**Verified by re-running the sweep before/after -- and found to have
zero measurable effect**: the post-fix sweep output was byte-identical to
the pre-fix baseline across all 400 fits. Kept anyway (a legitimate,
low-risk defensive fix against a literal zero denominator) and documented
as an honest non-result rather than omitted or oversold -- the same
"verify before trusting a derived fix" standard applied since Milestone
3's Stone-index bug. Likely explanation: the *true* DGP parameter can be
exactly zero, but the *estimated* `b0` mid-iteration is a continuous,
noisy value that essentially never lands on literal `0.0` in floating
point, so the guard's trigger condition rarely if ever fires in practice.

**Fix candidate (b): optional damping** (`aCtl.relax`, new
`quaidsControl` field, `(0,1]`, default `1` -- byte-identical to every
prior release unless a caller opts in). One new line in the iteration
loop, applied *before* the convergence check so `qOut.converged`/
`iterations` reflect the coefficients actually returned:
```gauss
b = solpd(gw, gg);
b = aCtl.relax*b + (1-aCtl.relax)*b0;
err = maxc(maxc(abs((b-b0)./(b0 + (b0 .== 0)))));
```
Verified by re-running the sweep at `aCtl.relax` in
`{1.0, .75, .5, .3}`: correct-convergence rate improves modestly at
`.75` (iterated AIDS 42%->43%, QUAIDS 24%->26.5%), is roughly a wash at
`.5`, and is clearly worse at `.3` (iterated AIDS down to 29%, QUAIDS
down to 19%). Damping substantially reduces the never-converged rate at
every setting tested, but mostly by converting "never converges" into
"converges to a wrong answer" rather than into a correct fit -- it does
not change which fixed point a bad price draw's iteration falls into,
just whether it settles at all. `relax=.75` is a real, modest,
evidence-backed improvement, documented as exactly that and no more.
`tests/quaids_reliability_regression_test.e` pins a concrete example:
seed 2 (QUAIDS) never-converges at the default `relax=1` but converges
correctly in 78 iterations at `relax=.75`.

**Honest scope, matching this project's established standard for partial
fixes** (the boundary-inference caveat on `quaidsCurvatureFit`'s standard
errors, the "no QUAIDS reference implementation" caveat on published-data
validation): this milestone does **not** claim the instability is
solved. Naive successive-substitution on this nonlinear FGLS system can
genuinely have multiple fixed points for a bad price draw -- no amount of
denominator-guarding or mild damping changes which basin of attraction
the iteration falls into. The exit criteria here are "characterized
precisely with real, reproducible numbers, fixed the one confirmed crash
bug, and added one modest, evidence-backed opt-in mitigation," not
"convergence guaranteed."

**Regression testing**: `tests/quaids_reliability_regression_test.e` (8
checks) -- `aCtl.relax=1` byte-identical to leaving it unset; the
previously-crashing seed no longer crashes and correctly reports
`qOut.converged==0`/`qOut.symValid==0`; the seed-2 damping example above
pinned exactly. The full existing 9-file suite (`seed=204` synthetic
fixtures, Blanciforti86 published-data cross-check, and all others)
re-ran clean with unchanged tolerances after all three changes -- that is
the evidence that nothing else moved; no separate pinned "golden" numbers
from before this milestone were needed.

**Version bump to `0.8.0`**: `aCtl.relax` is a new field on the public
`quaidsControl` struct -- real new public API surface, matching the
Milestone 10/11 precedent of bumping on any new public struct field or
proc regardless of whether its default changes existing behavior. (The
crash fix and denominator guard alone would not have warranted a bump,
per the Milestone 3/7-9 precedent for pure bugfixes with no new API
surface.)

## Milestone 13: QUAIDS curvature imposition

```gauss
library optmt, quaids;

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;         // QUAIDS -- supported since Milestone 13
aCtl.maxiter = 100;
aCtl.homogenous = 1;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

aCtl.relax = .25;    // effectively required for QUAIDS's curvature loop
                     // (see below) -- not needed for AIDS

struct quaidsCurvOut cOut;
cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
call printQuaidsCurvature(cOut);
```

Requested by the repo owner after Milestone 12 closed ("what's next" ->
recommended closing Milestone 10's own deferred QUAIDS scope: extending
`quaidsCurvatureFit()`, which previously hard-errored on any QUAIDS fit,
to accept QUAIDS too).

**The simplification, verified before implementation, not assumed**: the
Milestone 10 scoping note framed QUAIDS as "entangling three nonlinear
parameter blocks instead of two" -- true of the algebra, but the apparent
entanglement resolves via the *same* lag-then-solve trick `quaidsFit()`'s
own main iteration loop already uses everywhere else in this codebase.
Confirmed by direct code reading (not assumed) against both
`src/quaids.src` (the main loop reads `_beta` from the *previous* round's
`b` to build `b_p`/`lx2`, then re-estimates a *fresh* `beta`/`lambda`/
`gamma` jointly via one `solpd()`) and `src/quaidscurvature.src` (which
already lags `beta` by one outer round when building `K0` for AIDS, the
identical pattern). `a_p`/`lx` have no `lambda` dependence at all, so the
existing deflator plumbing needed no change -- `beta`/`lambda` never join
`optmt`'s searched parameter set, which stays `vech(A)`-only, unchanged
in size, for either model.

**Exact algebraic identity check, before touching any production code**
-- the same "verify before trusting a derived formula" discipline as
Milestone 11's welfare-formula check: hand-transcribed the K0-split
version of `quaidsSlutzky()`'s QUAIDS formula (`wepc = -diag(w)+w*w'+gama
+(beta'beta+beta'mu+mu'beta+2*mu'mu)*lx`, `mu=lambda*lx/b(p)`) and
confirmed it exactly reproduces `quaidsSlutzky()`'s own formula (copied
verbatim for the comparison side, not re-derived, to avoid making the
same transcription error twice) -- match to ~3.5e-15, floating-point
noise.

**Implementation** (`src/quaidscurvature.src`): `_quaidsCurvGivenA()`
gained an `lx2Fixed` regressor column (built from the previous round's
`beta`, mirroring `quaids.src:279-280`'s own `b_p`/`lx2` construction) so
the same one-shot OLS solve now also estimates `lambda`; the outer loop's
`K0`/warm-start construction gained the `mu`-cross-term, using the
previous round's `beta`/`lambda`. `_quaidsCurvRecoverFull()` gained a
third output (`lambda`, adding-up recovered the same way `beta` already
is). The Hessian's dimensionality (over `vech(A)`) stayed unchanged --
confirmed, not assumed, since this was the design's most load-bearing
"genuinely doesn't need to change" claim.

**A real numerical-instability finding, found only by running it, not by
re-reading the algebra**: the first real attempt (a generic QUAIDS
fixture, `seed=204`, default `aCtl.relax=1`) diverged -- `beta`/`gamma`
grew geometrically each outer round (0.36->2.1->6.5->66->NaN within 8
rounds) before crashing. Root cause: `b(p)=exp(prices*beta')` inside the
lagged `lx2Fixed`/`mu` construction creates a genuinely stronger feedback
loop than AIDS's own `beta'beta` term (no exponential amplification
there). Applying `aCtl.relax` (Milestone 12) to the curvature loop's own
`alpha`/`beta`/`lambda` update -- deliberately *not* `gama`, which must
stay exactly `-A*A'-K0` by construction, or the reparametrization's whole
guarantee breaks -- stabilized it; `relax=.25` converged cleanly
(`maxEig~8.5e-8`) but needed ~185 outer rounds, well past the AIDS-only
`maxOuterIter=50` cap -- bumped to 300 (a no-op for AIDS, which still
converges in ~10-20 rounds either way).

**A second, real, pre-existing bug found along the way, not introduced
by this milestone**: `_quaidsCurvRecoverFull()`'s adding-up recovery
applied the CONSTANT row's "sums to 1" formula (`1 - sumc(...)`) to
*every* intercept-related row, including any extra demographic-shifter
rows (`nint>0`), which must instead sum to 0 (a shifter reallocates
shares, it doesn't change their total -- the same distinction
`tests/quaidsfixtures.src`'s own `al`/`al1` construction already draws).
Latent since Milestone 10, never triggered because the AIDS-only
curvature fixture deliberately has `nint=0` -- surfaced only once a
`nint>0` QUAIDS fixture was tried, producing `sum(wPt)~2.99` instead of
`1` and a spurious `+23` Slutzky eigenvalue that took real isolation work
to trace to this one line (a staleness hypothesis and a "does NSD extend
from the block to the full matrix" hypothesis were both tested and ruled
out first, via direct diagnostic prints, before the actual cause was
found). Fixed with a one-line change; verified as a pure no-op for
`nint=0` (the AIDS fixture's exact case) -- zero regression risk.

**The QUAIDS-analog curvature-consistent fixture was attempted but not
shipped, documented honestly rather than forced**: mirroring
`_quaidsCurvatureSyntheticDGP()`'s AIDS-only self-consistent fixed-point
construction (extended with the `mu`-cross-term in its own per-round
`K0`), a broad screen (tens of thousands of seed/`aScale`/starting-point
combinations) found dozens of seeds where the construction is numerically
self-consistent *and* genuinely negative-semidefinite (`maxEig~1e-16`)
-- but **every single one** implied economically nonsensical mean budget
shares (large negative entries, entries exceeding 1), confirmed to be a
structural property of this fixed-point map for this DGP family (tried
multiple `aScale`/`lambda`-scale/starting-point variants), not a
"wrong seed" problem -- unlike AIDS's own construction, which found a
working seed (`500`) after screening "dozens," not tens of thousands.
`tests/quaids_curvature_test.e`'s QUAIDS block instead validates against
this library's own already-validated general QUAIDS fixture
(`_quaidsSyntheticDGP(seed=204)`), checking convergence, exact negative-
semidefiniteness at the reference point, non-vacuousness, and shape/
finiteness -- a real, deliberately weaker tier of evidence than the AIDS
block's "recovers a known true gamma" check, documented as such rather
than silently equated with it, matching this project's established
honesty standard for partial results (the curvature-SE boundary-
inference caveat, QUAIDS's own "no independent reference implementation"
published-data-validation gap).

**Testing**: `tests/quaids_curvature_test.e` extended in place (17 -> 31
checks) -- no shared AIDS/QUAIDS helper was forced, since the two
blocks' checks differ in kind (true-gamma recovery vs. convergence/NSD/
shape). `tests/package_public_api.e` gained a third inline dataset
(`seed=204` QUAIDS) exercising the QUAIDS curvature path against the
real installed package. The full existing 10-file source-tree suite
re-ran clean with unchanged tolerances.

**Version bump to `0.9.0`**: `quaidsCurvatureFit()`'s signature is
unchanged (no new proc, no new required `quaidsControl` field -- it
reuses `aCtl.relax` from Milestone 12), so this doesn't meet the letter
of the established "new proc/new field" bump policy. Bumped anyway,
since QUAIDS curvature support going from unsupported (hard error) to
supported is a real, user-facing capability addition, not a bugfix.

## Milestone 14: continuous integration

```yaml
# .github/workflows/tests.yml
on:
  push:
    branches: [master]
  workflow_dispatch:
jobs:
  test:
    runs-on: [self-hosted, Windows, gauss]
    steps:
      - uses: actions/checkout@v4
      - name: Run source-tree test suite
        shell: cmd
        run: powershell -ExecutionPolicy Bypass -File tests\run_source_tests.ps1 -SkipBootstrap
```

Requested by the repo owner after Milestone 13 closed ("what's next" ->
offered bootstrap standard errors for `quaidsCurvatureFit()` or CI; repo
owner chose CI first, bootstrap SEs as an explicit follow-up -- see
"Milestone 15" below).

**The one fact that shapes this milestone**: this repo's tests require
GAUSS (`tgauss.exe`), a commercial, licensed product installed only on
the maintainer's machine -- GitHub-hosted runners cannot run them, so CI
requires a **self-hosted runner on this machine**, confirmed with the
repo owner rather than assumed, given the standing system-change/blast-
radius considerations for anything installing persistent background
services.

**A real, load-bearing fact found during planning**: `gh repo view`
confirmed this repository is **public**. GitHub's own guidance is
explicit that self-hosted runners on public repos are a real security
risk -- a fork could open a pull request that runs arbitrary workflow
code on the runner machine. Mitigation: the workflow triggers on `push`
to `master` only (plus manual `workflow_dispatch`), **never**
`pull_request` -- push access is restricted to repo collaborators
(currently just the repo owner), closing the fork/PR attack vector
entirely. Branch protection requiring status checks to pass was
considered and explicitly **not** enabled: doing so would require the
workflow to also trigger on `pull_request` to gate merges, reopening
exactly the risk the push-only design avoids -- a genuine design
tension, not a simple add-on, surfaced when the repo owner asked about it
directly.

**Runner setup**: a GitHub Actions runner (v2.336.0) registered at the
repo level via a short-lived token (`gh api repos/.../actions/runners/
registration-token -X POST`), installed and started as a Windows service
via `config.cmd --runasservice` (not the older, separate `svc.cmd`/
`RunnerService.exe install` two-step, which this runner version does not
support the same way -- confirmed via a web search after `RunnerService.exe
install` produced a Windows Forms dialog rather than installing anything).
Installing a Windows service requires an elevated shell the assistant's
own shell did not have (confirmed via a `WindowsPrincipal`/
`IsInRole(Administrator)` check) -- this one step was handed to the repo
owner to run themselves in their own elevated PowerShell session, rather
than attempting any self-elevation or bypass.

**A real, non-obvious bug found and fixed by reading the actual failed
run's log, not by guessing**: the first real CI run failed with
`PSSecurityException: UnauthorizedAccess -- running scripts is disabled
on this system` (the runner service's account, `NT AUTHORITY\NETWORK
SERVICE`, has a restrictive default PowerShell execution policy). Adding
`-ExecutionPolicy Bypass` directly inside the workflow's `run:` block
under `shell: powershell` did **not** fix it -- confirmed by re-running
and seeing the identical error. Root cause, found by examining the
failed run's log closely: GitHub Actions' `shell: powershell` invokes the
entire `run:` block via `powershell -command ". 'tempfile.ps1'"` --
GitHub Actions writes the whole block (including the `-ExecutionPolicy
Bypass` text) into a temp `.ps1` file, then dot-sources THAT file via a
bare, unqualified call with no bypass flag of its own, so the OUTER
invocation hits the policy wall before the inner command is ever reached.
Fixed by switching to `shell: cmd` for the step (no PowerShell execution-
policy layer at all for the outer invocation) with the actual command
being an explicit `powershell -ExecutionPolicy Bypass -File
tests\run_source_tests.ps1` -- now the first and only PowerShell
invocation, and the bypass applies correctly. Verified via `gh run view
--log` showing genuine execution of all test files with correct PASS
summaries, not a silently-skipped success.

**Scope, deliberately**: only `tests/run_source_tests.ps1` (the fast,
read-only source-tree suite) runs automatically -- not the full release-
verification/rebuild-and-reinstall pipeline, which mutates the shared,
installed GAUSS package directory (`c:\gauss26\pkgs\quaids`) and stays a
manually-run command, consistent with "don't mutate shared machine state
on every push."

**No version bump**: CI tooling doesn't change GAUSS public API surface,
matching this project's established policy (version bumps track public
API changes, not every milestone).

## Milestone 15: bootstrap standard errors

```gauss
library optmt, quaids;

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;          // 1 for AIDS, 0 for QUAIDS (needs aCtl.relax)
aCtl.maxiter = 100;
aCtl.homogenous = 1;

struct quaidsCurvBootOut bootOut;
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, 42);
call printQuaidsCurvatureBootstrap(bootOut);
```

Requested by the repo owner after Milestone 14 closed, per the ordering
established when Milestone 14 itself was scoped ("CI first, then
bootstrap"). Closes a gap `quaidsCurvatureFit()`'s own header has
documented since Milestone 10: its classical delta-method standard error
is known to be unreliable whenever the estimated Cholesky factor sits at
the boundary (a near-zero entry) of the negative-semidefinite cone -- a
frequent, confirmed occurrence in both the AIDS and QUAIDS test fixtures.
A bootstrap (resample, refit, take the empirical spread) is the standard
fix for exactly this class of boundary-constrained-estimation inference
problem.

**A real timing measurement, gathered before designing anything, shaped
the whole milestone**: a single AIDS curvature fit measures ~0.87s (23
outer iterations); a single QUAIDS curvature fit measures ~7.26s (185
outer iterations, needs `aCtl.relax=.25` per Milestone 13's own finding).
At conventional bootstrap replication counts (200-1000), QUAIDS alone
would take 24 minutes to 2 hours. This is why `quaidsCurvatureBootstrapFit()`
has **no default `B`** -- the caller must choose deliberately, confirmed
with the repo owner via explicit scoping questions (also resolved this
way: no default progress printing during the loop, matching this
codebase's own silent-Fit-proc convention and this author's own sibling-
repo bootstrap precedents in `gauss-qardl`/`dccelib`; no in-proc
percentile confidence intervals in this pass, though raw draws are kept
in `bootOut.bBoot` for a caller or a later milestone to use; and the new
automated test stays opt-in via a `-SkipBootstrap` flag rather than
running in the default per-push CI, since it adds ~45-50s -- roughly
doubling `run_source_tests.ps1`'s baseline runtime even at a small `B`).

**Design, informed by this same author's existing GAUSS bootstrap code in
sibling repos** (`gauss-qardl`'s `blockBootstrapQARDLDiag`,
`dccelib`'s `mgBootstrap`/`mgBootstrapSE`): a nonparametric i.i.d. row
(pairs) bootstrap (the correct choice for this cross-sectional data, as
opposed to `gauss-qardl`'s block bootstrap, which exists there for
genuine time-series dependence not present here) -- one draw of row
indices applied identically to `w`/`intcpt`/`prices`/`totexp`/`instr`,
refitting the *whole* pipeline (`quaidsFit()` then `quaidsCurvatureFit()`)
on each resample rather than re-perturbing the already-fitted curvature
struct, so first-stage IV sampling variability is captured too. A
replication counts only if both stages converge and the resulting
coefficient vector is finite -- up to `5*B` total resamples are attempted
before giving up on reaching `B` completed replications, the same
attempt-cap convention `blockBootstrapQARDLDiag` already uses. `intcpt`'s
`== 0` ("no extra intercept shifters") convention is preserved during
resampling (checked directly -- indexing a scalar `0` with row indices
errors, confirmed by running it, not assumed).

**Two real, non-obvious GAUSS failure modes found and fixed while
building this, both by running real stress tests, not by reasoning about
the language spec**:

1. **A struct-returning proc call CAN be safely wrapped in this
   codebase's `trap`/`scalmiss` idiom** (Milestone 12's pattern,
   `src/quaids.src:552-563`) -- confirmed empirically with a small
   isolated probe before relying on it: GAUSS does not abort the whole
   statement when a trapped error occurs inside a called proc: it
   substitutes a scalar-missing value for the erroring expression and
   lets that proc's own execution continue to its `retp()`, so checking
   the specific fields a caller reads afterward (e.g. `converged`) is
   sufficient -- there is no way to `scalmiss()` a whole struct directly
   (it takes a matrix, not a struct). This let both
   `quaidsFit()` and `quaidsCurvatureFit()` be wrapped per bootstrap
   replication (a new private helper, `_quaidsCurvBootOneRep()`) without
   modifying either proc's own signature.
2. **`trap 1` does NOT catch every failure mode** -- a real crash
   (`error G0528: More returns than targets`) surfaced when stress-
   testing the QUAIDS path under resampling, inside
   `quaidsCurvatureFit()`'s own `{ va, ve } = eighv(...)` calls (the
   warm-start and the Hessian-based SE block), **even with the whole call
   wrapped in `trap 1`**. Root cause, confirmed by direct probing (not
   assumed): on a sufficiently degenerate input, `eighv()` itself returns
   *fewer* values than the destructuring assignment expects -- a call-
   arity mismatch, not a trappable numerical error, so `trap` does not
   intercept it. Narrowed the trigger with a focused probe script: a
   zero matrix and a singular (rank-deficient) matrix both still return 2
   values fine; a matrix containing a plain `Inf` does not. A `x .eq x`
   NaN-detection pre-check alone was **not enough**, since `Inf .eq Inf`
   is `TRUE` under IEEE 754 -- a separate `abs(x) < 1e100` magnitude bound
   was needed too. Fixed by pre-checking finiteness (both NaN and Inf)
   before ever calling `eighv()` in both places inside
   `quaidsCurvatureFit()` itself, falling back to a harmless identity
   decomposition (`va=1s, ve=I`) on failure -- a real, general robustness
   improvement to `quaidsCurvatureFit()`, not just a bootstrap-specific
   workaround, since the underlying gap (no internal guard around these
   two calls) predates this milestone and was explicitly flagged as a
   risk during this milestone's own design review.

**A test-design lesson, found while writing the plausibility check, not
assumed going in**: the delta-method SE's boundary-inference
unreliability does **not** stay confined to the specific gamma row/column
tied to a boundary Cholesky entry -- the classical NLS covariance is for
the *whole* `vech(A)` vector at once, and the numerical-Jacobian delta
method mixes every `vech(A)` parameter into every reported coefficient's
SE, so a single boundary entry can inflate `seDelta` anywhere in the
vector (confirmed directly: `seDelta` reached into the hundreds on rows
with no boundary-adjacent Cholesky entry of their own, on the very
fixture this test uses). An initial per-cell "same order of magnitude,
excluding boundary rows" check was written, run, and found to fail for
exactly this reason -- replaced with a plausibility check on the concrete
thing this milestone is actually for: `seBoot` itself stays well-behaved
and bounded where `seDelta` does not.

**Testing**: `tests/quaids_curvature_bootstrap_test.e` (26 checks) --
bootstrap run bookkeeping, shape/finiteness of the bootstrap SE, exact
echo of the base point estimate/delta-method SE, and the well-behaved-SE
plausibility check above, on both an AIDS (`B=15`, ~13s) and QUAIDS
(`B=5`, ~36s) fixture reusing `tests/quaidsfixtures.src`'s existing
`tobs=3000` datasets unchanged (only `B` was kept small, not `tobs` --
shrinking `tobs` risked changing convergence behavior already validated
at that size). Not part of `run_source_tests.ps1`'s default run (see
`-SkipBootstrap` above).

**A third real bug, found only by running the full local suite for
real**: `run_source_tests.ps1`'s `Invoke-GaussBatch` helper (Milestone 7)
read the child `tgauss.exe` process's stdout fully before touching its
stderr -- a classic .NET `Process` deadlock precondition (if a child
fills both OS pipe buffers before either is drained, it blocks mid-write
while the parent blocks reading the other stream, and neither proceeds).
No prior test file produced enough combined output to trigger it;
`quaids_curvature_bootstrap_test.e`'s QUAIDS block does, since a bad
resample routinely makes `optmt` print dozens to hundreds of harmless
"Optmt: function evaluation failed" lines to stderr. Confirmed directly:
running the same file straight via `tgauss -b -x
quaids_curvature_bootstrap_test.e` (a real console, no pipe to fill)
finished in under a minute, while running it through
`run_source_tests.ps1` left the child sitting idle for over eight hours
before the hang was noticed and the process killed. Fixed by draining
both streams asynchronously (`OutputDataReceived`/`ErrorDataReceived` +
`BeginOutputReadLine()`/`BeginErrorReadLine()`) instead of two sequential
`ReadToEnd()` calls -- the standard fix for this well-known .NET pitfall.
Re-ran the full local suite (11 files, no skips) afterward and confirmed
it completes cleanly with no hang.

**Version bump to `0.10.0`**: a new required public proc
(`quaidsCurvatureBootstrapFit`) and a new required public struct
(`quaidsCurvBootOut`), matching this project's established policy of
bumping on real new public API surface. No new package dependency --
pure GAUSS built-ins plus the already-required `optmt`.

## Milestone 16: predicted budget shares

```gauss
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
nint = qOut.nint;
intcptPt = meanc(qOut.intcptFull);
pricesPt = meanc(prices);
totexpPt = meanc(totexp);

struct quaidsSharesOut sharesOut;
sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call printQuaidsShares(sharesOut);
print "sum of predicted shares:" sumc(sharesOut.w);   // exactly 1
```

After Milestone 15 closed, the repo owner asked for a broader outline of
what a "full demand-system workflow" would still need, then asked to work
through the resulting list in order. First item: expose the model's
predicted budget share at an arbitrary point as a standalone public proc
-- useful for out-of-sample prediction and policy simulation, previously
only available by hand-deriving the share equation.

**Not a new formula -- a third, deliberately independent copy of an
existing one.** Direct exploration confirmed the exact share equation
already exists as an intermediate step inside `quaidsElas_()`
(`src/quaidselas.src:49-61`, on the way to computing elasticities), and
was independently duplicated a second time in
`tests/quaids_elasticities_test.e`'s own private `modelShareAt()` helper
(needed only so that test could check its Engel/Cournot/homogeneity
identities against the right `w`). Rather than refactor `quaidsElas_()`
to share code with a new proc, `quaidsSharesFit()` (`src/quaidsshares.src`)
implements its own small private helper (`_quaidsSharesAt()`) with the
identical formula, honestly documented as a third copy -- matching this
project's established preference (Milestone 2's scoping note) against
touching already-shipped, tested code without a strong reason. This was
an explicit scoping choice, confirmed with the repo owner via
`AskUserQuestion` alongside two other choices: the proc/struct naming
(`quaidsSharesFit`/`printQuaidsShares`, matching the existing noun-based
`quaidsCurvatureFit`/`quaidsElasFit`/`quaidsWelfareFit` convention) and
scope (budget shares only, not derived physical quantities -- computing
`q_i = w_i*exp(totexp)/exp(price_i)` would require assuming `prices`/
`totexp` are logs of levels in mutually consistent units, which nothing
else in this library requires of the caller).

**Standard errors differ in kind from `quaidsElasFit`/`quaidsWelfareFit`**:
those two only ever need marginal standard errors, so they accumulate
variances via a cheaper recursive formula. `quaidsSharesFit()` reports the
**full `n x n` covariance** of the predicted share vector (`sharesOut.v`),
built via an explicit numerical Jacobian propagated as `jacW*v*jacW'` --
the same matrix-form delta method `quaidsCurvatureFit()` already uses --
so a caller can correctly test hypotheses spanning more than one good
(e.g. whether two goods' shares differ significantly), not just read off
marginal SEs.

**Testing**: `tests/quaids_shares_test.e` (21 checks) -- the point
estimate matches a *fresh*, independently hand-evaluated version of the
share formula (written directly in the test, not calling any `src/` proc)
on both an AIDS (`nint=0`) and QUAIDS (`nint=1`) fixture; exact adding-up
(`sum(w)==1`, a direct consequence of `qOut.bestB`'s own adding-up
construction); shape/finiteness/non-negativity of `se`/`v`, and that `se`
exactly equals `sqrt(diag(v))`; a shifted evaluation point gives a
genuinely different predicted share (non-vacuous). Additionally,
`tests/quaids_elasticities_test.e`'s private `modelShareAt()` helper was
removed and replaced with direct calls to `quaidsSharesFit()` -- its
existing Engel/Cournot/homogeneity exact-identity checks (which only hold
given the *correct* model-implied `w`) re-ran clean, serving as a real
cross-check that the new proc's point estimate is right, not just that it
runs. `tests/package_public_api.e` gained a call to
`quaidsSharesFit()`/`printQuaidsShares()` too, exercising it through the
installed-package gate.

**Version bump to `0.11.0`**: a new required public proc
(`quaidsSharesFit`) and a new required public struct (`quaidsSharesOut`),
matching this project's established policy of bumping on real new public
API surface. No new package dependency.

## Milestone 17: AIDS-vs-QUAIDS specification test

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;         // QUAIDS -- required, see below
aCtl.homogenous = 0;     // unconstrained -- required, same as the other standalone tests

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

{ statQ, pvalQ, dfQ } = quaidsQuadraticTest(qOut);
if pvalQ > 0.05;
    print "Quadratic term not statistically justified -- consider AIDS.";
endif;
```

Second item of the full-demand-system-workflow outline the repo owner
asked to work through in order after Milestone 15 closed (first item was
Milestone 16's `quaidsSharesFit`). Closes a real workflow gap: a user
fitting QUAIDS had no formal way to check whether the extra quadratic-
term complexity was actually justified by the data, only informal
before/after comparison.

**Mirrors `quaidsHomogeneityTest()`/`quaidsJointTest()` exactly** (same
file, `src/quaidstests.src`; same `proc (3) = ...(struct quaidsOut qOut);`
returning a bare `(stat, pval, df)` tuple, not a struct; same
`qOut.homogenous /= 0` guard and `"ERROR - ..."; stop;` idiom) -- a
standard Wald test of `H0: lambda_i = 0` for every good, `df = n-1` for
the same reason the existing tests use `n-1` (equation `n`'s restriction
is implied by adding-up once the other `n-1` hold). Structurally simpler
than `quaidsHomogeneityTest`'s own `L` construction, since `lambda` is
already a single scalar coefficient per equation (one row per equation's
column block of `vec(qOut.b)`), not a row-sum across several gamma
columns needing summation.

**A second guard, not present on the sibling tests**: this test also
requires `qOut.linear == 0`. Confirmed by direct exploration of
`src/quaids.src` (the `qOut.b`/`qOut.v` population, guarded throughout by
`if not aCtl.linear`) that when `aCtl.linear = 1` (AIDS), `lambda` is not
a nuisance parameter fixed at zero -- it is **not estimated at all**, and
`qOut.b` is one row shorter with no lambda row to select. There is no way
to "test whether QUAIDS is needed starting from an AIDS fit"; this test
only makes sense as a check on an already-fitted QUAIDS model. The row
position itself (`lambdaRow = 1+nint+n+2`, immediately after beta and
immediately before the u-residual rows) was confirmed against
`src/quaids.src`'s own recovery step and independently corroborated by
`quaidsCurvOut`'s doc comment (`src/quaids.sdf`) describing the identical
"trailing lambda row" convention from Milestone 10/13's curvature work.

**A real test-design finding**: GAUSS's `stop` (used by this file's
existing guard idiom) halts the *entire batch job*, not just the current
proc call -- there is no way to assert "this call errors" and then
continue to the next check in the same script, confirmed by a direct
probe before attempting to write such a check. Neither
`quaidsHomogeneityTest`'s nor `quaidsJointTest`'s own guards have a
"confirm it errors" test in `tests/quaids_hypothesis_tests_test.e`
either -- consistent, not an oversight specific to this milestone.

**Testing**: `tests/quaids_hypothesis_tests_test.e` extended in place
(19 -> 22 checks) -- reuses the file's *existing* `qOutTrue` fixture
(`_quaidsSyntheticDGP(3000, 204, 1, 1)`, `quadratic=1`, true `lambda`
genuinely nonzero) directly for the POWER check (no new fit needed), plus
one fresh `quadratic=0` fixture (true `lambda` forced to exact zero) for
SIZE, both fit with `aCtl.linear=0` (required either way, per the guard
above).

**Version bump to `0.12.0`**: a new required public proc
(`quaidsQuadraticTest`), matching this project's established policy of
bumping on real new public API surface. No new struct, no new file, no
new package dependency.

## Milestone 18: percentile bootstrap confidence intervals

```gauss
struct quaidsCurvBootOut bootOut;
bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 200, 42);

{ ciLower, ciUpper } = quaidsCurvatureBootstrapCI(bootOut, 0.05);
print "95% CI for good 1's first coefficient:" ciLower[1,1] ciUpper[1,1];
```

Third item of the full-demand-system-workflow outline being worked
through in order (after Milestone 16's `quaidsSharesFit` and Milestone
17's `quaidsQuadraticTest`). Milestone 15's own header comment had
explicitly scoped this out at the time: "no in-proc percentile confidence
intervals in this pass, though raw draws are kept in `bootOut.bBoot` for
a caller or a later milestone to use" -- this milestone is exactly that
follow-on. No new resampling or refitting needed: `bootOut.bBoot` already
holds every completed replication's coefficient vector.

**New proc, not a signature change to `quaidsCurvatureBootstrapFit()`**:
changing that already-shipped (v0.10.0) proc's signature would be a
breaking change to released public API for a purely additive feature.
`quaidsCurvatureBootstrapCI(bootOut, alpha)` is a separate consumer proc,
returning a bare `(ciLower, ciUpper)` tuple (no struct, no paired
printer) -- mirrors this codebase's own `quaidsHomogeneityTest`/
`quaidsJointTest`/`quaidsQuadraticTest` convention for a derived-quantity
utility proc, not the heavier struct+printer "Fit" convention. `alpha` is
required (no default), matching `quaidsCurvatureBootstrapFit()`'s own
`B`/`seed` philosophy of never silently guessing an inference-affecting
parameter. Mechanics mirror `gauss-qardl`'s own
`blockBootstrapQARDLDiag` quantile-CI idiom (`quantile(boot_beta[.,jj],
q_lo)`, looped per column) rather than inventing a new one.

**A real, significant bug found and fixed, spanning two already-shipped
milestones**: building this proc's own ground-truth cross-check (does
`ciLower`/`ciUpper` actually bracket the point estimate at the *correct*
cell?) surfaced that GAUSS's `reshape()` fills row-major, not
column-major like `vec()` -- confirmed empirically with a small,
synthetic, hand-verifiable example (`reshape(vec(X), rows(X), cols(X))`
does **not** recover `X` in general; the correct inverse is
`reshape(v, cols(X), rows(X))'`, reshape into the *transposed* shape
then transpose back). `src/quaids.src`'s own pre-existing code already
used this correct transpose-based idiom consistently (e.g. `reshape(stderr,
n, ng)'`) -- but `quaidsCurvatureFit()`'s `seAll = reshape(seAll,
rows(bstack0), cols(alpha));` (Milestone 10, shipped since v0.6.0) and
`quaidsCurvatureBootstrapFit()`'s `seBoot = reshape(stdc(bBoot), bDim1,
bDim2);` (Milestone 15, shipped since v0.10.0) both dropped the
transpose when they were written, silently scrambling `cOut.se`/
`bootOut.seBoot`'s individual cells relative to `cOut.b`/`bootOut.b`
ever since. **Invisible to every existing test**: all of them check
shape, sign (non-negative), and finiteness (no NaN) -- properties that
are *permutation-invariant* and would pass identically whether or not
the cells were scrambled. Confirmed directly against real fitted data
(not just reasoned about): `cOut.se[i,j]` did **not** equal
`sqrt(cOut.v[k,k])` at the true vec-order position `k` for `cOut.b[i,j]`,
by a wide margin, until fixed. `cOut.v` (the full covariance) was never
affected, only the reshaped `se`/`seBoot` display convenience derived
from its diagonal. Fixed in all three places (`quaidsCurvatureFit`'s
`seAll`, `quaidsCurvatureBootstrapFit`'s `seBoot`, and this milestone's
own `quaidsCurvatureBootstrapCI` -- caught in the last before it ever
shipped).

**Regression guards added, not just a silent fix**: `tests/
quaids_curvature_test.e` and `tests/quaids_curvature_bootstrap_test.e`
both gained an explicit cell-position check (`cOut.se` vs. a fresh
`reshape(sqrt(diag(cOut.v)), cols(b), rows(b))'` recomputation, and the
same for `seBoot` vs. a fresh `stdc(bBoot)` recomputation) -- checking
the *specific* thing the existing shape/sign/finiteness checks could not
catch, so this exact bug can't silently return.

**Testing**: `tests/quaids_curvature_bootstrap_test.e` extended in place
(26 -> 37 checks) -- `quaidsCurvatureBootstrapCI()` shape/ordering
(`ciUpper >= ciLower`)/containment (`ciLower <= b <= ciUpper`) checks, and
a direct `quantile()` cross-check at one specific flattened index against
the correctly-mapped `(row, col)` cell of `ciLower` -- verifying the
reshape/index mapping itself, not just "it runs". `tests/
quaids_curvature_test.e` gained the `cOut.se` position-correctness
regression guard described above (17 -> 19 AIDS checks; 31 -> 33 total).
`tests/package_public_api.e` gained a call to
`quaidsCurvatureBootstrapCI()` too.

**Version bump to `0.13.0`**: a new required public proc
(`quaidsCurvatureBootstrapCI`), matching this project's established
policy of bumping on real new public API surface -- bundled with the
reshape bugfix in the same release, per this project's practice of not
holding a real correctness fix for a separate release when a version
bump is already warranted.

## Milestone 19: zero budget share correction (Shonkwiler-Yen)

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 0;          // 1 for AIDS, 0 for QUAIDS
aCtl.maxiter = 100;
aCtl.homogenous = 0;      // required -- see below

struct quaidsZeroOut zOut;
zOut = quaidsZeroFit(w, intcpt, prices, totexp, instr, aCtl);
call printQuaidsZero(zOut);
print "fraction of zero shares per good:" zOut.shareZeroFrac';
```

Fourth and, by direct exploration of `quaidsFit()`'s own body before
planning began, the **largest single addition to the library to date**
item of the full-demand-system-workflow outline being worked through in
order (after Milestone 16's `quaidsSharesFit`, Milestone 17's
`quaidsQuadraticTest`, and Milestone 18's `quaidsCurvatureBootstrapCI`).
Real survey/microdata routinely has corner solutions -- households
reporting zero expenditure on some goods -- which the linear/log-linear
AIDS/QUAIDS share equation has no mechanism for; fitting `quaidsFit()`
directly on such data is a known source of bias.

**The architectural obstacle, confirmed before any code was written**:
every stage of `quaidsFit()` (the coefficient GLS solve, the variance
formula, the homogeneity/symmetry minimum-distance restriction, the
overidentification test, the absolute-price recovery) is built on a
Kronecker-product identity (`S[1:n-1,1:n-1].*.gg`, `src/quaids.src`) that
holds only because every equation shares the *same* design matrix `X`. A
literal, textbook Shonkwiler & Yen (1999) correction rescales *every*
regressor in equation `i` by that equation's own first-stage probability
`F_i` -- which breaks the shared-`X` assumption outright and would need
rewriting the Kronecker-based core into a genuine block system. That is
not what this milestone does.

**The reformulation that avoids this, derived and confirmed
mathematically sound before implementation** (the same "verify before
trusting a derived formula" discipline as Milestone 3's Stone-index bug
and Milestone 11's welfare-formula check): Shonkwiler-Yen's equation

```
w_i = F_i*(alpha_i + sum_j gamma_ij*ln(p_j) + beta_i*lx [+ lambda_i*lx2]) + f_i*delta_i + e_i
```

divides cleanly by `F_i` (a *known*, first-stage-fitted quantity, held
fixed during the second stage):

```
w_i/F_i = alpha_i + sum_j gamma_ij*ln(p_j) + beta_i*lx [+ lambda_i*lx2] + (f_i/F_i)*delta_i + e_i/F_i
```

This turns the problematic *regressor* rescaling into a **dependent-
variable transform** (`wTilde_i = w_i/F_i`, computed once, before any GLS
solve) plus **one new shared regressor column per equation**
(`h_i = f_i/F_i`) -- structurally the same kind of addition as the `u`
(IV-residual) columns `quaidsFit()` already appends to its own shared `X`.
The shared-`X` trick, and everything built on it, survives untouched; the
whole translog-price-index outer iteration is reused essentially unchanged
(`src/quaidszerocorrect.src` mirrors `src/quaids.src`'s starting-value and
iteration blocks structurally), just fed `wTilde`/`h` instead of `w`.

**The one real complication the reformulation does not remove**: adding
`n` hazard columns to a shared `X` means the one-shot GLS solve estimates
a full `n x n` cross-equation `delta` block (every equation's response to
*every* good's hazard term), when Shonkwiler-Yen only wants the diagonal
(good `i`'s own hazard term in good `i`'s own equation). `_quaidsZeroDiagRestrict()`
imposes this via the same `design()`/selection-matrix minimum-distance
idiom `quaidsFit()`'s own symmetry-restriction stage already uses
(`src/quaids.src:526-614`), just with a diagonal (not
`gamma_ij=gamma_ji`) restriction pattern, and the same `trap 1,1;`/
`scalmiss()`-guarded graceful degradation (Milestone 12's idiom) if the
restriction's own matrix inversions fail on a badly-behaved fit.

**Decisions, confirmed with the repo owner via `AskUserQuestion` before
implementation**: full Shonkwiler-Yen implementation in one milestone, not
split into sub-steps; adding-up honestly documented as approximate for the
corrected coefficients, not forced exact via post-hoc renormalization
(this is a real, known property of the method itself -- each equation is
independently rescaled by its own good-specific `F_i` -- not a bug to
paper over); probit regressors reuse `intcpt`/`prices`/`totexp` already
passed to the estimator, no new required argument, matching this
library's "no separate exogenous-mode arguments" philosophy.

**Scope of this pass, deliberately limited** (matching the Milestone 10/
13 AIDS-then-QUAIDS curvature precedent for phased scope): **unconstrained
only** -- `quaidsZeroFit()` errors clearly if `aCtl.homogenous = 1`.
Imposing homogeneity/symmetry *on top of* the Shonkwiler-Yen correction is
real, additional work (combining two different minimum-distance
restrictions in the same pass), left for a follow-up. Standard errors use
a simplified `V(vec(b)) = S .*. inv(gg)` formula (the classic SUR-with-
shared-regressors covariance) -- honestly documented as **not** correcting
for the nonlinear translog-price-index feedback the way `quaidsFit()`'s
own Jacobian-corrected variance does, nor for first-stage probit/IV
generated-regressor uncertainty, matching the established precedent for
`quaidsCurvatureFit()`'s own "simplified, not a full sandwich" SE.

**Implementation** (`src/quaidszerocorrect.src`, new file -- a
self-contained sibling to `quaidsFit()`, not a modification of it,
matching Milestone 16's identical "don't touch already-shipped, tested
core code" reasoning): `quaidsZeroFit()` reuses `_quaidsIVFirstStage()`
(`src/quaidsiv.src`) for its own IV first stage, runs `n` independent
probits via GAUSS's **built-in** `glm()` (`c:\gauss26\src\glm.src`, base
runtime, no new package dependency) with `ctl.link="probit"` and
`ctl.constantFlag=-1` (the shared regressor block already carries its own
constant column), computes `F_i` (floored at `1e-3`) and the probit
density `f_i` (via `pdfn()`, since `glm()` does not expose it directly),
builds `wTilde`/`h`, then runs a `quaidsFit()`-mirroring Stone-index
starting value plus translog-price-index iteration loop targeting
`wTilde`/`h`. `_quaidsZeroDiagRestrict()` is applied once, after the main
iteration converges, since the restriction itself is linear and does not
need to sit inside the nonlinear loop.

**Real bugs found and fixed while building this, all via direct empirical
testing** (this project's standing discipline):

1. `struct glmControl gCtl;`/`struct glmOut gOut;` declared *inside* the
   per-good probit `do while` loop threw "Invalid structure redefinition"
   on the second pass -- GAUSS does not allow a struct type to be
   redeclared inside a loop body. Fixed by hoisting both declarations
   above the loop, only reassigning `gOut = glm(...)` inside it.
2. A local variable named `f` (the probit density) triggered "Duplicate
   definition of local 'f'" -- renamed to `fDens` throughout to eliminate
   any possible naming collision.
3. `F[., i] = maxc((gOut.yhat)'|(1e-3*ones(1, nobs)))';` had a stray
   trailing transpose: `maxc()` on a `2 x nobs` input already returns the
   correct `nobs x 1` column, and the extra `'` flipped it to `1 x nobs`,
   throwing "Rows don't match" on assignment into `F[.,i]`. Fixed by
   removing the trailing transpose.

**A known, unresolved limitation, found by running a seed screen, not
theorized**: GAUSS's built-in `glm()` can hard-crash on some inputs with
an uncatchable `Intel MKL ERROR: Parameter 5 was incorrect on entry to
DGELS` / "input 'x' contains missing values or estimation failed" --
confirmed this is **not** trappable via this codebase's usual
`trap 1,1;`/`scalmiss()` guard idiom (the same class of non-trappable
failure already documented for `eighv()`'s call-arity mismatch inside
`quaidsCurvatureBootstrapFit()`, Milestone 15). Some seeds (e.g. seed=2)
trigger it, others (the shipped seed=1 fixture) do not. Time-boxed
decision: pick a working seed rather than hardening `quaidsZeroFit()`
against this failure mode in this pass -- documented as a real, known
limitation in `docs/USAGE_GUIDE.md`/`docs/METHODOLOGY_NOTES.md`, not
silently absent.

**Fixture calibration, screened empirically rather than guessed**:
`tests/quaidsfixtures.src`'s new `_quaidsZeroSyntheticDGP(tobs, seed)`
generates a latent (uncensored) 5-good QUAIDS share the same way
`_quaidsSyntheticDGP()` does, then censors it the economically correct
way -- `w_i = max(0, latent_i) / sum_j max(0, latent_j)`, an accounting
identity (what real survey shares actually are: expenditure floored at
zero, divided by total *actual* expenditure), not an ad hoc
redistribution, so adding-up holds exactly in the *observed* data by
construction. Getting a genuine, non-degenerate censoring rate required
real experimentation: reducing the noise term alone (`_quaidsSyntheticDGP`'s
own scale of `2`, tried at `.12` then `.04`) had almost no effect on
zero-share fractions (confirmed empirically, not assumed); reducing price
variance instead made things *worse*, not better (also confirmed, then
reverted). The combination that worked -- kept full price variance,
scaled the *structural* coefficients (`gamma`/`beta`/`lambda`/`al1`) down
to `0.25x` their `_quaidsSyntheticDGP()` magnitudes -- confirmed that the
**deterministic** `price*gamma`/`expenditure*beta` structural swing, not
noise or price variance, is the dominant driver of negative latent shares
in this DGP family at the original scale. `seed=1` (found by direct
screening 1-30, not arbitrary) gives per-good zero-share fractions
`[0.843, 0.184, 0.817, 0.211, 0.171]` -- genuinely uneven (two goods
heavily censored, three moderately so) but non-degenerate (no good is
always or never zero), and both `quaidsZeroFit()` and, for comparison,
naive `quaidsFit()` converge cleanly on it.

**The core empirical validation, confirmed via a direct seed-level
comparison before writing the formal test suite**: on this fixture,
`quaidsZeroFit()`'s corrected coefficients recover the true *latent*
(uncensored) DGP parameters better than naively fitting `quaidsFit()` on
the same *censored* data, on **both** metrics -- max absolute difference
`2.0262317` (corrected) vs. `2.1677279` (naive); mean absolute difference
`0.40114995` (corrected) vs. `0.40335223` (naive). A real, if modest,
improvement, not a dramatic one -- consistent with Shonkwiler-Yen being a
known approximation to the fully efficient (but far more complex) Amemiya-
Tobin-style censored system estimator, and with the high per-good
censoring rates this specific fixture has. Documented honestly rather
than searching for a seed with a more dramatic-looking gap.

**Testing**: `tests/quaids_zero_test.e` (17 checks) -- the fixture's own
adding-up identity (exact, by construction); the diagonal-delta
restriction holds exactly (off-diagonal hazard-coefficient entries are
exactly `0`, on-diagonal entries are genuinely, nonzero-ly estimated);
shape/finiteness of `probitB`/`se`/`b`; all `n` first-stage probits
converged; and the core validation above (`quaidsZeroFit()` beats naive
`quaidsFit()` on both the max- and mean-absolute-difference metrics
against the true latent DGP). Added to `tests/run_source_tests.ps1`'s
default list (a single fit plus `n` probits, no heavy per-call cost like
the bootstrap tests). `tests/package_public_api.e` gained a fourth inline
dataset (mirroring `_quaidsZeroSyntheticDGP(3000, 1)`, duplicated inline
per this file's own "no dependency on tests/-only fixture code"
principle) exercising `quaidsZeroFit()`/`printQuaidsZero()` against the
real installed package.

**Version bump to `0.14.0`**: a new required public proc
(`quaidsZeroFit`, plus its paired printer `printQuaidsZero`) and a new
required public struct (`quaidsZeroOut`), matching this project's
established policy of bumping on real new public API surface. No new
package dependency -- `glm()` is part of GAUSS's own base runtime.

## Milestone 20: robust and cluster-robust standard errors

```gauss
struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.homogenous = 1;

struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

// Heteroskedasticity-robust:
struct quaidsRobustOut rOut;
rOut = quaidsRobustFit(qOut, w, prices, totexp, aCtl, 0);
call printQuaidsRobust(rOut);

// Cluster-robust (householdId is a Tx1 vector of group labels):
struct quaidsRobustOut rOutCluster;
rOutCluster = quaidsRobustFit(qOut, w, prices, totexp, aCtl, householdId);
call printQuaidsRobust(rOutCluster);

// Cluster-aware bootstrap alternative:
struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, householdId, 200, 42);
call printQuaidsRobustBootstrap(rbOut);
```

Fifth and final item of the full-demand-system-workflow outline being
worked through in order (after Milestone 16's `quaidsSharesFit`,
Milestone 17's `quaidsQuadraticTest`, Milestone 18's
`quaidsCurvatureBootstrapCI`, and Milestone 19's `quaidsZeroFit`). After
this milestone, the originally-outlined workflow is fully implemented.

**Confirmed before any code was written, not assumed**: every existing
covariance in this library (`qOut.homogV`/`symcV`, the IV first stage,
`quaidsCurvatureFit()`'s `cOut.v`, `quaidsZeroFit()`'s `zOut.v`) rests on
a single pooled, homoskedastic `S = sse/nobs` combined with the shared-
design-matrix `S.*.inv(gg)` Kronecker sandwich -- no heteroskedasticity-
robust or cluster-robust computation exists anywhere in this library.
Direct exploration of GAUSS's base runtime (`c:\gauss26\src\robust.src`:
`robustSE`/`clusterSE`/`hacSE`) and the installed `tsmt` package's near-
identical procs found neither generalizes to this library's **stacked
multi-equation system with a shared regressor matrix** -- both are
single-equation `(X'X)^-1 (...) (X'X)^-1` sandwiches for one dependent
variable and shared regressors. Reusing them would mean unpacking into
exactly the same per-cluster score-aggregation math this milestone builds
directly, for zero net simplification -- the same conclusion this project
reached about `gmmFitIV` at Milestone 2. This is genuinely new math, not
a reuse-an-existing-utility milestone.

**Decisions, confirmed with the repo owner via `AskUserQuestion` before
implementation**: (1) a new standalone `quaidsRobustFit()` sibling proc,
not a modification of `quaids.src` -- matching this project's dominant
precedent since Milestone 16 (curvature, welfare, shares, zero-
correction), over Milestone 12's narrower "modify shipped code"
exception; accepted trade-off: a **simplified bread** (`inv(gg)`-based,
not `quaidsFit()`'s own nonlinear-price-index-feedback Jacobian
correction), the same class of honestly-documented simplification
`quaidsCurvatureFit()`/`quaidsZeroFit()` already ship. (2) Robust and
cluster-robust **unified through one `clusterId` argument** rather than
two separate code paths: `clusterId = 0` means heteroskedasticity-robust
(every observation its own cluster, `G = nobs`); a supplied `Tx1` group-
label vector means cluster-robust with a CR1 small-sample correction --
marginal extra complexity over robust-only, since `clusterId = 0` is the
literal `G = nobs` special case of the same formula. (3) A **cluster-
aware bootstrap ships in this same milestone** -- the repo owner chose to
include this now rather than defer it as a separate follow-up, unlike the
Milestone 10/15 curvature/bootstrap split.

**The math**, given an already-fitted `qOut` and the raw sample:

1. Rebuild per-observation model-implied fitted shares at `qOut.bestB`
   using the same share formula already independently duplicated in
   `quaidsElas_()`/`quaidsSharesFit()`'s own `_quaidsSharesAt()` -- a
   *fourth* independent copy (`_quaidsRobustFittedW()`), vectorized across
   every sample row at once (the same `sumc()`-based idiom
   `tests/quaidsfixtures.src`'s own DGP construction and `quaidsFit()`'s/
   `quaidsZeroFit()`'s own iteration loops already use for whole-sample
   evaluation). Residuals `U = w[.,1:n1] - fittedW[.,1:n1]` (only the
   `n1` independently-estimated equations -- equation `n` is recovered
   via adding-up, never separately estimated, and gets no independent SE
   here either, matching every other diagnostic in this library).
2. Rebuild the shared regressor block `X = intcptFull~pricesHybrid[.,1:n1]
   ~endog~qOut.u` at the converged point (one evaluation, not a re-
   iteration). `pricesHybrid` mirrors `_quaidsIVFirstStage()`'s own
   mutation (columns `1:n-1` relative to the reference good `n`, column
   `n` left absolute) -- when `n1=n-1` (`aCtl.homogenous=1`) this drops
   the reference column entirely; when `n1=n` (`aCtl.homogenous=0`) it
   keeps all `n` hybrid columns, exactly mirroring `quaidsFit()`'s own
   STARTING VALUE block's `prices[.,1:n1]` regressor construction
   (`quaids.src:240`) for either case -- an honestly-documented
   simplification: this sandwich is built over the model's reduced-form
   moment-condition regressor set, not the exact post-minimum-distance-
   restriction regressor set `quaidsFit()`'s own homogeneity/symmetry
   stages additionally refine.
3. `Infl` (`nobs x (n1*k)`, `k=cols(X)`): for `i=1..n1`,
   `Infl[.,(i-1)*k+1:i*k] = X.*U[.,i]` -- the standard "per-observation
   score contribution" construction, via GAUSS's column-broadcast
   multiply (confirmed empirically that `X.*u` correctly broadcasts an
   `Tx1` vector against a `TxK` matrix column-wise, not assumed).
4. `clusterId=0`: `InflG=Infl`, `G=nobs`. Else: aggregate `Infl`'s rows by
   `clusterId` into `InflG` (`G x (n1*k)`) -- confirmed empirically (a
   small hand-verifiable example, not assumed from GAUSS's docs alone)
   that `un = unique(clusterId,1); D = clusterId .== un'; InflG = D'Infl;`
   correctly sums each group's rows via a single indicator-matrix
   multiply, no per-observation loop needed.
5. `Meat = c*(InflG'InflG)/nobs`, `c = (G/(G-1))*((nobs-1)/(nobs-K))`,
   `K=n1*k` -- the standard CR1 finite-sample correction, unifying the
   robust (`G=nobs`) and cluster (`G<nobs`) cases through the same
   formula.
6. `bread = eye(n1).*.inv(gg)`; `v = bread*Meat*bread`;
   `se = reshape(sqrt(diag(v)), n1, k)'` -- the exact reshape/transpose
   idiom Milestone 18 found broken twice already (Milestone 10/15). This
   milestone's own test writes the cell-position regression guard from
   day one, not after finding the bug a third time.

`rOut.b` is an **exact** (not approximate) re-expression of `qOut.bestB`'s
first `n1` columns in the `k`-row reduced-form basis `X`'s own columns
use: gamma's `(n1+1..n)`-th rows are dropped (redundant under exact
homogeneity/adding-up, and never a free regressor in `X`), giving `rOut.b`
exactly `k` rows with no approximation involved -- confirmed by direct
algebraic derivation (under homogeneity, `gamaFull[.,i]'*prices` collapses
exactly to `pricesRel[.,1:n1]*gama[1:n1,i]`, the same identity the
STARTING VALUE block's own regressor choice already relies on) before
writing the slicing code, not discovered by trial and error.

**Real bugs found and fixed while building this, all via direct empirical
testing** (this project's standing discipline):

1. **GAUSS identifiers are case-insensitive** -- a local named `K`
   silently collided with an already-declared local `k` in the same
   `local` statement (`error G0089: Duplicate definition of local 'K'`),
   and the subsequent `K = n1*k;` assignment then overwrote `k` itself,
   corrupting the `Infl[.,(i-1)*k+1:i*k]` column-range construction
   downstream (`error G0046: Columns don't match`). Fixed by renaming to
   `Kdof`, distinct in more than case from `k`.
2. **A genuine shape mismatch between the bootstrap and the closed-form
   sandwich.** An early version of `quaidsRobustBootstrapFit()` tracked
   `vec(qOut.bestB)` (the full, adding-up-recovered, `n`-column shape) as
   its bootstrap draws, while `quaidsRobustFit()`'s own `se` is in the
   `n1`-column reduced form -- one fewer row whenever `qOut.bestB`'s
   gamma block's redundant `n`-th row is dropped. Invisible until
   `printQuaidsRobustBootstrap()` was actually run against real data
   (`error G0058: Index out of range`, not a silently-wrong number).
   Fixed by extracting a shared private helper, `_quaidsRobustReduceB()`,
   called identically by both `quaidsRobustFit()` and the bootstrap's own
   `_quaidsRobustBootOneRep()`, so both report the identically-shaped
   reduced form rather than two independently-drifting copies of the same
   slicing logic.
3. **A genuinely flaky test, found only by running it repeatedly, not by
   reading it.** `tests/quaids_robust_bootstrap_test.e`'s core "cluster
   seBoot > naive seBoot" check passed on the first full-suite run but
   failed on the very next release-verification run with the identical
   code and seed -- confirmed, by running the file three times in a row
   directly (`tgauss -b -x`), to be genuinely nondeterministic (pass,
   fail, fail), not a fluke of one run. Root cause, isolated by first
   confirming `rndseed`+`rndu` themselves ARE fully reproducible in a
   standalone probe: `tests/quaidsfixtures.src`'s new
   `_quaidsClusterSyntheticDGP()` drew `clusterId` via bare
   `ceil(rndu(tobs,1)*nClusters)` -- unlike every other draw in that
   fixture (and every other fixture in this file), which all use
   `rndns(rows,cols,seed)`'s explicit-seed form. `rndu()` has no such
   seeded-call form (confirmed empirically: `rndu(r,c,seed)` throws
   `error G0136`), so this one draw depended on GAUSS's ambient global
   random state -- meaning *which rows shared a within-cluster shock*
   changed from run to run even though every structural DGP parameter was
   otherwise fixed. Fixed with one `rndseed seed;` call at the top of the
   proc, restoring full reproducibility (confirmed by re-running the
   bootstrap test four times in a row afterward, all passing identically).
4. **`tests/package_public_api.e` reused the wrong fixture.** The new
   `quaidsRobustBootstrapFit()` block initially reused this file's main
   seed=11 dataset (already documented in this same file as "a known
   non-converging seed for the iterated estimator") -- but the robust APIs
   require a converged base fit: `quaidsRobustFit()` requires the caller's
   `qOut` to have converged, and `quaidsRobustBootstrapFit()` explicitly
   requires its own internal base `quaidsFit()` call to converge (mirroring
   `quaidsCurvatureBootstrapFit()`'s identical requirement), throwing
   `"the base (unresampled) quaidsFit() did not converge"` against the
   real installed package -- caught by the release-verification gate
   itself, not by re-reading the file. Fixed by reusing the file's own
   already-converging seed=500 AIDS fixture
   (`wC`/`pricesC`/`totexpC`/`instrC`/`aCtlC`, the same one
   `quaidsCurvatureFit`'s own block already uses) instead.

**A real, empirically-confirmed finding about the simplified bread's
practical consequence, not a theoretical worry**: comparing
`quaidsRobustFit()`'s closed-form `se` (with `clusterId=0`) against
`qOut.homogSE` on a real fitted dataset showed the closed-form SE was
**dramatically more conservative** -- often more than an order of
magnitude larger, especially for the IV-residual coefficient (whose own
column has deliberately low marginal variance by this fixture family's
own strong-instrument design, `totexp = .85*instr + u` with
`u = .1*rndns(...)`). Rather than assume this was a bug, it was checked
directly: an independent hand-derivation of a classical (non-robust)
`S.*.gg`-style formula, built from the *same* `X`/`U` `quaidsRobustFit()`
itself uses (not `qOut`'s own, more efficient, Jacobian-corrected
machinery), landed within roughly a factor of 1-3 of the robust `se` --
confirming the huge gap against `qOut.homogSE` is entirely attributable
to comparing a simple equation-by-equation sandwich against the full
cross-equation-efficient FGLS system, not a defect in the sandwich
formula. `quaidsRobustBootstrapFit()`'s bootstrap SE, which resamples and
refits the actual efficient estimator, confirmed this: it landed much
closer to `qOut.homogSE` than the closed-form sandwich did on the same
data. This is documented prominently (this file, `docs/USAGE_GUIDE.md`,
`docs/METHODOLOGY_NOTES.md`, and `quaidsRobustFit()`'s own command-
reference page) as a real, expected property of the simplified-bread
design choice, not silently left for a user to discover.

**Testing**: `tests/quaids_robust_test.e` (26 checks after Milestone 22) -- the point
estimate matches a fresh, independent hand-evaluation of the sandwich
formula; the exact-identity regression guard (`clusterId=0` vs. an
explicit `seqa(1,1,nobs)` per-row label give byte-identical output); the
"stays in the same order of magnitude as an independently-derived
classical formula using the same `X`/`U`" check described above; the
reshape/cell-position regression guard; shape/finiteness/non-negativity;
and the core non-vacuous check -- on a new fixture,
`_quaidsClusterSyntheticDGP()` (`tests/quaidsfixtures.src`, a 5-good AIDS
dataset with a genuine within-cluster-correlated noise component, needed
because every other fixture in this file has plain i.i.d. noise and
cannot distinguish a correct cluster-robust SE from an incorrect one that
ignores clustering) -- cluster-robust `se` is measurably larger than the
naive `se`, the standard textbook consequence of ignoring real
clustering; and full-basis covariance expansion drives shares,
elasticities, and welfare robust delta-method SE without changing point
estimates. `tests/quaids_robust_bootstrap_test.e` (17 checks after
Milestone 22, added to
the existing `-SkipBootstrap`-gated group rather than a new flag) --
bootstrap bookkeeping, shape/finiteness, the reshape regression guard for
`seBoot`, exact echo of the base point estimate/`seRobust`, full-basis
bootstrap covariance expansion for downstream shares, and the
identical cluster-vs-naive non-vacuous check for the bootstrap variant.
`tests/package_public_api.e` gained calls to all four new procs
(`clusterId=0` only -- cluster-specific behavior is already thoroughly
validated in the two dedicated test files above; this is a release gate,
not a re-validation).

**Version bump to `0.15.0`**: two new required public procs
(`quaidsRobustFit`, `quaidsRobustBootstrapFit`, plus their paired
printers) and two new required public structs (`quaidsRobustOut`,
`quaidsRobustBootOut`), matching this project's established policy of
bumping on real new public API surface. No new package dependency --
`tsmt`'s `robustSE`/`clusterSE` were evaluated and not adopted (see
above), so this milestone needed no external package beyond GAUSS's own
base runtime primitives (`moment`, `invpd`, `unique`, `.*.`).

**This completes the originally-outlined five-item full-demand-system-
workflow**: predicted shares (Milestone 16), a quadratic-term
specification test (Milestone 17), bootstrap percentile confidence
intervals (Milestone 18), zero-share correction (Milestone 19), and now
robust/cluster standard errors (Milestone 20). One still-unrequested
follow-up remains, flagged at Milestone 19: homogeneity/symmetry
imposition on top of the Shonkwiler-Yen zero-share correction.

## Milestone 21: applied workflow driver (complete)

Post-20 development shifts from isolated post-estimation procs toward
one-call applied workflows that fit a system and return the most commonly
needed downstream quantities together. The first seed is
`quaidsWorkflowFit()` (`src/quaidsworkflow.src`), which deliberately
composes existing public APIs rather than adding another estimator:
`quaidsFit()` for the core fit, `quaidsSharesFit()` and `quaidsElasFit()`
at the sample mean, `quaidsRobustFit()` for robust/cluster-robust SE, and
numeric model/restriction summaries (`symPval`, `overidPvf`,
`quadraticPval` when the fit is an unconstrained QUAIDS comparison).
`quaidsWorkflowScenarioFit()` preserves the shorter base signature and adds
explicit welfare scenario inputs (`intcptPt`, `pricesPt0`, `pricesPt1`,
`totexpPt0`) to fill the same `quaidsWorkflowOut` struct's CV/EV fields
via `quaidsWelfareFit()`.

Current scope is intentionally conservative: mean-point post-estimation
and robust inference only when the base fit converges, and one explicit
welfare scenario at a time. Follow-ups tracked in `GOLD_STANDARD_TODO.md`
include export-ready result bundles/adapters, examples, and broader
workflow validation.

## Milestone 22: robust inference propagation (complete)

Milestone 20 correctly reported robust/cluster-robust coefficient SE in
`quaidsRobustFit()`'s reduced basis, but that basis was not directly
usable by `quaidsSharesFit()`, `quaidsElasFit()`, or `quaidsWelfareFit()`,
which all expect the full `qOut.bestB` layout. Milestone 22 adds two
public expansion helpers in `src/quaidsrobust.src`:
`quaidsRobustCovariance(qOut, rOut, aCtl)` and
`quaidsRobustBootstrapCovariance(qOut, rbOut, aCtl)`. They return
`qOut.bestB` plus a robust or bootstrap covariance transformed into the
full `vec(qOut.bestB)` basis by applying the same linear adding-up and
homogeneity recoveries used by the coefficient point estimates. The
helpers do not mutate `qOut`.

`quaidsWorkflowFit()` now uses `quaidsRobustCovariance()` internally and
fills robust/cluster-robust post-estimation fields:
`robustBestB`, `robustBestV`, `postRobustValid`, `sharesRobustSE`,
`sharesRobustV`, `incomeElasRobustSE`, `priceElasRobustSE`, and
`compPriceElasRobustSE`. `quaidsWorkflowScenarioFit()` propagates the
same robust covariance into welfare SE as `welfareRobustValid`,
`seCVRobust`, and `seEVRobust`.

Testing added source-tree parity in `tests/quaids_robust_test.e` and
`tests/quaids_workflow_test.e`, bootstrap expansion coverage in
`tests/quaids_robust_bootstrap_test.e`, and installed-package API calls in
`tests/package_public_api.e`. The important invariant is that robust
propagation changes standard errors/covariance only; predicted shares,
elasticities, and CV/EV point estimates stay identical to the classical
post-estimation calls.

## Milestone 23: preflight data/design diagnostics (complete)

Milestone 23 adds `src/quaidsdiagnostics.src` with
`quaidsPreflight()`/`printQuaidsPreflight()`. This is deliberately a
diagnostic/reporting layer, not a new estimator and not a hidden call to
`quaidsFit()`: it inspects the raw matrices before estimation and returns a
flat `quaidsPreflightOut` struct with machine-readable flags and metrics.

The first slice covers dimensions/control validity, finite-value checks
(`x .== x` plus `scalmiss()`), share adding-up, zero and negative budget
shares, price/total-expenditure/instrument variation, a first-stage design
invertibility guard, the same first-stage IV F statistic produced by
`_quaidsIVFirstStage()`, cluster counts/minimum cluster size/singletons,
and a simple convergence-risk screen (`0` low, `1` elevated, `2` high).
Hard failures set `ok=0`; zero shares, low variation, weak IV, few/singleton
clusters, and elevated convergence risk are warnings.

Testing: `tests/quaids_preflight_test.e` validates the clean path, zero-share
warning, negative-share hard failure, adding-up failure, explicit cluster
summaries, low price variation warning, and dimension-mismatch return. The
installed package smoke test also calls
`quaidsPreflight()` and `printQuaidsPreflight()` against the stable seed=500
fixture, asserting diagnostic shapes/design/IV/cluster fields rather than
requiring `pOut.ok == 1` on that noisy synthetic demand sample.

## Milestone 24: workflow preflight summary (complete)

Milestone 24 wires the preflight layer into the applied workflow path without
turning the workflow into a guard wrapper. `quaidsWorkflowFit()` now calls
`quaidsPreflight()` before `quaidsFit()` and echoes compact flat fields such
as `preflightOk`, `preflightWarnings`, `preflightIVFstat`,
`preflightDesignInvOk`, and `preflightNClusters` in `quaidsWorkflowOut`.
The preflight summary is diagnostic and non-gating: callers who want to stop
before estimation should still call `quaidsPreflight()` directly and branch
on `pOut.ok`.

Testing: `tests/quaids_workflow_test.e` now compares the workflow preflight
summary against a direct `quaidsPreflight()` call on the same fixture, and
`tests/package_public_api.e` asserts the installed-package workflow exposes
matching preflight fields. The full nested `quaidsPreflightOut` is
intentionally not embedded in `quaidsWorkflowOut`; the workflow struct stays
flat by convention.

## Milestone 25: sampling-weighted workflow evaluation (complete)

Milestone 25 starts the survey/microdata roadmap without altering the core
estimator. `quaidsSurveyWorkflowFit()` (`src/quaidssurvey.src`) calls the
existing `quaidsWorkflowFit()`, validates a `Tx1` nonnegative sampling-weight
vector, replaces the default sample-mean evaluation point with the weighted
mean of `qOut.intcptFull~prices~totexp`, and recomputes predicted shares,
elasticities, and robust propagated post-estimation SE at that point.

This is deliberately post-estimation support only. `quaidsFit()` still uses
its existing unweighted moment conditions and covariance formulas; there is
no claim of full survey-weighted estimation, replicate-weight variance,
strata handling, or design-based covariance. Those remain separate roadmap
items. **Update, Milestone 26**: the "estimator itself stays unweighted"
part of this scoping is resolved -- see "Milestone 26" below, which also
changes `quaidsSurveyWorkflowFit()`'s own behavior (its `weight` argument
now fits the estimator too, not just the evaluation point).

Testing (as of Milestone 25): `tests/quaids_survey_workflow_test.e` proved
that the weighted evaluation point matches a manual weighted mean, the
underlying fit was unchanged, direct `quaidsSharesFit()`/`quaidsElasFit()`
parity held at the weighted point, robust post-estimation SE were
recomputed there too, and constant weights reproduced the default workflow
up to numerical tolerance. `tests/package_public_api.e` also exercised the
installed public proc. See "Milestone 26" below for how this test file's
own assertions changed to match the new estimator-weighting behavior.

## Milestone 26: sampling-weighted estimation and design-based SE (complete)

Milestone 26 closes the biggest, most consequential piece of Milestone 25's
own "Follow-ups" note: making the *estimator itself* genuinely sampling-
weighted, with a matching weighted/clustered SE, not just a weighted
post-estimation evaluation point.

**Scope, confirmed with the repo owner via `AskUserQuestion` before
implementation**: (1) weighted point estimates plus a weighted/clustered
(Horvitz-Thompson-style linearized) sandwich SE only -- formal strata as a
concept distinct from clustering, and replicate-weight (BRR/jackknife)
variance, are explicitly deferred, matching this project's established
phased-scope precedent (Milestone 10's AIDS-only curvature before
Milestone 13's QUAIDS extension; Milestone 19's unconstrained-only zero-
correction). (2) Architecture: extend `quaidsFit()`/`_quaidsIVFirstStage()`
**in place** with a new optional `weight` argument, rather than
duplicating the ~700-line estimation core into a new sibling file --
the biggest deviation yet from the "new sibling proc, don't touch shipped
code" convention every milestone since 16 has followed, justified here
because weighting touches essentially the whole estimation core (starting
value, iteration loop, variance, overID test), not one separable phase.
Mirrors the one existing precedent for touching shipped code (Milestone
12's `aCtl.relax`): opt-in, default = current unweighted behavior, verified
byte-identical when omitted.

**The math**: standard survey-weighted least squares pre-multiplies both
sides of every cross-product by `sqrt(weight)` before the product, exact
via `(sqrt(w).*A)'(sqrt(w).*B) = A'diag(w)B`. Applied at every site in
`src/quaids.src`/`src/quaidsiv.src` where raw sample rows enter a
`moment()` call or cross-product: `_quaidsIVFirstStage()`'s `m1`/`sse`
moments, `quaidsFit()`'s STARTING VALUE moment, the in-loop coefficient
re-estimation, the chunked Jacobian-correction products, the `D` term, and
the overidentification block's `instr'w`/`endog'instr` terms. The
symmetry-restriction stage needed **no changes** -- confirmed a pure
function of already-computed `b`/`v`, so it inherits correct weighting
automatically. `weight` is renormalized internally so it sums to `nobs`
(the point estimate is invariant to any overall rescaling of `weight`),
keeping every existing `nobs`-denominated formula unaffected. New
`quaidsOut` fields `weighted`/`weightSum`/`effN` (Kish's effective sample
size, `(sum w)^2 / sum(w^2)`) report the weighting.

**A different, easy-to-conflate convention for the robust/cluster
sandwich, flagged prominently in code comments and given its own dedicated
test**: `quaidsRobustFit()`/`quaidsRobustBootstrapFit()` get the same
optional `weight`, but the standard pweight-robust-sandwich convention
scales the per-observation score contribution (`Infl`) by PLAIN `weight`,
not `sqrt(weight)` -- the Horvitz-Thompson convention (matching, e.g.,
Stata's `vce(robust)` + `pweight`). The bread (`gg`) keeps the
`sqrt(weight)` convention, matching the weighted design the point estimate
was fit under. `tests/quaids_survey_test.e` checks `quaidsRobustFit()`'s
output against two independent hand-evaluations -- one using the correct
convention, one deliberately using the wrong one -- and confirms the
implementation matches only the correct one, not just derives the formula
on paper.

**Real findings from building this, all confirmed empirically, not
assumed**:

1. GAUSS's dynargs/variadic mechanism (`proc (...) = name(fixedArgs, ...);`
   + `dynargsGet(n, default)`) works cleanly for a plain optional trailing
   matrix argument, not just structs -- confirmed via a scratch smoke test
   against `c:\gauss26\src\between.src`'s identical pattern before touching
   any real source file, giving `quaidsFit()`/`quaidsRobustFit()`/
   `quaidsRobustBootstrapFit()`/`quaidsWorkflowFit()`/
   `quaidsWorkflowScenarioFit()` the desired "zero blast radius when
   omitted" property with no fallback needed. `quaidsPreflight()` uses a
   REQUIRED positional `weight` instead, mirroring that proc's own
   established `clusterId` convention (a diagnostic pass, not an opt-in
   estimator extension) -- and because `_quaidsIVFirstStage()`, which
   `quaidsPreflight()` also calls directly, needed a required `weight` too
   (a private helper with few callers, all now updated:
   `quaidsFit()`, `quaidsPreflight()`, `quaidsZeroFit()`).
2. **Designing a genuinely non-vacuous "weighted beats naive" fixture took
   real iteration, not a first-try success.** The first design selected a
   biased sample based on an INCLUDED covariate (the demographic shifter
   `intcptPop`) -- and found, by actually fitting both naive and weighted
   `quaidsFit()` on the resulting sample, that naive estimation was barely
   biased at all, with no consistent max/mean-abs-diff improvement from
   weighting. This is the standard Heckman selection-bias result: selecting
   on a variable already included as a regressor does not bias OLS/GLS
   coefficients, since the model already conditions on it. The fix:
   `_quaidsSurveyWeightedDGP()` (`tests/quaidsfixtures.src`) selects on the
   idiosyncratic ERROR of good 1 instead (an unobserved-in-practice
   quantity, though known to this fixture generator, which is exactly what
   Horvitz-Thompson weighting needs to be given the true inclusion
   probability) -- a real self-selection story (e.g., households with
   unusually high, unmodeled good-1 consumption are more likely to respond
   to a survey). A second finding: a continuously-varying (logistic)
   inclusion probability produced a few very large individual weights
   (rare-selected rows), which inflated the weighted GLS solve's own
   variance enough to make the naive-vs-weighted comparison noisy and
   inconsistent; switching to a fixed TWO-STRATUM sampling rate (weights
   take only two values) fixed this. Seed screening (the same discipline
   as Milestone 3/10's own seed sensitivity work) found `seed=11` gives a
   dramatic, unambiguous example: excluding the IV-residual row (this
   project's own well-documented noisiest row), naive max/mean abs diff
   against the true population parameters is 1.58/0.23, vs. 0.21/0.03
   weighted -- roughly a 7x improvement on both metrics.
3. **`quaidsSurveyWorkflowFit()`'s own existing test asserted the opposite
   of what this milestone requires.** `tests/quaids_survey_workflow_test.e`
   (Milestone 25) explicitly checked `wfSurvey.bestB == wfBase.bestB`
   ("survey workflow leaves the underlying estimator unchanged") --
   correct under Milestone 25's scoping, now WRONG under Milestone 26's,
   since `quaidsSurveyWorkflowFit()`'s `weight` argument now also fits the
   estimator. Running the **full** existing test suite (18 files, no
   `-SkipBootstrap`) immediately after the core `quaids.src`/
   `quaidsiv.src`/`quaidsdiagnostics.src` changes -- the plan's explicit
   regression gate, before writing any new weighted-specific tests --
   confirmed this was the *only* file that failed; every other file's
   unweighted behavior stayed byte-identical. Fixed by replacing that
   check with its logical opposite (`wfSurvey.bestB` now matches a direct
   weighted `quaidsWorkflowFit()`/`quaidsFit()` call exactly, and a
   genuinely unequal weight measurably changes it relative to the
   unweighted baseline) -- a deliberate, documented behavior change to an
   already-shipped proc, not a bug, consistent with this milestone's own
   "biggest deviation yet from the don't-touch-shipped-code convention"
   framing.

**Version bump to `0.21.0`**: no new `.src` file and no new package
dependency, but real new required public API surface on six already-shipped
procs (`quaidsFit`, `_quaidsIVFirstStage`, `quaidsRobustFit`,
`quaidsRobustBootstrapFit`, `quaidsPreflight`, `quaidsWorkflowFit`,
`quaidsWorkflowScenarioFit`) plus three new `quaidsOut` fields and six new
`quaidsWorkflowOut`/`quaidsPreflightOut` fields, matching this project's
established policy of bumping on real new public API surface regardless of
whether a new file was added.

**Testing**: `tests/quaids_survey_test.e` (new, 19 checks -- see "Testing
status" below); `tests/quaids_workflow_test.e` (27 -> 32 checks);
`tests/quaids_preflight_test.e` (13 -> 16 checks);
`tests/quaids_survey_workflow_test.e` (15 -> 17 checks, with the behavior-
change fix above); `tests/package_public_api.e` gained real weight
exercises (an explicit uniform weight, reproducing the unweighted fit
exactly) for `quaidsFit`/`quaidsPreflight`/`quaidsRobustFit`/
`quaidsRobustBootstrapFit`/`quaidsWorkflowFit` against the real installed
package. Full local suite (18 files, no skips) re-ran clean after every
change.

## What GAUSS already provides — do not duplicate

Full detail and evaluation status is in `GOLD_STANDARD_TODO.md` under "What
GAUSS Already Provides." Summary:

- **No built-in SUR/systems-of-equations estimator and no existing
  AIDS/QUAIDS/demand-system implementation** anywhere in the GAUSS runtime or
  shipped packages (`pkgs/`). The cross-equation-restricted iterated FGLS
  core here is genuinely new — this is the library's reason to exist.
- **`gmmFitIV`** (`gmm.sdf`/`gmm_est.src`/`gmm_hac.src`) — evaluated at
  Milestone 2, **decision: not adopted**. It's a single-equation IV-GMM
  estimator (`y` is `N x 1`); under its `"onestep"`/`"unadj"` settings it's
  mathematically identical to the classical 2SLS the hand-rolled first stage
  already computes, so adopting it wouldn't change any estimates. It also
  doesn't expose `zzi` (`inv(Z'Z)`) or the full `[Z, endog]` moment matrix in
  the layout `quaidsFit()`'s later system-wide covariance and
  overidentification-test formulas need (`D*zzi*D'`, `z1z1i`, `zz1`, ...) —
  adopting it would mean unpacking a `gmmOut` struct to reconstruct those
  building blocks anyway, adding a package dependency for zero net
  simplification. The extraction that *did* happen (`quaidsiv.src`,
  Milestone 2) keeps the hand-rolled computation but gives it its own file.
- **`pubtable`** (`pkgs/pubtable`) is a complete LaTeX/HTML/RTF/CSV/XLS/
  Markdown table-export engine, already used by `gauss-qardl` via a small
  adapter bundled inside pubtable itself (`pubtable_qardl.src`, guarded by
  `#ifDef QARDL_SDF_INCLUDED`). **Adopted at Milestone 6**: this repo's own
  `src/pubtable_quaids.src` follows the same pattern (guard, `ptModelFromX`/
  `ptFromX` shape) but lives in this repo rather than physically inside the
  installed `pubtable` package — see "Milestone 6: reporting via pubtable"
  above for why.
- **`optmt`** (`pkgs/optmt`) is GAUSS's general-purpose nonlinear
  optimization package (PV-struct parameter packing, `pvPacki`/`pvUnpack`,
  minimizes a scalar objective). Flagged as the right tool for curvature
  imposition since Milestone 0's original scoping; **adopted at
  Milestone 10** for `quaidsCurvatureFit()`'s profiled nonlinear IV step
  (searching only `vech(A)`, with the remaining coefficients exactly
  identified by OLS at each candidate `A`) — this library's first real
  external package dependency (`package.json`'s `deps` array).
- **`loadd()`/`asdf()`/dataframe column selection** — used at Milestone 2 for
  `quaidsFull()` (see below). **Not** a `"y ~ x1 + x2"` formula string:
  AIDS/QUAIDS is a multi-equation system (N shares against N parallel
  prices), which doesn't fit GAUSS's single-equation formula grammar. Column
  name lists (`data[., stringArrayOfNames]`), matched positionally, are the
  natural fit instead.
- **`glm()`** (`c:\gauss26\src\glm.src`, GAUSS's own base-runtime
  generalized-linear-models proc, not a separate package) — **adopted at
  Milestone 19** for `quaidsZeroFit()`'s per-good first-stage probit
  (`glm(y, x, "binomial", ctl)` with `ctl.link="probit"`), no new
  `package.json` dependency. Confirmed empirically (not assumed) that it
  can hard-crash on some degenerate inputs in a way not caught by this
  codebase's usual `trap`/`scalmiss` guard — see "Milestone 19" above.
- **`robustSE`/`clusterSE`/`hacSE`** (`c:\gauss26\src\robust.src`, GAUSS
  base runtime) and the separately-licensed `tsmt` package's near-
  identical procs — evaluated at Milestone 20, **decision: not adopted**.
  Both are single-equation `(X'X)^-1 (...) (X'X)^-1` sandwiches for one
  dependent variable and shared regressors; neither generalizes to this
  library's stacked multi-equation system (`n1` share equations sharing
  one design matrix `X`). Adopting either would mean unpacking into
  exactly the same per-cluster score-aggregation math `quaidsRobustFit()`
  builds directly anyway, for zero net simplification — the same
  reasoning as the `gmmFitIV` evaluation above. `quaidsRobustFit()`
  reuses `unique()` (GAUSS base runtime) for cluster-label grouping
  instead, no new package dependency.
- **Primitives already used correctly and worth keeping**: `moment()`,
  `invpd()`, `solpd()`, `design()`/`vech()`/`xpnd()` (symmetric-matrix
  restriction algebra), `eigh()` (Slutzky check), `cdfchic`/`cdftc`/`cdffc`/
  `cdfnc`, `printfm()`, `quantile()`, `pdfn()` (probit density, Milestone
  19), `unique()` (cluster-label grouping, Milestone 20).

## GAUSS language conventions observed in this codebase

- **Variable naming**: short lowercase names (`w`, `n`, `nz`, `ge`, `gg`,
  `b`, `u`), matching the original author's terse econometrics-code style —
  preserve it rather than renaming to verbose identifiers.
- **Locals**: all local variables declared in one `local` statement at the
  top of each proc.
- **Struct declarations inside procs**: `struct quaidsControl aCtl;` appears
  as a formal parameter, not in the `local` list — standard GAUSS syntax.
  `struct quaidsOut qOut;` (a local, not a parameter) is declared the same
  way, separately from the `local` statement, inside `quaidsFit()`/`quaids()`.
- **Character-matrix name vectors vs. the native `string array` type**: name
  vectors built with the classic `0$+"X"$+ftocv(...)` idiom (`xnam`, `wnam`,
  `znam`, `unam`, `enam`) are legacy character matrices, not the newer
  `string array` type — struct fields holding them must be declared `matrix`,
  not `string array`, or you get a `G0071 Type mismatch` at assignment. Hit
  and fixed this exact error while building `quaidsOut` at Milestone 1. The
  *reverse* direction (legacy char-matrix into a native `string`/
  `string array`-typed field, e.g. a `pubtable` `ptModel`/`ptTable`) throws
  the same error and has no dedicated conversion builtin — see "Milestone
  6: reporting via pubtable" above for the `$|`-concatenation idiom that
  does work.
- **Loop style**: `do while ok; ... endo;` / `do while i<=n; ... i=i+1;
  endo;` — not GAUSS's newer `for` loop syntax.
- **Symmetric-restriction idiom**: `design(vec(xpnd(seqa(1,1,k*(k+1)/2))))`
  builds the selection matrix `R` such that `vec(G) = R*vech(G)` for
  symmetric `G` — this is how homogeneity/symmetry constraints get imposed
  via minimum distance. Reuse this idiom rather than re-deriving it.
- **Matrix concatenation**: `~` horizontal, `|` vertical, as usual in GAUSS.
- **Relative vs. absolute prices**: `prices` is converted to relative
  (`prices[.,1:n-1] - prices[.,n]`) near the top of `quaidsFit()`, then
  converted back to absolute before the recovery step at the end — mutating
  the input matrix in place *within `quaidsFit()`'s own local scope*. Because
  GAUSS passes matrix arguments by value, this mutation does **not** leak
  back into `quaids()`'s own `prices` local — the legacy wrapper's `prices`
  stays in its original absolute-price form throughout, which is exactly
  what lets `quaids()` reuse its own `prices`/`totexp`/`w`/`instr` arguments
  unchanged for the elasticities/descriptive-stats/Slutzky calls after
  `quaidsFit()` returns. `intcpt` is different: `quaidsFit()`'s internal
  constant-column-prepend is a real structural change (not a round-trip), so
  `quaidsFit()` returns the mutated version as `qOut.intcptFull` for
  `quaids()` to reuse.

## Testing status

Seventeen non-bootstrap-gated routine source-tree tests exist (plus a
`tests/guard_error_cases/` directory of six standalone expected-error
scripts, see below), all run from `tests/` as the working directory:

```
tgauss -b -x quaids_schema_test.e
tgauss -b -x quaids_formula_parity_test.e
tgauss -b -x quaids_synthetic_validation_test.e
tgauss -b -x quaids_published_validation_test.e
tgauss -b -x quaids_hypothesis_tests_test.e
tgauss -b -x quaids_elasticities_test.e
tgauss -b -x quaids_shares_test.e
tgauss -b -x quaids_pubtable_test.e
tgauss -b -x quaids_curvature_test.e
tgauss -b -x quaids_welfare_test.e
tgauss -b -x quaids_reliability_regression_test.e
tgauss -b -x quaids_zero_test.e
tgauss -b -x quaids_robust_test.e
tgauss -b -x quaids_preflight_test.e
tgauss -b -x quaids_workflow_test.e
tgauss -b -x quaids_survey_workflow_test.e
tgauss -b -x quaids_survey_test.e
```

`quaids_curvature_bootstrap_test.e` and `quaids_robust_bootstrap_test.e`
are two more, both gated behind `-SkipBootstrap` (see below) -- listed
with the other bootstrap-cost caveats further down, not in the block
above. `tests/guard_error_cases/` (six scripts: `robust_nonconverged_qout.e`,
`robust_bad_cluster_length.e`, `robust_one_cluster.e`,
`curvature_nonconverged_qout.e`, `curvature_invalid_sym.e`,
`quaids_bad_b0_shape.e`) each confirm one specific validation guard fails
loudly and clearly on bad input, added alongside the "Unreleased"
fail-fast fixes to `quaidsRobustFit()`/`quaidsRobustBootstrapFit()`/
`quaidsCurvatureFit()`/`quaidsFit()`'s `aCtl.b0` handling. Twenty-five
files run in total with no flags skipped (17 + 6 + 2 bootstrap-gated).

- `quaids_schema_test.e` (Milestone 1, 36 checks): asserts `quaidsOut` field
  values/shapes, that `quaidsFit()` prints nothing between call and return,
  and that the legacy `quaids()` wrapper's four returned matrices are
  exactly (not approximately) equal to the `quaidsOut` fields they're drawn
  from. Grew by 2 checks as part of the Milestone 21-25 handoff's `aCtl.b0`
  and degrees-of-freedom fixes (see the `## Unreleased`-turned-`0.20.0`
  entry in `CHANGELOG.md`).
- `quaids_formula_parity_test.e` (Milestone 2, 17 checks): asserts
  `quaidsFull(dataframe, ...)` produces numerically identical output to
  `quaidsFit(matrices...)` on the same underlying data, including the
  `extraVars == 0` path.
- `quaids_synthetic_validation_test.e` (Milestone 3, 22 checks): asserts
  `quaidsFit()` recovers true DGP parameters within a documented tolerance
  across LA-AIDS/iterated-AIDS/QUAIDS x with/without-IV. See "Milestone 3:
  synthetic validation findings" above for the tolerances and the seed
  sensitivity finding this surfaced.
- `quaids_published_validation_test.e` (Milestone 3, 19 checks): asserts
  `quaidsFit()` on real published data (`Blanciforti86`) matches an
  independent R reference within tolerance, plus adding-up/homogeneity/
  symmetry sanity checks. This is the test that caught the Stone-index
  starting-value bug — see "Milestone 3: real bug found and fixed" above.
  Milestone 9 extended it (11 -> 19 checks) to also validate iterated
  AIDS against R `aidsEst(method="IL", ...)` — see "Milestone 9" above.
- `quaids_hypothesis_tests_test.e` (Milestone 4, 22 checks): size and power
  for `quaidsHomogeneityTest()`/`quaidsJointTest()`, a power check for the
  existing symmetry-given-homogeneity test, and the first-ever exercise of
  the overidentification test. See "Milestone 4: new hypothesis tests"
  above. Milestone 17 added 3 checks for `quaidsQuadraticTest()` (size and
  power) -- see "Milestone 17: AIDS-vs-QUAIDS specification test" above.
- `quaids_elasticities_test.e` (Milestone 5, 17 checks): parity between
  `quaidsElasFit()`/`printQuaidsElas()` and `quaidsElas_()`, plus three
  exact algebraic identities checked at points other than the four
  standard ones. See "Milestone 5: elasticities generalization" above.
  Milestone 16 replaced its private `modelShareAt()` helper with a direct
  call to the new `quaidsSharesFit()`.
- `quaids_shares_test.e` (Milestone 16, 21 checks): `quaidsSharesFit()`'s
  point estimate matches a fresh, independently hand-evaluated share
  formula on both an AIDS and QUAIDS fixture; exact adding-up
  (`sum(w)==1`); shape/finiteness/non-negativity of `se`/`v` and
  `se==sqrt(diag(v))`; non-vacuousness at a shifted point. See "Milestone
  16: predicted budget shares" above.
- `quaids_pubtable_test.e` (Milestone 6, 37 checks; requires the `pubtable`
  package installed): exact numeric parity between `pubtable`
  `ptModel.estimates`/`stdErrors` and the `qOut`/`elasOut` values they're
  built from, shape/title checks, the `ptFromQuaidsFamily` dispatcher, and
  an end-to-end export smoke test that writes real `.tex`/`.md`/`.csv`
  files and reads them back. See "Milestone 6: reporting via pubtable"
  above. Milestone 21 added checks for `ptTablesFromQuaidsWorkflow()`.
- `quaids_curvature_test.e` (Milestone 10, 35 checks total; requires the
  `optmt` package installed): AIDS block (19 checks) -- recovery against
  a known-curvature-consistent true gamma, exact adding-up/homogeneity/
  symmetry, near-exact negative-semidefiniteness at the reference point,
  and a non-vacuousness check. Milestone 13 added a QUAIDS block (14
  checks) -- convergence, exact NSD, non-vacuousness, and shape/
  finiteness against the general QUAIDS fixture, deliberately not
  true-gamma recovery (see "Milestone 13: QUAIDS curvature imposition"
  for why). See "Milestone 10: curvature imposition" above. Milestone 18
  added a `cOut.se` cell-position regression guard to both blocks (2
  checks) -- see "Milestone 18: percentile bootstrap confidence
  intervals" above.
- `quaids_welfare_test.e` (Milestone 11, 20 checks): exact zero-price-
  change identity, an exact round-trip identity (feeding the computed
  expenditure function's output back into the indirect utility function
  returns the original utility level), a first-order/Marshallian limiting-
  case numerical check, SE finiteness/non-negativity, and a CV/EV sign-
  agreement check — run once against a QUAIDS fit and once against an
  AIDS fit. See "Milestone 11: welfare measures" above.
- `quaids_reliability_regression_test.e` (Milestone 12, 11 checks):
  `aCtl.relax=1` reproduces byte-identical output to leaving `aCtl.relax`
  unset (the new field is a true no-op at its default); a previously-
  crashing seed (QUAIDS, seed 43 — an unguarded `invpd()` in the
  symmetry-test block) now degrades gracefully (`qOut.converged==0`,
  `qOut.symValid==0`) instead of aborting the call; a concrete example
  (QUAIDS, seed 2) confirms `aCtl.relax=.75` measurably changes a
  never-converged fit into a correctly-converged one. See "Milestone 12:
  numerical reliability" above. Grew by 3 checks as part of the Milestone
  21-25 handoff's `aCtl.b0`/degrees-of-freedom fixes.

- `quaids_curvature_bootstrap_test.e` (Milestone 15, 37 checks; requires
  the `optmt` package installed): bootstrap run bookkeeping (requested/
  completed/failed/attempts), shape/finiteness of the bootstrap SE, exact
  echo of the base (unresampled) point estimate and delta-method SE, and
  a plausibility check that the bootstrap SE stays well-behaved where the
  delta-method SE does not — on both an AIDS and a QUAIDS fixture, with a
  small `B` (15/5 respectively) chosen to bound this file's own runtime.
  **Not** run by `run_source_tests.ps1`'s default invocation — see
  `-SkipBootstrap` and "Milestone 15: bootstrap standard errors" above.
  Milestone 18 added `quaidsCurvatureBootstrapCI()` checks and a
  `seBoot` cell-position regression guard (11 checks) — see "Milestone
  18: percentile bootstrap confidence intervals" above.
- `quaids_zero_test.e` (Milestone 19, 19 checks): the fixture's own
  adding-up identity (exact, by construction); the diagonal-delta
  restriction holds exactly (off-diagonal hazard-coefficient entries
  exactly `0`, on-diagonal entries genuinely estimated); shape/finiteness
  of `probitB`/`se`/`b`; all `n` first-stage probits converged; and the
  core validation — `quaidsZeroFit()`'s corrected coefficients recover
  the true latent (uncensored) DGP parameters better than naively fitting
  `quaidsFit()` on the same censored data, on both a max- and mean-
  absolute-difference basis. See "Milestone 19: zero budget share
  correction (Shonkwiler-Yen)" above. Grew by 2 checks validating
  Milestone 21-25's `aCtl.b0` support for `quaidsZeroFit()`.
- `quaids_robust_test.e` (Milestone 20/22, 26 checks): the point estimate
  matches a fresh, independent hand-evaluation of the sandwich formula;
  an exact-identity regression guard (`clusterId=0` vs. an explicit
  `seqa(1,1,nobs)` per-row label give byte-identical output); the sandwich
  stays in the same order of magnitude as an independently-derived
  classical formula built from the same regressors/residuals; the
  reshape/cell-position regression guard (written from day one); shape/
  finiteness/non-negativity; and the core non-vacuous check — cluster-
  robust `se` is measurably larger than the naive `se` on
  `_quaidsClusterSyntheticDGP`'s genuinely clustered data. See "Milestone
  20: robust and cluster-robust standard errors" above. Milestone 22 adds
  full-basis covariance expansion checks for robust shares, elasticities,
  and welfare.
- `quaids_robust_bootstrap_test.e` (Milestone 20/22, 17 checks): bootstrap
  run bookkeeping, shape/finiteness, the reshape regression guard for
  `seBoot`, exact echo of the base point estimate/`seRobust`, full-basis
  bootstrap covariance expansion for downstream shares, and a
  plausibility check that a cluster-aware bootstrap's `seBoot` exceeds a
  plain-row bootstrap's on the same genuinely clustered data. **Not** run
  by `run_source_tests.ps1`'s default invocation — added to the same
  `-SkipBootstrap`-gated group as `quaids_curvature_bootstrap_test.e`
  rather than a new flag.
- `quaids_preflight_test.e` (Milestone 23, 16 checks): clean preflight,
  zero-share warning, negative-share hard failure, adding-up failure,
  explicit cluster summaries, low price variation warning, and
  dimension-mismatch return. Milestone 26 added 3 checks: a genuinely
  unequal valid weight's `weightValid`/`weightSum`/reduced `effN`, and a
  non-finite weight entry's hard preflight error.
- `quaids_workflow_test.e` (Milestone 21 seed, 32 checks): checks that
  `quaidsWorkflowFit()` returns the same core fit, sample-mean evaluation
  point, predicted shares, elasticities, robust coefficient SE, and robust
  propagated shares/elasticity SE as
  explicit calls to `quaidsFit()`/`quaidsSharesFit()`/`quaidsElasFit()`/
  `quaidsRobustFit()`/`quaidsRobustCovariance()` on the same fixture. It
  also checks symmetry/overidentification/quadratic summary fields,
  Milestone 24's compact preflight summary fields against a direct
  `quaidsPreflight()` call, and verifies
  `quaidsWorkflowScenarioFit()`'s classical and robust CV/EV fields against
  `quaidsWelfareFit()`. Milestone 26 added 5 checks: `quaidsWorkflowFit()`'s
  own optional estimator `weight` argument matches direct
  `quaidsFit()`/`quaidsPreflight()`/`quaidsRobustFit()` calls with the same
  weight, and a genuinely unequal weight changes the workflow's point
  estimate (non-vacuous).
- `quaids_survey_workflow_test.e` (Milestone 25, 17 checks): checks that
  `quaidsSurveyWorkflowFit()` validates and echoes sampling-weight
  diagnostics, computes the same weighted evaluation point as a manual
  weighted mean, recomputes shares/elasticities and robust propagated SE
  at that point, and reproduces the default workflow under constant
  weights up to numerical tolerance. **Milestone 26 behavior change**: this
  file's own "leaves the underlying estimator unchanged" check (a
  Milestone 25 assertion) was replaced with the opposite claim, matching
  `quaidsSurveyWorkflowFit()`'s own new behavior -- `wfSurvey.bestB` now
  matches a direct weighted `quaidsFit()` call exactly, and a genuinely
  unequal weight measurably changes it relative to the unweighted
  baseline. Full-suite regression testing (`tests/run_source_tests.ps1`
  with no flags skipped, run immediately after the core `quaids.src`/
  `quaidsiv.src`/`quaidsdiagnostics.src` changes, before writing any new
  weighted-specific tests) confirmed this was the *only* existing test
  file whose assertions needed updating -- every other file's unweighted
  behavior stayed byte-identical.
- `quaids_survey_test.e` (Milestone 26, 19 checks): the exact-identity
  regression guard (weight omitted / an explicit uniform weight reproduce
  the pre-Milestone-26 unweighted behavior byte-for-byte across
  `quaidsFit()`/`quaidsPreflight()`/`quaidsRobustFit()`/
  `quaidsRobustBootstrapFit()`); an invalid weight's hard preflight error;
  the sqrt(weight)-bread vs. plain-weight-meat sandwich convention split,
  checked directly against two fresh hand-evaluations (one matching the
  documented convention, one deliberately using the wrong one -- see
  "Milestone 26" below for why this is the single easiest detail here to
  get backwards); and the core non-vacuous check -- weighted estimation
  recovers `_quaidsSurveyWeightedDGP`'s true population parameters
  measurably better (both max- and mean-abs-diff) than naive estimation on
  the same informatively-sampled data.

All eighteen routine source-tree files print one
`PASS`/`FAIL` line per check and a final `...: ALL N CHECKS PASSED` (or a
failure count) summary line — check that line, since `tgauss`'s exit code
is not currently a reliable pass/fail signal for this harness. `tests/
run_source_tests.ps1` (Milestone 7) runs `verify_package_manifest.ps1`
plus all of these in one shot (except `quaids_curvature_bootstrap_test.e`/
`quaids_robust_bootstrap_test.e` when `-SkipBootstrap` is passed, as the
automatic CI workflow does) and checks this same summary-line convention
(not just GAUSS-level compile/execute errors).

`tests/quaids_convergence_sweep.e` (Milestone 12) is a real, committed
diagnostic tool, but is deliberately **not** one of the ten above — it
prints no "ALL N CHECKS PASSED" line and is not run by
`run_source_tests.ps1` (there is no convergence guarantee to gate on).
Run it manually via `tests/run_convergence_sweep.ps1` whenever you want
to re-measure the iterated estimator's convergence-failure rate — see
"Milestone 12: numerical reliability" above.

`tests/package_public_api.e` (Milestone 7) is different in kind from the
source-tree files above: it loads `library quaids;` against a real
*installed* copy of the package (currently `c:\gauss26\pkgs\quaids`) rather
than `#include`-ing the source tree, so it only runs correctly after
`scripts/run_release_verification.ps1 -InstallArtifact` (or equivalent)
has installed the package — it is a release gate, not a routine
`run_source_tests.ps1` member. See "Milestone 7: package build and release
tooling" above.

`examples/quaids_example.e` remains a manual, eyeball-comparison smoke
script (no assertions) — superseded for correctness-checking purposes by
`quaids_synthetic_validation_test.e`, but kept as a simple, readable
end-to-end usage example. `examples/pubtable_export_example.e` (Milestone
6) is the same style for the reporting layer — exports a coefficient table
and elasticity tables to `.tex`/`.md`/`.csv`, requires `pubtable` installed.
Both use `library quaids;` against the installed package (Milestone 9),
matching README.md/USAGE_GUIDE.md's documented usage pattern — see
"Milestone 9" below.

To run the examples from a GAUSS 26 console/batch shell:

```
tgauss -b -x examples/quaids_example.e
tgauss -b -x examples/pubtable_export_example.e
```

`quaids_example.e` no longer depends on its working directory, since it
loads the installed package via `library quaids;` rather than `#include`.
`pubtable_export_example.e` still needs `examples/` as the working
directory (or paths adjusted), since it `#include`s
`../src/pubtable_quaids.src` by relative path — that adapter is not part
of the installed package (see "Milestone 6" above).

## Package manifest

`package.json` lists (relative to `src/`, in load order): `quaids.sdf`,
`quaidsutil.src`, `quaidsiv.src`, `quaidszerocorrect.src`,
`quaidselas.src`, `quaidsshares.src`, `quaidsslutzky.src`, `quaids.src`,
`quaidsformula.src`, `quaidstests.src`, `quaidscurvature.src`,
`quaidswelfare.src`, `quaidsrobust.src`, `quaidsdiagnostics.src`,
`quaidsworkflow.src`, `quaidssurvey.src`.
`quaidsshares.src` (Milestone 16) has no load-order dependency on
anything beyond `quaids.sdf` (its private `_quaidsSharesAt()` helper is
a fresh, independent implementation, not a call into `quaidselas.src`)
and adds no new entry to `deps` — pure closed-form algebra, same
footprint as `quaidswelfare.src`.
`quaidszerocorrect.src` (Milestone 19) loads right after `quaidsiv.src`,
which its `quaidsZeroFit()` calls directly (`_quaidsIVFirstStage()`) —
no other load-order dependency, and no new `deps` entry, since `glm()` is
part of GAUSS's own base runtime, not a separate package.
`quaids.src` must load after `quaidsiv.src`/`quaidselas.src`/
`quaidsslutzky.src` since it calls procs they define; `quaidsformula.src`
must load after `quaids.src` since `quaidsFull()` calls `quaidsFit()`.
`quaidstests.src` has no load-order dependency on the others beyond
`quaids.sdf` (it only reads an already-computed `quaidsOut`).
`quaidscurvature.src` (Milestone 10) loads after `quaids.src` (needs
`quaidsFit()`'s output) -- unlike `pubtable_quaids.src`, it IS required in
`src` (real public API), so `package.json`'s `deps` array now lists
`optmt` (this library's first real external package dependency).
`quaidswelfare.src` (Milestone 11) loads last; it has no load-order
dependency on anything beyond `quaids.sdf` (it only reads already-fitted
`b`/`v` coefficient/covariance arguments, the same footprint as
`quaidselas.src`) and adds no new entry to `deps` — pure closed-form
algebra, no external package needed. Milestone 12 added no new `src`
file — `quaids.sdf`/`quaidsutil.src`/`quaids.src` were modified in place
(the new `aCtl.relax` field, the crash-guard, and the denominator guard)
— but still bumped the version to `0.8.0`, since `relax` is new public
`quaidsControl` API surface regardless of which file it lives in.
Milestone 13 also added no new `src` file — `quaidscurvature.src` itself
was extended in place to accept QUAIDS fits, reusing `aCtl.relax` rather
than adding another field — but bumped the version to `0.9.0` anyway
(real new capability, not a bugfix; see "Milestone 13: QUAIDS curvature
imposition" above).
`quaidsrobust.src` (Milestone 20) loads last; it calls `quaidsFit()`
directly (inside `quaidsRobustBootstrapFit()`'s per-replication refit), so
it must load after `quaids.src`, but adds no new `deps` entry —
`unique()` (cluster-label grouping) is part of GAUSS's own base runtime,
not a separate package; `robustSE`/`clusterSE` (base runtime `robust.src`)
and `tsmt`'s near-identical procs were evaluated and not adopted (see
"What GAUSS already provides" above).
`quaidsdiagnostics.src` (Milestone 23) loads after the core public APIs
because `quaidsPreflight()` reuses `_quaidsIVFirstStage()` for first-stage
diagnostics and is itself called by `quaidsworkflow.src`.
`quaidsworkflow.src` (Milestone 21 seed) loads after diagnostics because it is a
composition layer over the existing public APIs: it calls `quaidsFit()`,
`quaidsSharesFit()`, `quaidsElasFit()`, `quaidsRobustFit()`, and, since
Milestone 24, `quaidsPreflight()` from `quaidsdiagnostics.src`; it adds no
new package dependency.
`quaidssurvey.src` (Milestone 25) loads after `quaidsworkflow.src` because
`quaidsSurveyWorkflowFit()` calls the base workflow and then recomputes
post-estimation shares/elasticities at a sampling-weighted evaluation
point. **Update, Milestone 26**: it no longer leaves the estimator
moments inside `quaidsFit()` unweighted -- its own `weight` argument is
now also forwarded into `quaidsWorkflowFit()`'s new optional `weight`
argument, fitting a genuine weighted estimator. Milestone 26 itself added
no new `src` file -- `quaids.sdf`/`quaidsiv.src`/`quaids.src`/
`quaidsrobust.src`/`quaidsdiagnostics.src`/`quaidsworkflow.src`/
`quaidssurvey.src` were all modified in place (a new optional `weight`
argument on `quaidsFit()`/`quaidsRobustFit()`/
`quaidsRobustBootstrapFit()`/`quaidsWorkflowFit()`/
`quaidsWorkflowScenarioFit()`, a required `weight` argument on
`quaidsPreflight()`/`_quaidsIVFirstStage()`, and new `quaidsOut`/
`quaidsPreflightOut`/`quaidsWorkflowOut` fields) — but still bumped the
version to `0.21.0`, matching the Milestone 12/13 precedent of bumping on
real new public API surface regardless of whether a new file was added.
No new `deps` entry — pure GAUSS built-ins, the same primitives every
other covariance/moment computation in this library already uses.
`src/pubtable_quaids.src` (Milestone 6) is deliberately **not** in this
array — it has a hard dependency on `pubtable.sdf`'s struct types, and
adding it would make `pubtable` a hard dependency for the whole package to
compile; see "Milestone 6: reporting via pubtable" above.
**Buildable/installable as of Milestone 7**: `scripts/
build_lcg.ps1`/`build_package.ps1`/`run_release_verification.ps1` build a
release `.zip` and can install it to a real GAUSS package directory
(`c:\gauss26\pkgs\quaids` on this machine) so `library quaids;` works — see
"Milestone 7: package build and release tooling" above. If you add a new
*required* `.src` file, add it to the `"src"` array (respecting load
order), bump the version, and rebuild/reinstall (`scripts/
run_release_verification.ps1 -BuildArtifact -ForceArtifact
-InstallArtifact`) so the installed copy and its `.lcg` catalog stay in
sync with the source tree; optional/adapter files like
`pubtable_quaids.src` still warrant a version bump (real new public API
surface) but stay out of `src`.

## References

- Deaton, A., Muellbauer, J. (1980). "An Almost Ideal Demand System."
  *American Economic Review*, 70(3), 312–326.
- Banks, J., Blundell, R., Lewbel, A. (1997). "Quadratic Engel Curves and
  Consumer Demand." *Review of Economics and Statistics*, 79(4), 527–539.
- Blanciforti, L., Green, R., King, G. (1986). *U.S. Consumer Behavior Over
  the Postwar Period: An Almost Ideal Demand System Analysis*. Giannini
  Foundation Monograph No. 40. Published-replication dataset (as
  `Blanciforti86` in the R package `micEconAids`), used by
  `tests/quaids_published_validation_test.e` — see
  `tests/fixtures/published/SOURCE.md`.
- Aptech GAUSS coding conventions: https://github.com/aptech/gauss-llm-reference
- Sibling library for structure/process conventions: `gauss-qardl`
  (`CLAUDE.md`, `GOLD_STANDARD_TODO.md`, `docs/command-reference/`,
  `pubtable_qardl.src` adapter pattern).
