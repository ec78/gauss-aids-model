# QUAIDS Public Release Roadmap

## Objective

Prepare the QUAIDS library for a trustworthy first public release with a
reproducible installation, an explicit stability contract, consistent
documentation, and a complete path from customer data to validated results.

This roadmap is intentionally narrower than `GOLD_STANDARD_TODO.md`. It tracks
only work that affects public-release readiness and customer implementation.

## Recommended Release Posture

Ship the first public version as a **public alpha** unless iterated AIDS and
QUAIDS wrong-solution detection is materially improved before release.

- Treat LA-AIDS as the stable baseline workflow.
- Label iterated AIDS, QUAIDS, curvature imposition, and their bootstrap paths
  according to their measured numerical limitations.
- Do not imply that `qOut.converged == 1` alone establishes a trustworthy
  iterated solution.
- Keep the first release API small and explicit. Resolve compatibility-sensitive
  naming and dependency decisions before promising API stability.

## Priority and Effort Scale

- **P0**: blocks the public release.
- **P1**: required for effective customer implementation; may proceed in
  parallel after P0 design decisions.
- **P2**: desirable release polish or post-release hardening.
- **S**: up to one day.
- **M**: roughly two to four days.
- **L**: roughly one to two weeks and likely iterative.

Effort estimates are relative engineering sizes, not calendar commitments.

## Phase 0: Define the Public Contract

**Status: complete** (PR-001, PR-002, PR-003) as of the `0.1.0 -
Unreleased` changelog entry. `scripts/verify_public_api.ps1` (run as
part of `tests/run_source_tests.ps1`/CI) enforces the release-metadata
version consistency and public-API-inventory acceptance evidence for
PR-001/PR-003 going forward -- it is not a one-time manual check. Not yet
done, deliberately, since it is a release action rather than a code
change: the actual Git tag, release title, and `CITATION.cff`
`date-released` are still pending until the repo owner is ready to cut
the release itself (PR-001's "byte-identifiable with the tagged commit"
evidence).

### PR-001 — Select the release identity and maturity level

- **Priority / effort:** P0 / S
- **Work:**
  - Choose the first public version. Recommended: `0.1.0` for a public alpha,
    or explicitly document why the internal milestone version `0.24.0` will be
    retained.
  - Apply the same version and release date to `package.json`, `CITATION.cff`,
    the changelog, artifact name, Git tag, and release title.
  - Move every change included in the artifact out of `Unreleased` and into the
    selected release entry.
  - State whether semantic-versioning compatibility begins with this release.
- **Acceptance evidence:**
  - One automated metadata check confirms the version across all release files.
  - The artifact is built from, and byte-identifiable with, the tagged commit.
  - The README, changelog, and citation metadata use the same maturity language.

### PR-002 — Declare the supported environment

- **Priority / effort:** P0 / S
- **Work:**
  - State the minimum GAUSS version and supported operating systems.
  - Separate "supported and tested" platforms from "expected to work" platforms.
  - Define whether command-line build/test scripts are Windows-only even if the
    installed GAUSS package is portable.
- **Acceptance evidence:**
  - The README contains a compact support matrix.
  - Every supported environment appears in the release test matrix or has a
    documented manual validation record.

### PR-003 — Freeze or repair compatibility-sensitive public names

- **Priority / effort:** P0 / M
- **Work:**
  - Decide how to handle the misspelled public field `homogenous` before users
    depend on it. Recommended: add a correctly spelled API path if GAUSS struct
    compatibility permits it, retain the old spelling as a deprecated alias for
    a defined period, and document precedence if both are set.
  - Decide whether `quaidsElas_` is a supported public procedure. Rename or hide
    it if it is internal; otherwise add it to the command index and compatibility
    policy.
  - Identify every public struct field and procedure covered by compatibility
    guarantees.
- **Acceptance evidence:**
  - A machine-readable public API inventory is committed and checked in CI.
  - The command reference covers every supported public procedure.
  - Deprecated names have tests and a documented removal policy.

## Phase 1: Make Installation Reproducible

**Status: PR-101 complete.** Implemented the "Recommended" option:
`quaidscurvature.src` is no longer in `package.json`'s `src` array or
`deps`; it is an opt-in adapter requiring `library optmt, quaids;
#include quaidscurvature.src`, mirroring the existing `pubtable` adapter
exactly. `docs/public-api.json`'s new `optional_modules` array and
`scripts/verify_public_api.ps1` enforce this going forward. See
CLAUDE.md's "Public Release Phase 1" section for the real GAUSS quirks
found while implementing and testing it.

**PR-102 partially complete.** Everything independently verifiable from
this environment is done (fresh install-directory delete/recreate every
release-verification run, `library quaids;` loading in a brand-new GAUSS
job, the installed-package gate -- now two separate files, one proving
core needs no `optmt` -- running against a freshly reinstalled copy, not
a stale development install). The real GAUSS Tools > Install Application
/ Package Manager acceptance evidence needs a clean machine and GUI
access this environment doesn't have; the repo owner will run that pass
themselves. Suggested checklist for that pass:
1. Install `quaids <version>.zip` via Tools > Install Application on a
   machine that has never had this package installed.
2. In a fresh GAUSS session, confirm `library quaids;` loads with no
   manual source-path configuration, then run the core quick start from
   the README.
3. Separately install `optmt`, then confirm
   `library optmt, quaids; #include quaidscurvature.src` plus a
   `quaidsCurvatureFit()` call works.
4. Install the immediately preceding internal artifact (`0.24.0`) first,
   then install `0.1.0` over it, and confirm the result matches a clean
   `0.1.0` install (no stale `.lcg` entries from the old catalog).
5. Uninstall, then reinstall `0.1.0` from scratch, and re-confirm step 2.
6. Record the exact GAUSS version, OS, and pass/fail outcome of each step
   for this release candidate.

### PR-101 — Resolve the `optmt` dependency model

- **Priority / effort:** P0 / M
- **Recommended decision:** Keep core estimation dependency-free by separating
  curvature support from the package's always-loaded source catalog. Treat
  curvature like the existing optional `pubtable` adapter, or publish a small
  companion package if GAUSS cannot conditionally compile the required structs.
- **Alternative:** If separation is not practical, declare `optmt` mandatory
  everywhere and load `library optmt, quaids;` in the quick start and all examples.
- **Work:**
  - Obtain an authoritative answer from Aptech or a verified public package about
    the `package.json.deps` schema and valid version constraint format.
  - Confirm the published `optmt` package name, availability, license terms, and
    current installable version.
  - Remove the current conflict between `package.json`, the README, the command
    reference, examples, and the installed-package test.
- **Acceptance evidence:**
  - A clean machine without QUAIDS or `optmt` can install and run the core quick
    start using the documented procedure.
  - A clean curvature installation succeeds using only documented steps.
  - The dependency behavior is tested through the official GAUSS installation
    path, not only through the repository's custom extraction/build scripts.

### PR-102 — Add clean-install acceptance tests

- **Priority / effort:** P0 / M
- **Work:**
  - Test **Tools > Install Application** against a clean GAUSS 26 installation.
  - Test the Package Manager path if the package will be distributed through a
    channel.
  - Test an upgrade over the immediately preceding internal artifact and a full
    uninstall/reinstall.
  - Record the exact environment, inputs, and output for each release candidate.
- **Acceptance evidence:**
  - Installation creates a usable catalog without manual source-path changes.
  - `library quaids;` or the deliberately documented alternative loads in a new
    GAUSS session.
  - The installed-package public API test passes against the newly installed
    artifact, not a pre-existing development installation.

### PR-103 — Make the release artifact self-contained

- **Priority / effort:** P0 / S
- **Work:**
  - Remove customer-document references to files omitted from the archive, or
    replace them with a concise shipped limitations/release-notes document.
  - Do not require customers to consult `CLAUDE.md` or
    `GOLD_STANDARD_TODO.md` for operational guidance.
  - Add an archive-level local-link check; the existing repository-level link
    check is not sufficient.
- **Acceptance evidence:**
  - Every relative link and every explicitly named guidance file in the shipped
    documentation exists inside the archive.
  - The installed package contains all files needed to run its shipped examples.

## Phase 2: Establish a Numerical Reliability Contract

**Status: PR-201 complete.** Added a "Model & Feature Support Tiers" table
directly in README.md, right after the Quick Start (not buried in a deep
reference doc): LA-AIDS is labeled Stable; iterated AIDS, QUAIDS,
zero-share correction, and curvature imposition are labeled Experimental
(QUAIDS marked highest risk, since it is `quaidsControlCreate()`'s actual
shipped default); bootstrap/replicate procedures are labeled as inheriting
their base model's tier. The same table, with the full "why" column, was
added as a new "Support Tier Summary" section at the top of
`docs/FEATURE_SUPPORT_MATRIX.md` (previously the convergence-rate data
existed only in a deep "## Notes" prose section). `docs/USAGE_GUIDE.md`'s
"Choosing A Model" table gained an explicit tier column, and
`docs/command-reference/quaidsControlCreate.md` now warns directly on the
`linear`/`maxiter` defaults that they select the highest-risk combination.
All three explicitly state what `qOut.converged == 1` does and does not
prove (tolerance convergence only, not solution uniqueness or recovery of
the intended fixed point) using the sweep's own "converged-but-wrong"
bucket as the concrete counterexample. The README's own Quick Start
example was also changed to lead with LA-AIDS (`aCtl.maxiter=1`) instead
of QUAIDS and now checks `qOut.converged` before proceeding -- this also
satisfies PR-203's second acceptance-evidence bullet ("the default quick
start does not silently enter the highest-risk estimator") without
touching `quaidsControlCreate()`'s actual coded defaults, per the repo
owner's explicit choice via `AskUserQuestion` (documentation-only fix;
PR-203's "change the coded default" option was declined). PR-202 (a
solution-stability diagnostic) was explicitly deferred per the same
decision -- it is not required for the First Public Alpha Exit Criteria,
which only lists PR-201 as required from this phase. Full source-tree
test suite re-ran clean after these docs-only edits (no `src/`/`tests/`
files touched).

### PR-201 — Define support tiers by estimator

- **Priority / effort:** P0 / S
- **Work:**
  - Publish separate support tiers for LA-AIDS, iterated AIDS, QUAIDS, zero-share
    correction, curvature imposition, and bootstrap procedures.
  - Put the measured convergence and wrong-fixed-point rates next to model
    selection guidance and the quick start, not only in deep reference notes.
  - Explain that `qOut.converged` detects tolerance convergence but does not prove
    solution uniqueness or recovery of the intended fixed point.
- **Acceptance evidence:**
  - A new user can identify the stable baseline model and the additional checks
    required for experimental estimators without reading internal roadmap files.

### PR-202 — Add a solution-stability diagnostic

- **Priority / effort:** P0 for a production-readiness claim; P1 for a clearly
  labeled public alpha / L
- **Work:**
  - Add a supported multi-start or restart workflow using documented starting
    values and damping settings.
  - Compare converged solutions using coefficients, objective/criterion values,
    fitted shares, and economic identities.
  - Return a stability flag distinct from `converged`, plus the number of distinct
    candidate solutions found.
  - Prevent high-level workflow helpers from presenting post-estimation results as
    routine output when stability has not been established.
  - Re-run the committed 200-seed sweep and preserve the results as a versioned
    release artifact.
- **Acceptance evidence:**
  - The diagnostic flags the previously documented wrong-fixed-point cases.
  - Tests cover non-convergence, a unique stable solution, and multiple distinct
    converged solutions.
  - Documentation defines exactly what each convergence/stability flag proves.

### PR-203 — Choose defaults based on evidence

- **Priority / effort:** P1 / M
- **Work:**
  - Evaluate whether model defaults should remain QUAIDS with `relax=1`, change to
    a safer baseline, or require explicit model selection.
  - Do not treat `relax=.75` as a complete fix; its measured benefit is modest.
  - If defaults change, document the compatibility impact and validate them across
    the full benchmark set.
- **Acceptance evidence:**
  - Default choices are justified by committed benchmark results.
  - The default quick start does not silently enter the highest-risk estimator.

## Phase 3: Correct and Simplify Customer Documentation

**Status: PR-301 complete.** A search-based review found three real,
confirmed contradictions (not the full list of Work bullets -- the
`optmt` install/load instructions and the `aCtl.b0`/`zOut.bRaw` warm-start
description were already reconciled by Phase 1/Phase 0's own edits):
`docs/command-reference/quaidsZeroFit.md` claimed `aCtl.homogenous = 0` is
the default (the real coded default, set in `src/quaidsutil.src`, is
`1`); `docs/USAGE_GUIDE.md`'s Zero Budget Shares section claimed
homogeneity/symmetry imposition is "out of scope in this first pass" for
`quaidsZeroFit()` even though Milestone 30 added it (its own code example
two paragraphs above already correctly showed `aCtl.homogenous = 1`);
and the same file's survey-workflow section claimed replicate-weight
(BRR/jackknife) variance "remain[s] roadmap items" even though
`quaidsReplicateWeightFit()` (Milestone 27) already ships it and this
same file documents it in its own dedicated section further down. All
three fixed. Added `scripts/verify_docs_consistency.ps1` (wired into
`tests/run_source_tests.ps1` and therefore CI) as the "targeted
documentation-consistency test" the acceptance evidence calls for: it
cross-checks `docs/command-reference/quaidsControlCreate.md`'s defaults
table against `quaidsControlCreate()`'s actual coded defaults for every
field, checks every doc page for any other stale `aCtl.homogenous`
default claim (not just the one file that happened to be wrong), and
carries regression guards for the two prose-only contradictions above.
Verified each guard actually fails by deliberately reintroducing each bug
in turn and re-running the script (one genuine bug in the check itself
was found and fixed this way -- a line-wrapped markdown phrase broke a
literal-space regex match) before confirming the final, corrected state
passes clean. Full source-tree suite re-ran clean afterward.

### PR-301 — Reconcile known contradictions

- **Priority / effort:** P0 / S
- **Work:**
  - Correct the stale statement that replicate-weight variance remains a roadmap
    item.
  - Correct the stale statement that zero-share homogeneity/symmetry is out of
    scope.
  - Correct the documented `quaidsZeroFit` default for `aCtl.homogenous`.
  - Correct `quaidsControlCreate` so zero-model warm starts require `zOut.bRaw`,
    not `zOut.b`.
  - Reconcile all `optmt` installation and load instructions.
- **Acceptance evidence:**
  - A targeted documentation-consistency test covers defaults, dependency status,
    and feature availability.
  - Search-based review finds no conflicting current-state claims.

### PR-302 — Separate customer guidance from engineering history

- **Priority / effort:** P1 / M
- **Work:**
  - Remove milestone chronology, bug archaeology, and implementation diary text
    from primary customer pages.
  - Preserve that material in the changelog, contributor documentation, or an
    architecture/history appendix.
  - Keep customer pages organized around decisions, inputs, outputs, limitations,
    and actions.
- **Acceptance evidence:**
  - The README supports installation and a first result without referencing an
    internal milestone.
  - Each limitations section states current behavior directly.

### PR-303 — Add documentation quality gates

- **Priority / effort:** P1 / M
- **Work:**
  - Check repository and archive links, heading structure, required command-page
    sections, and public-procedure coverage.
  - Validate code snippets where practical, especially signatures with keyword
    arguments.
  - Check release metadata and documented defaults against source definitions.
- **Acceptance evidence:**
  - Documentation checks run in CI and fail on a deliberately introduced stale
    default, missing command page, or archive-only broken link.

## Phase 4: Build the Customer Implementation Path

### PR-401 — Replace the quick start with the recommended workflow

- **Priority / effort:** P1 / S
- **Work:**
  - Lead with the shortest supported path from validated inputs to fitted results.
  - Show preflight gating, estimation, convergence/stability checks, predicted
    shares, elasticities, and robust uncertainty.
  - Move the legacy `quaids()` wrapper and manual mean-point assembly to an
    advanced or compatibility section.
- **Acceptance evidence:**
  - A new customer can obtain and validate a first result by copying one compact,
    tested block.

### PR-402 — Add a real-data, file-to-results example

- **Priority / effort:** P1 / M
- **Work:**
  - Promote the published Blanciforti dataset, or another redistributable dataset,
    into a customer-facing example.
  - Load an actual CSV with `loadd()`, map share and price columns, run preflight,
    fit the supported baseline, interpret diagnostics, compute elasticities, and
    export a table.
  - Include expected headline results and tolerances so customers can verify their
    environment.
  - Keep synthetic examples for feature isolation, but explain their purpose and
    unrealistic-share limitation once rather than throughout the suite.
- **Acceptance evidence:**
  - The example runs from the distributed artifact in a clean installation.
  - Its expected results are asserted in automated tests.

### PR-403 — Add a data-preparation guide

- **Priority / effort:** P1 / M
- **Required topics:**
  - Constructing expenditure shares and total expenditure.
  - Required price and expenditure transformations and unit consistency.
  - Good/category ordering across shares and prices.
  - Missing values, invalid observations, zeros, and corner solutions.
  - Instrument selection requirements and weak-instrument diagnostics.
  - Demographic intercept shifters.
  - Sampling weights, clusters, replicate weights, strata, and unsupported survey
    features.
  - Minimum sample/design-size considerations and recommended preflight checks.
- **Acceptance evidence:**
  - The real-data example links to each relevant preparation decision.
  - The guide contains a final input-contract checklist.

### PR-404 — Add a troubleshooting and interpretation guide

- **Priority / effort:** P1 / M
- **Work:**
  - Create a symptom-to-action table for installation failures, undefined structs,
    missing includes, shape mismatches, weak instruments, invalid shares,
    non-convergence, multiple solutions, failed bootstrap replicates, and optional
    reporting dependencies.
  - Explain which result fields establish validity before coefficients,
    elasticities, welfare measures, or curvature results are reported.
  - Explain the robust-sandwich versus bootstrap tradeoff without requiring the
    reader to reconstruct it from methodology notes.
- **Acceptance evidence:**
  - Each known customer-visible failure mode has one recommended next action and a
    link to the relevant command page.

### PR-405 — Make every example location-independent and smoke-tested

- **Priority / effort:** P1 / M
- **Work:**
  - Make includes resolve relative to the program/artifact rather than the caller's
    current working directory.
  - Ensure the `pubtable` adapter example works when launched from the GAUSS GUI,
    from the examples directory, and by absolute script path.
  - Add a runner for all examples, with explicit optional-dependency skips and
    checks for GAUSS compile/execute errors.
  - Clean generated exports after tests and assert that they are absent from the
    release archive.
- **Acceptance evidence:**
  - All examples pass from an arbitrary working directory.
  - Example smoke tests run in CI and in the release gate.

## Phase 5: Add Support and Project Operations

### PR-501 — Publish customer support expectations

- **Priority / effort:** P1 / S
- **Work:**
  - Add `SUPPORT.md` with the supported GAUSS versions/platforms, where to report
    bugs, what diagnostic information to include, and expected support boundaries.
  - Add a structured bug-report issue template capturing package version, GAUSS
    version, OS, model controls, convergence fields, and a minimal reproduction.
- **Acceptance evidence:**
  - README links to support instructions.
  - A submitted issue template contains enough information to reproduce a typical
    install or convergence problem.

### PR-502 — Add contributor and security guidance

- **Priority / effort:** P2 / S
- **Work:**
  - Add `CONTRIBUTING.md` with test commands, style expectations, and the release
    verification workflow.
  - Add a concise security/private-reporting policy appropriate for a numerical
    research library.
  - Document how release artifacts and tags are produced.
- **Acceptance evidence:**
  - A new contributor can run the fast suite and understand the full release gate
    without consulting internal AI context files.

## Phase 6: Release Candidate Gate

### PR-601 — Automate a single go/no-go command

- **Priority / effort:** P0 / M
- **The gate must verify:**
  1. Clean worktree and synchronized release metadata.
  2. Manifest and public API inventory consistency.
  3. Full source suite, including bootstrap tests.
  4. Documentation and archive-link checks.
  5. All example smoke tests.
  6. Artifact build and content verification.
  7. Clean official installation.
  8. Installed-package public API test against that exact artifact.
  9. Numerical benchmark/sweep results for supported estimator tiers.
  10. Artifact checksum and final release notes.
- **Acceptance evidence:**
  - The command exits nonzero for any failed gate.
  - The release record captures tool versions, test summaries, benchmark results,
    artifact checksum, and the source commit.

## First Public Alpha Exit Criteria

All of the following must be true:

- PR-001 through PR-003 are complete.
- PR-101 through PR-103 are complete.
- PR-201 and PR-301 are complete.
- PR-401 through PR-405 are complete.
- PR-501 is complete.
- PR-601 passes against a clean installation of the exact release artifact.
- If PR-202 is incomplete, iterated AIDS and QUAIDS are explicitly experimental,
  the quick start uses the stable baseline, and the limitations explain why a
  converged flag alone is insufficient.

## Production-Readiness Exit Criteria

These are intentionally stricter than the public-alpha gate:

- PR-202 is complete and detects the known wrong-fixed-point cases.
- The convergence sweep is repeated on multiple representative DGP families and
  applied datasets, not a single synthetic family.
- Supported iterated estimators meet a predeclared success threshold with no known
  undetected wrong-solution cases in the release benchmark.
- Curvature and zero-share inference limitations are either resolved or retained
  under an explicitly narrower support tier.
- At least one independent QUAIDS cross-implementation validation is available, or
  the absence of one remains a prominent evidence limitation.

## Suggested Execution Order

1. **Release contract:** PR-001, PR-002, PR-003.
2. **Installation spike:** PR-101 and PR-102 before further packaging work.
3. **Truthful numerical posture:** PR-201, then begin PR-202 in parallel with
   customer documentation.
4. **Documentation correction:** PR-103, PR-301, PR-302, PR-303.
5. **Implementation enablement:** PR-401 through PR-405.
6. **Support operations:** PR-501 and PR-502.
7. **Release gate and candidate:** PR-601, followed by a clean-machine release
   rehearsal and final go/no-go review.

## Deferred Follow-Up Backlog

The following work is valuable but should not displace the first-release blockers:

- Full nonlinear-feedback-corrected robust sandwich covariance.
- Formal survey strata and finite-population correction.
- Built-in JK/BRR/Fay design constructors rather than caller-supplied replicate
  columns only.
- Independent published-data QUAIDS and curvature-imposition validation.
- Additional real datasets and domain-specific recipes.
- Cross-platform automation beyond the declared first-release support matrix.
