# Contributing

Everything a new contributor needs to run the test suite and understand
the release gate, without consulting any internal AI-assistant context
file (`CLAUDE.md`, `dev/GOLD_STANDARD_TODO.md`) -- those exist for AI
coding assistants working session-to-session on this repo, not as a
contributor prerequisite.

## Prerequisites

- GAUSS 26 or later (`tgauss.exe` on your `PATH`, or note its full path).
- PowerShell (Windows; the test/build scripts are Windows-only tooling,
  not part of the library's own portable runtime API).
- Optional, only if you're touching the corresponding adapter: the
  `optmt` package (curvature imposition) and/or the `pubtable` package
  (table export), each installed separately via GAUSS's own package
  manager.

## Running the test suite

The fast, source-tree suite (`#include`-based, does not require
installing the package) -- run this before every commit:

```powershell
powershell -ExecutionPolicy Bypass -File tests\run_source_tests.ps1
```

Useful flags: `-SkipBootstrap` (skips two bootstrap-refit test files that
add ~45-50s), `-SkipCurvature`/`-SkipPubtable` (skip if you don't have
`optmt`/`pubtable` installed locally).

Smoke-test every example (needs the package installed -- see below):

```powershell
powershell -ExecutionPolicy Bypass -File tests\run_examples_smoke.ps1
```

## Style expectations

This codebase follows the original author's own established conventions
-- match them rather than introducing a different style:

- **Short, lowercase variable names** (`w`, `n`, `gg`, `b`, `u`), matching
  this codebase's terse econometrics-code style, not verbose renaming.
- **All locals declared in one `local` statement** at the top of each
  proc.
- **Struct-returning procs** are declared `proc (struct TypeName) =
  name(...);` -- no variable name in the return slot -- so a caller can
  skip pre-declaring the target variable (`qOut = quaidsFit(...);` works
  with no prior `struct quaidsOut qOut;` line). Every existing
  struct-returning proc follows this; new ones should too.
- **Loop style**: `do while ok; ... endo;`, not GAUSS's newer `for` loop.
- **Character-matrix name vectors** (built with the `0$+"X"$+ftocv(...)`
  idiom) are legacy character matrices, not the native `string array`
  type -- struct fields holding them must be declared `matrix`. Match the
  existing field's type when adding a new one to an existing struct
  pattern.
- **Comments explain WHY, not WHAT** -- a hidden constraint, a subtle
  invariant, a workaround for a specific GAUSS quirk. Avoid narrating
  implementation history (which milestone added something, dates,
  internal chronology) in customer-facing `README.md`/`docs/**/*.md`
  pages -- that belongs in `CHANGELOG.md` instead; see
  `scripts/verify_docs_quality.ps1`'s own heading/link checks and
  `scripts/verify_docs_consistency.ps1` for the automated guards on this.
- **Never trust a derived formula, fix, or script without actually
  running it.** This codebase's own history has repeatedly found real
  bugs (in already-shipped code, not just new code) purely by executing
  something and reading the output carefully, not by re-reading source.
  When you fix a real bug, add a regression test and verify it actually
  fails without your fix (temporarily revert it, confirm the test fails,
  then restore the fix) before considering the fix done.

## Documentation quality gates

Three scripts run as part of `tests/run_source_tests.ps1` (and therefore
CI) and will fail your PR if violated:

- `scripts/verify_public_api.ps1` -- release-metadata version consistency
  across `package.json`/`CITATION.cff`/`docs/public-api.json`/
  `CHANGELOG.md`, and that every documented procedure/struct exists in
  `src/`.
- `scripts/verify_docs_consistency.ps1` -- documented defaults
  (`docs/command-reference/quaidsControlCreate.md`'s table) match
  `quaidsControlCreate()`'s actual coded defaults.
- `scripts/verify_docs_quality.ps1` -- every `docs/command-reference/*.md`
  page follows the required Purpose/Format/Parameters/Returns/Remarks/
  Examples/Source/See Also heading structure, every internal doc link and
  `#anchor` resolves, and every keyword argument in a ` ```gauss ` code
  snippet matches the real procedure signature.

If you add a new public procedure or struct field, add or update its
`docs/command-reference/*.md` page and its entry in
`docs/public-api.json` in the same change.

## Release verification workflow

The full release gate, run before tagging a release (and available to
run at any time to sanity-check the pipeline):

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_release_verification.ps1 -BuildArtifact -ForceArtifact -InstallArtifact
```

This runs, in order: the full source-tree suite (no skips), builds a
release `.zip` (`scripts\build_package.ps1`, which stages `package.json`
plus the root files (`README.md`/`CHANGELOG.md`/`CITATION.cff`/`LICENSE`/
`llms.txt`) and the `src`/`docs`/`examples`/`scripts`/`tests` directories,
strips known generated test/example artifacts, and self-verifies the
result via `scripts\verify_release_artifact.ps1` -- including an
archive-level check that every relative doc link and `#anchor` resolves
*inside the built archive itself*, not just the git working tree, since
a few maintainer-only files (`CLAUDE.md`, `dev/GOLD_STANDARD_TODO.md`,
`dev/PUBLIC_RELEASE_ROADMAP.md`) are never shipped), installs the artifact
into a real GAUSS package directory (`<GaussHome>\pkgs\quaids` by
default), then runs the installed-package public API test and the
example smoke tests against that exact installed copy.

Each individual step is also runnable on its own -- see that script's own
parameters (`-BuildArtifact`, `-InstallArtifact`, `-SkipInstalledPackageTest`).

## The release gate (single go/no-go command)

`scripts\run_release_gate.ps1` wraps the workflow above into a single
command with a clean exit-code contract -- run this before tagging a
release:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_release_gate.ps1
```

In addition to everything `run_release_verification.ps1` already does
(including bootstrap tests -- this always runs the full suite, no
skips), it checks the git worktree is clean, runs the 200-seed
convergence sweep and captures its per-model summary, computes the
built artifact's SHA256 checksum, and extracts the current release notes
from `CHANGELOG.md`. It exits nonzero on the first failed check and
writes a structured JSON release record to `release_records/` (gitignored
-- regenerate per run) either way, capturing tool versions, the source
commit, the convergence-sweep results, and the artifact checksum. Pass
`-AllowDirtyWorktree` to run it against uncommitted changes (for
sanity-checking only, not before an actual release) or
`-SkipConvergenceSweep` to skip the slowest step during iteration.

## How release artifacts and tags are produced

1. Bump `package.json`'s `"version"`, `CITATION.cff`'s `version:`
   (plus `date-released:`), and `docs/public-api.json`'s `"version"`
   together, and give `CHANGELOG.md`'s top entry a matching
   `## <version> - <date>` heading -- `scripts/verify_public_api.ps1`
   enforces all four stay in sync.
2. Commit the version bump. `scripts\run_release_gate.ps1` (see above)
   requires a clean worktree as its own first check, so this has to
   happen *before* the gate runs, not after.
3. Run the release gate (`powershell -ExecutionPolicy Bypass -File
   scripts\run_release_gate.ps1`) against that clean, committed state and
   confirm it reports `GO`.
4. Tag the release commit (`git tag v<version>`, e.g. `git tag v0.1.0`)
   and push the tag (`git push origin v<version>`).
5. Create a GitHub release from that tag (`gh release create v<version>
   "quaids <version>.zip" --notes-file <changelog excerpt>`), attaching
   the built `.zip` as a release asset. The release record the gate wrote
   to `release_records/` has the artifact's SHA256 checksum, tool
   versions, and convergence-sweep results worth carrying into the
   release notes or keeping alongside the tag.

## Security

This is a numerical research library that operates on user-supplied
local matrices/dataframes -- it has no network surface, does not handle
authentication or credentials, and does not parse untrusted remote input.
The realistic security-relevant surface is narrow: malformed/adversarial
local input data crashing a GAUSS session (a robustness bug, reportable
as a normal issue -- see [SUPPORT.md](SUPPORT.md)), and the integrity of
the distributed release artifact itself.

If you believe you've found a genuine security issue (as opposed to a
robustness bug), please report it privately rather than opening a public
GitHub issue: email [eric.clower78@gmail.com](mailto:eric.clower78@gmail.com)
with a description and, if possible, a reproduction. Expect an
acknowledgment on a best-effort basis (see [SUPPORT.md](SUPPORT.md) for
this project's support boundaries) -- there is no dedicated security team
or formal disclosure SLA for a project at this stage.
