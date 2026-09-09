# Support

## Supported environment

| Component | Release status |
| --- | --- |
| GAUSS 26 or later | Supported; release tests currently use GAUSS 26.1.4 |
| Windows | Supported and tested release platform |
| macOS and Linux | Expected to work for installed-library use, but not release-tested or supported for `0.1.0` |

See [README's Requirements section](README.md#requirements) for the full
table, including which optional packages (`optmt`, `pubtable`) are needed
for which features.

## Where to report a problem

Open an issue on GitHub:
[github.com/ec78/gauss-aids-model/issues](https://github.com/ec78/gauss-aids-model/issues).
The bug-report template is selected automatically when you open a new
issue and asks for the information below directly.

**Before filing**, check whether your question is already answered:

- [Troubleshooting and Interpretation Guide](docs/TROUBLESHOOTING_GUIDE.md) --
  a symptom-to-action table for installation, include, shape, and
  convergence problems, plus which result fields establish validity
  before trusting output.
- [Data Preparation Guide](docs/DATA_PREPARATION_GUIDE.md) -- getting raw
  data into this library's expected shape.
- [Model & Feature Support Tiers](README.md#model--feature-support-tiers) /
  [Feature Support Matrix](docs/FEATURE_SUPPORT_MATRIX.md) -- some
  limitations (non-convergence rates for the iterated estimator, curvature
  standard-error boundary issues, no independent QUAIDS reference
  implementation, and others) are already known, measured, and documented
  -- encountering one of these is not itself a new bug to report, though a
  case that contradicts the documented behavior is.

## What to include in a bug report

A useful report includes:

- **Package version** -- `package.json`'s `"version"` field (source repo),
  or the installed package's own `package.json` (`<GaussHome>\pkgs\quaids\package.json`).
- **GAUSS version** -- printed at the top of every batch run (e.g.
  `GAUSS 26.1.4`).
- **Operating system.**
- **The model controls you set** -- every non-default `quaidsControl`
  (`aCtl`) field: `linear`, `maxiter`, `homogenous`, `relax`, `alpha0`,
  `err`, and any `weight`/`clusterId`/`replicateWeights`/`scaleFactor`
  arguments passed to the specific procedure you called.
- **For a convergence or unexpected-result report**: the relevant
  diagnostic fields -- `qOut.converged`/`iterations`/`finalErr` (or the
  equivalent fields on whichever struct/proc you used -- `zOut`/`cOut`/
  `wfOut` all carry their own). See the [Troubleshooting
  Guide](docs/TROUBLESHOOTING_GUIDE.md#what-establishes-a-results-validity)
  for the full list per proc.
- **A minimal, runnable reproduction.** Prefer synthetic data (this
  library's own `examples/example_data.src` generator, or a small
  hand-built fixture) over your real data; if real data is required to
  reproduce the problem, describe its shape (rows, goods, whether zero
  shares/missing values are present) rather than attaching it.

The issue template captures all of this directly -- filling it in
completely is usually enough for a maintainer to reproduce a typical
install or convergence problem without back-and-forth.

## Support boundaries

This is a public alpha (`0.1.0`) maintained by a single author, not a
commercially supported product -- see [Compatibility
Policy](README.md#compatibility-policy) for what "alpha" means for this
library's own API stability. Issues are addressed on a best-effort basis;
there is no service-level agreement or guaranteed response time.

For security-sensitive reports, see
[CONTRIBUTING.md](CONTRIBUTING.md#security)'s private-reporting
instructions instead of a public issue.
