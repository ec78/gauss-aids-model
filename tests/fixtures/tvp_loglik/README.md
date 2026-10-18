# TVP-AIDS fixed-Q log-likelihood cross-check (promotion P5)

This is the independent second check that CLAUDE.md requires for the
TVP-AIDS likelihood. Q and H are held fixed (no optimization). The exact
diffuse log-likelihood from the GAUSS production path is compared with
**statsmodels**' exact diffuse Kalman filter, a separately authored,
widely used implementation, on a committed deterministic fixture.

Python is **validation only**. Nothing in `src/` or the installed `quaids`
package uses it, and the GAUSS test reads the committed reference without
running Python.

## Files

| File | What it is |
|---|---|
| `tvp_loglik_fixture.txt` | Data: 120 periods of `t w1 w2 w3 prel1 prel2 lx` (`%.17g`). This is the P4 Stage 6 DGP, `_quaidsTVPRawSyntheticDGP()` with n = 3, n1 = 2, k_states = 7, tobs = 120, seed 7, and the `state0`/`trueQ`/`H` of `tests/quaidstvpfit_test.e`. `prel_j` = log price j minus log price 3; `lx` = Stone-deflated log expenditure. |
| `export_tvp_loglik_fixture.e` | GAUSS script that wrote the fixture. It documents provenance; the committed `.txt` is authoritative. |
| `tvp_loglik_cases.txt` | The fixed `(diag Q, H)` cases, hand-written and read by both sides. |
| `statsmodels_reference.py` | The independent reference generator. `--check` recomputes and compares to the committed output. |
| `tvp_loglik_reference.txt` | statsmodels' per-period log-likelihood, one column per case (`repr` floats). |
| `tvp_loglik_reference.json` | Provenance: versions, platform, totals, number of diffuse periods. |
| `sslib_nondiag_h_repro.py` | Historical: hand verification of the case 3 disagreement under sslib v2.1.0, since fixed upstream (see below). |

The GAUSS check is `tests/quaidstvp_loglik_xcheck_test.e`, registered in
`tests/run_source_tests.ps1`'s sslib group (skipped by `-SkipTVPKalman`).

## Cases

| Case | Q diagonal (alpha1, alpha2, g11, g12, g22, beta1, beta2) | H | Role |
|---|---|---|---|
| 1 | 1e-5, 1e-5, 1e-4 x 5 | 2.5e-5 I | DGP truth |
| 2 | 2.5e-6 x 7 | 2.5e-5 I | the q0 start `quaidsTVPFit` test uses (0.1 mean diag H) |
| 3 | 3e-5, 7e-6, 2e-4, 5e-5, 1.5e-4, 8e-5, 4e-5 | [3e-5 1e-5; 1e-5 2e-5] | heterogeneous Q, **non-diagonal** H |

## Why the reference is independent

- **Filter:** statsmodels'
  `statsmodels.tsa.statespace.kalman_filter.KalmanFilter` with
  `initialize_diffuse()`, which is compiled Cython code. It shares no code
  with sslib's `kalmanFilterDiffuseTVP()` (gauss-state-space).
- **Design matrices:** the Python side rebuilds every Z_t from `prel`/`lx`
  and the documented state layout `[alpha | gamma (unique, row-major upper
  triangle) | beta]`. It does not read `_quaidsTVPBuildZ()` output, so Z
  construction is checked as well.
- **Inputs:** data, Q and H come from the committed text files. Neither
  side computes the other's numbers.

**Convention.** Both libraries use the univariate exact diffuse
likelihood (Durbin & Koopman 2012, sections 5.2 and 6.4):
- a diffuse element contributes `-0.5 (ln 2pi + ln F_inf)`;
- an ordinary element contributes `-0.5 (ln 2pi + ln F + v^2/F)`;
- every period is summed.

The script sets `loglikelihood_burn = 0`, because statsmodels' `MLEModel`
would otherwise drop the diffuse periods. The convention was checked on a
local-level toy model, where period 1 gives `-0.5 ln 2pi` and period 2
matches the hand formula. statsmodels reports `nobs_diffuse = 4` for every
case: 7 diffuse states, 2 observations per period.

## Result (2026-10-17, sslib v2.1.1)

Environment:
- GAUSS 26.1.4.5565 with `GAUSS26_CFG=tests\gauss26_cfg_override`;
- sslib v2.1.1 (`bf827ae`), `verify_sslib_pin.ps1 -Strict` passing;
- Python 3.12.0, numpy 2.5.1, statsmodels 0.15.0 (Windows 11).

| Case | GAUSS llf | statsmodels llf | max per-period diff | total diff | Status |
|---|---|---|---|---|---|
| 1 | 690.6901320251270 | 690.6901320251268 | 4.97e-14 | 2.27e-13 | PASS |
| 2 | 281.2743884578091 | 281.2743884578097 | 3.97e-13 | 6.25e-13 | PASS |
| 3 | 674.2112678504466 | 674.2112678504466 | 1.71e-13 | 0.00e+00 | PASS |

All three cases are counted full-sample checks: every period and the total
within the tolerance below. Cases 1 and 2 are bit-identical to the v2.1.0
result. Case 3 is the non-diagonal-H case; it failed on v2.1.0 (see
"Case 3 history" below).

For every case, check (c) holds exactly (difference 0). It is a
**wiring check, not a numerical oracle.** It runs sslib's `_tvpLFNImpl()`
with `transformed = 0`, through `_quaidsTVPQUpdate()` and the diffuse
init `_quaidsTVPMLEFit()` uses, and gets the same per-period vector as
the filter path.
- That is the objective CMLMT would maximize for an untransformed fit.
- `quaidsTVPFit()` sets `positive_vars`, so its one CMLMT pass is
  transformed (`sstransformParms()` runs first). That path is not
  exercised here.
- Up to v2.1.0, `ssFitTVP()` also re-ran CMLMT untransformed for the
  covariance. Since v2.1.1 it converts the covariance by the delta method
  instead, so `quaidsTVPFit()` runs no untransformed CMLMT pass.
- Both sides call the same `kalmanFilterDiffuseTVP()`.

The statsmodels comparisons in cases 1-3 are the independent numerical
check.

**Tolerance: 1e-9 absolute**, per period and on the total. That is about
1.5e-12 relative to totals of 281-691.
- The worst measured difference (6.3e-13) is double-precision rounding in
  a 7-state, 120-period recursion with `F ~ 1e-5`, summed in a different
  order.
- 1e-9 allows more than three orders of magnitude of headroom for other
  platforms and BLAS builds.
- It is still about nine orders below a real convention or logic error:
  a missing `ln 2pi` term moves a period by 0.92, and the v2.1.0
  spurious diffuse update moved one by about 14.
- No tolerance or reference value was changed to get this result. The
  statsmodels reference (`--check`: 0.000e+00 vs committed) is the same
  file that previously showed the case 3 disagreement.

### Case 3 history: a sslib v2.1.0 defect, fixed in v2.1.1

In sslib v2.1.0, `kalmanFilterDiffuseTVP()` handled a non-diagonal H as
follows:
1. It decorrelated by `Linv = inv(chol(H)')` row scaling.
2. It then applied an **absolute** per-element diffuse test,
   `F_inf > 1e-12`, in those scaled units.
3. It ended the diffuse phase only at the start of a period
   (`_kfDiffusePhaseOver`).

On this fixture, element 1 of period 4 removes the last diffuse direction.
Element 2's rounding residue was inflated by the `1/sqrt(H)` scaling to
about 2e-11, which passed the 1e-12 test. Element 2 therefore took a
spurious diffuse update (about `-0.5 ln F_inf`), and the state path
diverged from that point. GAUSS gave 688.0306020430203, a total gap of
13.8193341925737059 against statsmodels, first off at period 4.

`sslib_nondiag_h_repro.py` reproduces the v2.1.0 behaviour in plain numpy,
using the same algorithm and threshold. It does not call sslib, so it
still shows the old behaviour:
- the period 4 element 2 residue is `F_inf = 2.334e-11`;
- the total is 688.26. sslib v2.1.0 gave 688.03; the exact jump depends
  on how the residue rounds.

The same recursion matches statsmodels to about 1.7e-13 when either:
- it uses LDL decorrelation (residue -7e-16); or
- it keeps Cholesky but uses a threshold of 1e-10 or 1e-8.

The defect was reported upstream and fixed in gauss-state-space `004c1be`,
released as v2.1.1 (`bf827ae`). All four exact diffuse filters now test
`F_inf * L[j,j]^2 > 1e-12`, i.e. in LDL-decorrelation units (statsmodels'
convention, and the original units for a diagonal H). They also count
diffuse directions from `rank(P_inf_0)` and set `P_inf` to exactly zero
once the count reaches zero within a period. gauss-state-space added its
own regression for this (`test/sskalmandiffusecorrh.inc`).

Before the repin, this test pinned the defect with two counted
characterization checks (first period off = 4; total gap =
13.8193341925737059 within 1e-9), designed to fail on any upstream change.
On v2.1.1 both failed as intended, and case 3 was converted to the same
full-sample checks as cases 1 and 2 (13 checks became 12).

**Scope of what this validates.** With this non-diagonal H, the fixed-Q
exact diffuse log-likelihood matches an independent implementation over
the full sample. It does not by itself validate estimated Q, or settle the
support tier of non-diagonal H in `quaidsTVPFit()`. That judgment belongs
to the TVP promotion release gate (P8, see `dev/TVP_PROMOTION_ROADMAP.md`).

## Commands

```powershell
# GAUSS check (from tests/)
$env:GAUSS26_CFG = (Resolve-Path tests\gauss26_cfg_override).Path
cd tests; C:\gauss26\tgauss.exe -b -x quaidstvp_loglik_xcheck_test.e   # ALL 12 CHECKS PASSED

# Reference (from tests/fixtures/tvp_loglik/; needs numpy + statsmodels)
python statsmodels_reference.py --check   # 0.000e+00 vs committed
python statsmodels_reference.py           # regenerate (rewrites .txt/.json)
python sslib_nondiag_h_repro.py           # case 3 hand verification (v2.1.0 behaviour)

# Fixture provenance (no sslib needed)
C:\gauss26\tgauss.exe -b -x export_tvp_loglik_fixture.e
```
