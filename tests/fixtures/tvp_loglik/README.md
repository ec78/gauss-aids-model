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
| `sslib_nondiag_h_repro.py` | Hand verification of the case 3 disagreement (see below). |

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

## Result (2026-10-16)

Environment:
- GAUSS 26.1.4.5565 with `GAUSS26_CFG=tests\gauss26_cfg_override`;
- sslib v2.1.0 (`7d0e042`), `verify_sslib_pin.ps1 -Strict` passing;
- Python 3.12.0, numpy 2.5.1, statsmodels 0.15.0 (Windows 11).

| Case | GAUSS llf | statsmodels llf | max per-period diff | total diff | Status |
|---|---|---|---|---|---|
| 1 | 690.6901320251270 | 690.6901320251268 | 4.97e-14 | 2.27e-13 | PASS |
| 2 | 281.2743884578091 | 281.2743884578097 | 3.97e-13 | 6.25e-13 | PASS |
| 3 | 688.0306020430203 | 674.2112678504466 | 14.19 (period 4) | 13.82 | **known upstream defect, not validated** |

For every case, check (c) holds exactly (difference 0). It is a
**wiring check, not a numerical oracle.** It runs sslib's `_tvpLFNImpl()`
with `transformed = 0`, through `_quaidsTVPQUpdate()` and the diffuse
init `_quaidsTVPMLEFit()` uses, and gets the same per-period vector as
the filter path.
- That is the untransformed path `ssFitTVP()` takes on its second
  (covariance) CMLMT pass.
- It does not exercise the transformed first pass, which does the
  maximizing (`positive_vars` is set, so `sstransformParms()` runs
  first).
- Both sides call the same `kalmanFilterDiffuseTVP()`.

The statsmodels comparisons in cases 1 and 2 are the independent
numerical check.

**Tolerance: 1e-9 absolute**, per period and on the total. That is about
1.5e-12 relative to totals of 281-691.
- The worst measured difference (6.3e-13) is double-precision rounding in
  a 7-state, 120-period recursion with `F ~ 1e-5`, summed in a different
  order.
- 1e-9 allows more than three orders of magnitude of headroom for other
  platforms and BLAS builds.
- It is still about nine orders below a real convention or logic error:
  a missing `ln 2pi` term moves a period by 0.92, and a spurious diffuse
  update moves one by about 14.
- No existing tolerance or reference value was changed to get this result.

### Case 3: non-diagonal H, a known sslib v2.1.0 defect

`kalmanFilterDiffuseTVP()` handles a non-diagonal H as follows:
1. It decorrelates by `Linv = inv(chol(H)')` row scaling.
2. It then applies an **absolute** per-element diffuse test,
   `F_inf > 1e-12`.
3. It ends the diffuse phase only at the start of a period
   (`_kfDiffusePhaseOver`).

On this fixture, element 1 of period 4 removes the last diffuse direction.
Element 2's rounding residue is inflated by the `1/sqrt(H)` scaling to
about 2e-11, which passes the 1e-12 test. Element 2 therefore takes a
spurious diffuse update (about `-0.5 ln F_inf`), and the state path
diverges from that point.

`sslib_nondiag_h_repro.py` reproduces this in plain numpy, using the
same algorithm and threshold:
- the period 4 element 2 residue is `F_inf = 2.334e-11`;
- the total is 688.26. sslib gives 688.03; the exact jump depends on how
  the residue rounds.

The same recursion matches statsmodels to about 1.7e-13 when either:
- it uses LDL decorrelation (residue -7e-16); or
- it keeps Cholesky but uses a threshold of 1e-10 or 1e-8.

So this is a defect in the upstream implementation, not a convention
ambiguity.

Diagonal H is never row-scaled and is unaffected. That covers cases 1 and
2, every other committed TVP test, example 14 and the documented
`diagrv(...)` H recipe. Case 3 of this test is the only non-diagonal-H
case, and it stays excluded and unvalidated.

**When it can trigger.** The mechanism needs the diffuse phase to end
**mid-period**. The common way this happens is `k_states mod n1 != 0`,
i.e. an even `n1`, since each observation removes one diffuse direction.
That is a sufficient, common condition, not a necessary one: a period
whose `Z_t` is rank-deficient can also leave the phase ending partway
through. It is also not guaranteed to fire, because whether the residue
passes `1e-12` depends on rounding. Treat any non-diagonal H as
unvalidated, whatever `n1` is.

**What the test asserts for case 3.**
- Periods 1-3 match statsmodels within 1e-9 (before the divergence).
- Two counted checks **pin** the defect so it cannot change silently:
  - the first period off by more than 1e-9 is 4;
  - the GAUSS-minus-statsmodels total is `13.8193341925737059` within
    1e-9. That is the gap GAUSS 26.1.4 with sslib v2.1.0 gives on this
    machine, printed with `%.16f`.
- These two checks characterize the defect. They do **not** validate
  any numbers.
- The gap's size comes from a rounding residue, so it reproduces only on
  a fixed build. The tight pin is on purpose: any change to the upstream
  baseline should force a re-review rather than pass.
- If the gap disappears, the test prints an explicit `FAIL` saying the
  upstream baseline changed. The first-period pin then fails, and this
  characterization and the promotion gate must be updated.
- The full-sample gap is still printed as a `KNOWN-UPSTREAM-DEFECT` line.

**Non-diagonal H is not validated.** It remains a TVP promotion blocker
until gauss-state-space fixes the threshold and a regression is added.
The committed statsmodels column already holds the correct answer for
that regression.

## Commands

```powershell
# GAUSS check (from tests/)
$env:GAUSS26_CFG = (Resolve-Path tests\gauss26_cfg_override).Path
cd tests; C:\gauss26\tgauss.exe -b -x quaidstvp_loglik_xcheck_test.e   # ALL 13 CHECKS PASSED

# Reference (from tests/fixtures/tvp_loglik/; needs numpy + statsmodels)
python statsmodels_reference.py --check   # 0.000e+00 vs committed
python statsmodels_reference.py           # regenerate (rewrites .txt/.json)
python sslib_nondiag_h_repro.py           # case 3 hand verification

# Fixture provenance (no sslib needed)
C:\gauss26\tgauss.exe -b -x export_tvp_loglik_fixture.e
```
