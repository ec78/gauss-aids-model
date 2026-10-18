/*
** quaidstvp_loglik_xcheck_test.e
**
** TVP-AIDS promotion P5: the second, independent check CLAUDE.md requires
** for the TVP-AIDS likelihood. At FIXED Q/H (no optimization), the exact
** diffuse log-likelihood GAUSS computes must match an independent
** implementation -- statsmodels' exact diffuse Kalman filter -- on a
** committed deterministic fixture.
**
** Inputs, all committed under tests/fixtures/tvp_loglik/ (see README.md
** there for provenance and how to regenerate):
**   tvp_loglik_fixture.txt   -- the P4 Stage 6 DGP data (n = 3 goods,
**                               n1 = 2, k_states = 7, tobs = 120, seed 7):
**                               shares, relative log prices, Stone-
**                               deflated log expenditure.
**   tvp_loglik_cases.txt     -- three fixed (Q diagonal, H) cases: the DGP
**                               truth; the flat q0 start quaidsTVPFit's
**                               test uses; heterogeneous Q with a
**                               NON-diagonal H (exercises sslib's
**                               Cholesky decorrelation branch).
**   tvp_loglik_reference.txt -- statsmodels' per-period log-likelihood
**                               for each case (statsmodels_reference.py).
**
** GAUSS side, the production path only: _quaidsTVPBuildZ() builds Z_t
** from the raw regressors, _quaidsTVPBuildModel() packages the model and
** _quaidsTVPKalmanFit(.., 1) runs sslib's kalmanFilterDiffuseTVP(). The
** statsmodels side rebuilds Z_t itself, so Z construction is checked too.
**
** Checks per case (1, 2: diagonal H; 3: NON-diagonal H): (a) every
** period's log-likelihood contribution, and (b) the total, match
** statsmodels within TOL. These are the independent numerical checks.
** Case 3 exercises sslib's decorrelation branch for a non-diagonal H.
**
** Check (c), every case, is a WIRING check, not a numerical oracle:
** sslib's _tvpLFNImpl() with transformed = 0, driven through this
** package's _quaidsTVPQUpdate() with the diffuse init _quaidsTVPMLEFit()
** uses, must return exactly the per-period vector of (a). It confirms the
** Q update and the y orientation reach the same sslib filter. It is the
** objective CMLMT would maximize for an UNtransformed fit (no
** positive_vars). quaidsTVPFit() sets positive_vars, so its one CMLMT
** pass is transformed (sstransformParms() runs first); that path is not
** exercised here. Since sslib v2.1.1, ssFitTVP() no longer re-runs CMLMT
** untransformed for the covariance (it uses the delta method instead).
** Both sides call the same kalmanFilterDiffuseTVP().
**
** Case 3 history. sslib v2.1.0 decorrelated a non-diagonal H by
** inv(chol(H)') row scaling and then applied an absolute per-element
** diffuse threshold (F_inf > 1e-12) in those units. On this fixture the
** last diffuse direction is removed by element 1 of period 4, and element
** 2's rounding residue passed 1e-12, so it took a spurious diffuse update:
** a GAUSS-minus-statsmodels total of 13.8193341925737059, first off at
** period 4 (hand-verified in fixtures/tvp_loglik/sslib_nondiag_h_repro.py).
** Up to the P5 audit follow-up this test pinned that defect. sslib v2.1.1
** (gauss-state-space 004c1be) tests F_inf in LDL-decorrelation units and
** counts diffuse directions from rank(P_inf_0); case 3 now matches
** statsmodels over the full sample and is counted like cases 1 and 2.
** The committed statsmodels reference column was not changed.
**
** Tolerance. TOL = 1e-9 absolute, per period and on the total. Totals are
** 281-691 in magnitude, so this is about 1.5e-12 relative. Measured
** agreement (sslib v2.1.1 / GAUSS 26.1.4 vs statsmodels 0.15.0, see
** fixtures/tvp_loglik/README.md) is at most 4.0e-13 per period and
** 6.3e-13 on a total, all three cases. That is double-precision rounding in a 120-period,
** 7-state recursion whose v^2/F terms divide by F ~ 1e-5, summed in a
** different order. TOL leaves more than three orders of magnitude of
** headroom for other BLAS builds and platforms. It is still about nine
** orders below any real convention or logic difference: dropping one
** ln(2*pi) term moves a period by 0.92, and the v2.1.0 spurious diffuse
** update moved one by about 14.
**
** Requires `library cmlmt, tsmt, sslib;` and the GAUSS26_CFG override
** documented in CLAUDE.md (run_source_tests.ps1 sets it for this file).
**
** Run from the tests/ directory:
**   tgauss -b -x quaidstvp_loglik_xcheck_test.e
*/

new;
library cmlmt, tsmt, sslib;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidstvp.src;
#include ../src/quaidstvpkalman.src;
#include ../src/quaidstvpmle.src;

nfail = 0;
ncheck = 0;

proc (0) = check(cond, label);
    local i;
    i = ncheck + 1;
    ncheck = i;
    if cond;
        print "PASS  " $+ label;
    else;
        print "FAIL  " $+ label;
        nfail = nfail + 1;
    endif;
endp;

TOL = 1e-9;
n1 = 2;
k_states = 2*n1 + n1*(n1+1)/2;

dat = csvReadM("fixtures/tvp_loglik/tvp_loglik_fixture.txt", 2, 1, " ");
cases = csvReadM("fixtures/tvp_loglik/tvp_loglik_cases.txt", 2, 1, " ");
ref = csvReadM("fixtures/tvp_loglik/tvp_loglik_reference.txt", 2, 1, " ");

tobs = rows(dat);
ncase = rows(cases);

call check(tobs == 120 and cols(dat) == 7 and dat[tobs, 1] == tobs,
    "fixture: 120 periods x (t, w1..w3, prel1, prel2, lx) read");
call check(ncase == 3 and cols(cases) == 11,
    "cases: 3 fixed (Q, H) cases x (case, q1..q7, h11, h12, h22) read");
call check(rows(ref) == tobs and cols(ref) == ncase + 1 and ref[., 1] == seqa(1, 1, tobs),
    "reference: 120 periods x (t, one column per case) read");

y = dat[., 2:1+n1];
prel = dat[., 5:4+n1];
lx = dat[., 7];
Zarr = _quaidsTVPBuildZ(prel, lx, n1);

struct tvpModel tvpm;
struct kalmanResult rslt;
struct ssControl ssCtl;
struct modelResults mm;

c = 1;
do while c <= ncase;
    q = cases[c, 2:1+k_states]';
    H = (cases[c, 9]~cases[c, 10])|(cases[c, 10]~cases[c, 11]);
    llRef = ref[., 1+c];

    /* (a), (b): production filter path. */
    tvpm = _quaidsTVPBuildModel(Zarr, diagrv(eye(k_states), q), H, n1);
    rslt = _quaidsTVPKalmanFit(tvpm, y', 1);
    dPer = maxc(abs(rslt.loglikelihood - llRef));
    dTot = abs(sumc(rslt.loglikelihood) - sumc(llRef));
    if H[1, 2] == 0;
        hLabel = " (diagonal H)";
    else;
        hLabel = " (non-diagonal H)";
    endif;

    print "case " $+ ftocv(c, 1, 0) $+ ": GAUSS llf = " $+ sprintf("%.13f", sumc(rslt.loglikelihood))
        $+ ", statsmodels llf = " $+ sprintf("%.13f", sumc(llRef));
    print "        max |per-period diff| = " $+ sprintf("%.3e", dPer)
        $+ ", |total diff| = " $+ sprintf("%.3e", dTot);

    call check(dPer < TOL,
        "case " $+ ftocv(c, 1, 0) $+ hLabel $+ ": every per-period log-likelihood matches statsmodels within 1e-9");
    call check(dTot < TOL,
        "case " $+ ftocv(c, 1, 0) $+ hLabel $+ ": total log-likelihood matches statsmodels within 1e-9");

    /* (c): wiring only. sslib's untransformed objective: diffuse init,
       _quaidsTVPQUpdate writing diag(q) into Q, transformed = 0 so p is q
       itself, no analytic gradient (ssControlCreate's gradproc = 0).
       Compared to the same sslib filter, so it is not a numerical oracle;
       quaidsTVPFit()'s transformed (positive_vars) CMLMT pass is not
       exercised. */
    ssCtl = ssControlCreate(k_states, n1, k_states, tobs);
    ssCtl.init = "diffuse";
    ssCtl.transformed = 0;
    ssCtl.tvp = _quaidsTVPBuildModel(Zarr, eye(k_states), H, n1);
    mm = _tvpLFNImpl(q, y, &_quaidsTVPQUpdate, &init_diffTVP, 0, 0, ssCtl);
    call check(maxc(abs(mm.function - rslt.loglikelihood)) == 0,
        "case " $+ ftocv(c, 1, 0) $+ ": wiring -- untransformed objective (sslib _tvpLFNImpl via _quaidsTVPQUpdate, transformed = 0) equals the filter path exactly");

    c = c + 1;
endo;


/* ==========================================================================
   Summary
   ========================================================================== */

if nfail == 0;
    print "ALL " $+ ftocv(ncheck, 1, 0) $+ " CHECKS PASSED";
else;
    print "" $+ ftocv(nfail, 1, 0) $+ " CHECKS FAILED";
endif;
