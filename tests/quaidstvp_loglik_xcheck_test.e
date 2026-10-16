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
** Checks per DIAGONAL-H case (1, 2): (a) every period's log-likelihood
** contribution, and (b) the total, match statsmodels within TOL; (c) the
** objective ssFitTVP()'s CMLMT pass actually maximizes -- sslib's
** _tvpLFNImpl() driven through this package's _quaidsTVPQUpdate() exactly
** as _quaidsTVPMLEFit() configures it -- returns the same per-period
** vector as (a), so the cross-checked number is the MLE's own objective.
**
** Case 3 (NON-diagonal H) is a KNOWN UPSTREAM DEFECT characterization, NOT
** a validation. sslib v2.1.0's kalmanFilterDiffuseTVP() decorrelates a
** non-diagonal H by inv(chol(H)') row scaling and then applies an absolute
** per-element diffuse threshold (F_inf > 1e-12). On this fixture the last
** diffuse direction is removed by element 1 of period 4, and element 2's
** rounding residue (inflated by the 1/sqrt(H) scaling) passes 1e-12, so
** it takes a spurious diffuse update: about +14 at period 4, then a
** diverging state path. Hand-verified in
** fixtures/tvp_loglik/sslib_nondiag_h_repro.py. Only the pre-divergence
** periods 1-3 are asserted, plus check (c). The full-sample discrepancy
** is printed as KNOWN-UPSTREAM-DEFECT and is not counted. Non-diagonal H is
** therefore NOT validated. It stays a promotion blocker until
** gauss-state-space fixes it and a regression is added (see
** dev/TVP_PROMOTION_ROADMAP.md). If an sslib fix lands, statsmodels
** already holds the correct full-sample answer for case 3 in the
** committed reference.
**
** Tolerance. TOL = 1e-9 absolute, per period and on the total. Totals are
** 281-691 in magnitude, so this is about 1.5e-12 relative. Measured
** agreement (sslib v2.1.0 / GAUSS 26.1.4 vs statsmodels 0.15.0, see
** fixtures/tvp_loglik/README.md) is at most 4.0e-13 per period and
** 6.3e-13 on a total. That is double-precision rounding in a 120-period,
** 7-state recursion whose v^2/F terms divide by F ~ 1e-5, summed in a
** different order. TOL leaves more than three orders of magnitude of
** headroom for other BLAS builds and platforms. It is still about nine
** orders below any real convention or logic difference: dropping one
** ln(2*pi) term moves a period by 0.92, and a spurious diffuse update
** moves one by about 14.
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
    isDiag = H[1, 2] == 0;

    print "case " $+ ftocv(c, 1, 0) $+ ": GAUSS llf = " $+ sprintf("%.13f", sumc(rslt.loglikelihood))
        $+ ", statsmodels llf = " $+ sprintf("%.13f", sumc(llRef));
    print "        max |per-period diff| = " $+ sprintf("%.3e", dPer)
        $+ ", |total diff| = " $+ sprintf("%.3e", dTot);

    if isDiag;
        call check(dPer < TOL,
            "case " $+ ftocv(c, 1, 0) $+ " (diagonal H): every per-period log-likelihood matches statsmodels within 1e-9");
        call check(dTot < TOL,
            "case " $+ ftocv(c, 1, 0) $+ " (diagonal H): total log-likelihood matches statsmodels within 1e-9");
    else;
        nPre = 3;
        call check(maxc(abs(rslt.loglikelihood[1:nPre] - llRef[1:nPre])) < TOL,
            "case " $+ ftocv(c, 1, 0) $+ " (non-diagonal H): pre-divergence periods 1-3 match statsmodels within 1e-9");
        dFlag = abs(rslt.loglikelihood - llRef) .> TOL;
        if sumc(dFlag) > 0;
            firstBad = minc(selif(seqa(1, 1, tobs), dFlag));
            print "KNOWN-UPSTREAM-DEFECT  case " $+ ftocv(c, 1, 0) $+ " (non-diagonal H, sslib kalmanFilterDiffuseTVP): first period off = "
                $+ ftocv(firstBad, 1, 0) $+ ", GAUSS - statsmodels total = " $+ sprintf("%.10f", sumc(rslt.loglikelihood) - sumc(llRef))
                $+ ", max |per-period diff| = " $+ sprintf("%.6e", dPer) $+ " (not counted; NOT validated)";
        else;
            print "NOTE  case " $+ ftocv(c, 1, 0) $+ " (non-diagonal H) now matches statsmodels over the full sample -- the upstream defect appears fixed; revalidate and promote this case to a counted check.";
        endif;
    endif;

    /* (c): the objective ssFitTVP()'s untransformed CMLMT pass evaluates,
       configured exactly as _quaidsTVPMLEFit() does (diffuse init,
       _quaidsTVPQUpdate writing diag(q) into Q). transformed = 0 so p is
       q itself; no analytic gradient (ssControlCreate's gradproc = 0). */
    ssCtl = ssControlCreate(k_states, n1, k_states, tobs);
    ssCtl.init = "diffuse";
    ssCtl.transformed = 0;
    ssCtl.tvp = _quaidsTVPBuildModel(Zarr, eye(k_states), H, n1);
    mm = _tvpLFNImpl(q, y, &_quaidsTVPQUpdate, &init_diffTVP, 0, 0, ssCtl);
    call check(maxc(abs(mm.function - rslt.loglikelihood)) == 0,
        "case " $+ ftocv(c, 1, 0) $+ ": the MLE objective (sslib _tvpLFNImpl via _quaidsTVPQUpdate) equals the filter path exactly");

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
