/*
** quaids_curvature_bootstrap_test.e
**
** Milestone 15: validates quaidsCurvatureBootstrapFit()/
** printQuaidsCurvatureBootstrap() (src/quaidscurvature.src) -- the i.i.d.
** row (pairs) bootstrap standard error for a quaidsCurvatureFit()
** coefficient vector, added to close a documented gap: the existing
** classical delta-method SE is known to be unreliable whenever the
** estimated Cholesky factor lands on the boundary of the negative-
** semidefinite cone (see quaids_curvature_test.e's own boundary-entry
** check, and src/quaidscurvature.src's header).
**
** Reuses this repo's existing curvature fixtures at their existing
** tobs=3000 (does NOT shrink tobs -- that risks changing convergence
** behavior already validated at that size). Only B is kept small, since a
** single AIDS curvature fit measures ~0.87s and a single QUAIDS curvature
** fit ~7.26s (Milestone 13's own timing) -- B=15 for AIDS (~13s) and B=5
** for QUAIDS (~36s) keeps this file's own runtime bounded while still
** exercising the resample/refit/aggregate loop for real. This is why the
** file is NOT in run_source_tests.ps1's default test list -- see that
** script's -SkipBootstrap switch, added alongside this file, matching the
** existing -SkipPubtable/-SkipCurvature opt-out pattern.
**
** Does not check exact algebraic identities (there are none to check for
** a bootstrap SE) -- instead checks completion-count bookkeeping,
** shape/finiteness of the reported SEs, and a same-order-of-magnitude
** plausibility comparison against the existing delta-method SE, EXCLUDING
** entries at a boundary (near-zero) Cholesky-factor position, reusing
** quaids_curvature_test.e's own `minc(abs(vech(cOut.cholA))) < 1e-6`
** boundary-detection idiom -- boundary entries are exactly where the two
** SEs are expected to diverge (that divergence is the whole reason this
** milestone exists), so they are excluded from the "same order of
** magnitude" comparison rather than silently making the comparison
** vacuous by including them.
**
** Uses a fixed seed (42) throughout for full reproducibility.
**
** Milestone 18 adds quaidsCurvatureBootstrapCI() percentile-CI checks
** (shape, ordering, containment of the point estimate, and a direct
** quantile() cross-check at a specific flattened index), plus a
** regression guard confirming seBoot's individual cells are positioned
** correctly relative to bootOut.b -- this caught a real, silent bug
** (GAUSS's reshape() fills row-major, not column-major like vec()) that
** had scrambled cOut.se/bootOut.seBoot's cell positions since Milestone
** 10/15, invisible to the shape/sign/finiteness checks already in place
** (those are permutation-invariant). See GOLD_STANDARD_TODO.md's
** Milestone 18 section.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_curvature_bootstrap_test.e
*/

new;
library optmt;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidscurvature.src;
#include quaidsfixtures.src;

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

/* ===================================================================
   AIDS block
   =================================================================== */

tobs = 3000;
{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsCurvatureSyntheticDGP(tobs);

aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

cOut = quaidsCurvatureFit(qOut, w, prices, totexp, aCtl);
call check(cOut.converged == 1, "AIDS: base curvature fit converged (precondition)");

B = 15;
seed = 42;

bootOut = quaidsCurvatureBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, B, seed);

call check(bootOut.nRequested == B, "AIDS: nRequested echoes B");
call check(bootOut.nCompleted >= 1, "AIDS: at least one bootstrap replication completed");
call check(bootOut.nCompleted == B, "AIDS: all requested replications completed (no failures on this well-behaved fixture)");
call check(bootOut.nFailed + bootOut.nCompleted <= bootOut.nAttempts, "AIDS: completed+failed <= attempts");
call check(bootOut.seed == seed, "AIDS: seed is echoed back");

call check(rows(bootOut.b) == rows(cOut.b) and cols(bootOut.b) == cols(cOut.b), "AIDS: bootOut.b matches cOut.b's shape");
call check(sumc(sumc(bootOut.b .== cOut.b)) == rows(cOut.b)*cols(cOut.b), "AIDS: bootOut.b is exactly the base (unresampled) point estimate");
call check(sumc(sumc(bootOut.seDelta .== cOut.se)) == rows(cOut.se)*cols(cOut.se), "AIDS: bootOut.seDelta is exactly the base delta-method SE");

call check(rows(bootOut.seBoot) == rows(cOut.b) and cols(bootOut.seBoot) == cols(cOut.b), "AIDS: seBoot matches cOut.b's shape");
call check(minc(minc(bootOut.seBoot)) >= 0, "AIDS: seBoot are all non-negative");
call check(sumc(sumc(bootOut.seBoot .== bootOut.seBoot)) == rows(bootOut.seBoot)*cols(bootOut.seBoot), "AIDS: seBoot contains no NaN/missing values");

call check(rows(bootOut.bBoot) == bootOut.nCompleted and cols(bootOut.bBoot) == rows(cOut.b)*cols(cOut.b),
    "AIDS: bBoot has one row per completed replication and one column per coefficient");

/* Milestone 18 regression guard: seBoot's INDIVIDUAL CELLS must match
   an independent per-column stdc(bBoot) at the correct vec(bootOut.b)-
   order position -- not just have the right shape/sign/finiteness
   (permutation-invariant checks that would still pass even if every SE
   were scrambled to the wrong coefficient, exactly what a real, silent
   reshape bug did here before it was caught). */
seBootFromStdc = reshape(stdc(bootOut.bBoot), cols(bootOut.b), rows(bootOut.b))';
call check(maxc(maxc(abs(bootOut.seBoot - seBootFromStdc))) < 1e-10,
    "AIDS: seBoot's cells are correctly positioned relative to bootOut.b (not scrambled by row/column-major reshape mismatch)");

/* Plausibility check -- NOT a per-cell "same order of magnitude"
   comparison against seDelta. Building this test found that the
   delta-method SE's boundary-inference unreliability does not stay
   confined to the specific gamma row/column tied to a boundary
   (near-zero) Cholesky entry: the classical NLS covariance is for the
   WHOLE vech(A) vector at once, and the numerical-Jacobian delta method
   mixes every vech(A) parameter into every reported coefficient's SE, so
   a single boundary entry can inflate seDelta anywhere in the vector
   (confirmed directly on this fixture: seDelta reaches into the hundreds
   on rows with no boundary-adjacent Cholesky entry of their own). A
   per-cell magnitude comparison excluding only the directly-tied rows is
   therefore not a reliable check. Instead, check the concrete thing this
   milestone is actually for: seBoot itself stays well-behaved (small,
   finite, non-negative) where seDelta does not -- i.e. the bootstrap SE
   does not inherit the delta method's boundary blow-up. */
call check(maxc(maxc(bootOut.seBoot)) < 2.0,
    "AIDS: seBoot stays well-behaved (bounded) even though seDelta reaches into the hundreds on this fixture");

/* Percentile confidence intervals (Milestone 18). */
alphaCI = 0.05;
{ ciLower, ciUpper } = quaidsCurvatureBootstrapCI(bootOut, alphaCI);
call check(rows(ciLower) == rows(bootOut.b) and cols(ciLower) == cols(bootOut.b), "AIDS: ciLower matches bootOut.b's shape");
call check(rows(ciUpper) == rows(bootOut.b) and cols(ciUpper) == cols(bootOut.b), "AIDS: ciUpper matches bootOut.b's shape");
call check(minc(minc(ciUpper - ciLower)) >= 0, "AIDS: ciUpper >= ciLower elementwise");
call check(sumc(sumc((ciLower .<= bootOut.b) .and (bootOut.b .<= ciUpper))) == rows(bootOut.b)*cols(bootOut.b),
    "AIDS: the point estimate falls within its own percentile CI at every cell (expected on this well-behaved fixture)");

/* Direct cross-check against a hand-computed quantile at one specific
   flattened (vec-order) position, confirming the reshape/index mapping
   -- not just that the proc runs. Position k=5 in vec(bootOut.b)-order
   maps to row ((k-1) mod rowsB)+1, col floor((k-1)/rowsB)+1. */
rowsB1 = rows(bootOut.b);
kCheck = 5;
iCheck = ((kCheck-1) % rowsB1) + 1;
jCheck = floor((kCheck-1)/rowsB1) + 1;
qLoDirect = quantile(bootOut.bBoot[., kCheck], alphaCI/2);
call check(abs(qLoDirect - ciLower[iCheck, jCheck]) < 1e-10,
    "AIDS: ciLower at a specific flattened index matches a direct quantile() call at the correctly-mapped (row, col)");

call printQuaidsCurvatureBootstrap(bootOut);
call check(1, "AIDS: printQuaidsCurvatureBootstrap() runs without error");


/* ===================================================================
   QUAIDS block
   =================================================================== */

tobsQ = 3000;
{ wQ, intcptQ, pricesQ, totexpQ, instrQ, trueParamsQ } = _quaidsSyntheticDGP(tobsQ, 204, 1, 1);

aCtlQ = quaidsControlCreate();
aCtlQ.linear = 0;
aCtlQ.maxiter = 100;
aCtlQ.homogenous = 1;
aCtlQ.err = .0001;
aCtlQ.relax = .25;

qOutQ = quaidsFit(wQ, intcptQ, pricesQ, totexpQ, instrQ, aCtlQ);

cOutQ = quaidsCurvatureFit(qOutQ, wQ, pricesQ, totexpQ, aCtlQ);
call check(cOutQ.converged == 1, "QUAIDS: base curvature fit converged (precondition)");

BQ = 5;

bootOutQ = quaidsCurvatureBootstrapFit(wQ, intcptQ, pricesQ, totexpQ, instrQ, aCtlQ, BQ, seed);

call check(bootOutQ.nRequested == BQ, "QUAIDS: nRequested echoes B");
call check(bootOutQ.nCompleted >= 1, "QUAIDS: at least one bootstrap replication completed");
call check(bootOutQ.nFailed + bootOutQ.nCompleted <= bootOutQ.nAttempts, "QUAIDS: completed+failed <= attempts");

call check(rows(bootOutQ.b) == rows(cOutQ.b) and cols(bootOutQ.b) == cols(cOutQ.b), "QUAIDS: bootOutQ.b matches cOutQ.b's shape (includes lambda row)");
call check(sumc(sumc(bootOutQ.b .== cOutQ.b)) == rows(cOutQ.b)*cols(cOutQ.b), "QUAIDS: bootOutQ.b is exactly the base (unresampled) point estimate");

call check(rows(bootOutQ.seBoot) == rows(cOutQ.b) and cols(bootOutQ.seBoot) == cols(cOutQ.b), "QUAIDS: seBoot matches cOutQ.b's shape");
call check(minc(minc(bootOutQ.seBoot)) >= 0, "QUAIDS: seBoot are all non-negative");
call check(sumc(sumc(bootOutQ.seBoot .== bootOutQ.seBoot)) == rows(bootOutQ.seBoot)*cols(bootOutQ.seBoot), "QUAIDS: seBoot contains no NaN/missing values");

call check(rows(bootOutQ.bBoot) == bootOutQ.nCompleted and cols(bootOutQ.bBoot) == rows(cOutQ.b)*cols(cOutQ.b),
    "QUAIDS: bBoot has one row per completed replication and one column per coefficient");

seBootFromStdcQ = reshape(stdc(bootOutQ.bBoot), cols(bootOutQ.b), rows(bootOutQ.b))';
call check(maxc(maxc(abs(bootOutQ.seBoot - seBootFromStdcQ))) < 1e-10,
    "QUAIDS: seBoot's cells are correctly positioned relative to bootOutQ.b (not scrambled by row/column-major reshape mismatch)");

{ ciLowerQ, ciUpperQ } = quaidsCurvatureBootstrapCI(bootOutQ, alphaCI);
call check(rows(ciLowerQ) == rows(bootOutQ.b) and cols(ciLowerQ) == cols(bootOutQ.b), "QUAIDS: ciLower matches bootOutQ.b's shape");
call check(rows(ciUpperQ) == rows(bootOutQ.b) and cols(ciUpperQ) == cols(bootOutQ.b), "QUAIDS: ciUpper matches bootOutQ.b's shape");
call check(minc(minc(ciUpperQ - ciLowerQ)) >= 0, "QUAIDS: ciUpper >= ciLower elementwise");

rowsBQ = rows(bootOutQ.b);
qLoDirectQ = quantile(bootOutQ.bBoot[., kCheck], alphaCI/2);
iCheckQ = ((kCheck-1) % rowsBQ) + 1;
jCheckQ = floor((kCheck-1)/rowsBQ) + 1;
call check(abs(qLoDirectQ - ciLowerQ[iCheckQ, jCheckQ]) < 1e-10,
    "QUAIDS: ciLower at a specific flattened index matches a direct quantile() call at the correctly-mapped (row, col)");

call printQuaidsCurvatureBootstrap(bootOutQ);
call check(1, "QUAIDS: printQuaidsCurvatureBootstrap() runs without error");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "CURVATURE BOOTSTRAP TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "CURVATURE BOOTSTRAP TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
