/*
** quaids_robust_bootstrap_test.e
**
** Milestone 20: validates quaidsRobustBootstrapFit()/
** printQuaidsRobustBootstrap() (src/quaidsrobust.src) -- a cluster-aware
** nonparametric bootstrap alternative to quaidsRobustFit()'s closed-form
** sandwich, resampling whole clusters (or, when clusterId is unset, plain
** i.i.d. rows -- the same G=nobs special case quaidsRobustFit() itself
** uses).
**
** Not run by tests/run_source_tests.ps1's default invocation -- added to
** the same -SkipBootstrap-gated group as quaids_curvature_bootstrap_test.e
** (refitting quaidsFit() B times per case adds real runtime, the same
** reason that file is gated).
**
** Checks, in order: bootstrap run bookkeeping (requested/completed/
** failed/attempts); shape/finiteness of the bootstrap SE; the reshape/
** cell-position regression guard for seBoot (written from day one, the
** same class of bug Milestone 18 found and fixed twice already in
** quaidsCurvatureFit()/quaidsCurvatureBootstrapFit()); exact echo of the
** base (unresampled) point estimate and closed-form seRobust; and the
** core plausibility check -- on tests/quaidsfixtures.src's
** _quaidsClusterSyntheticDGP fixture, a cluster-aware bootstrap's seBoot
** is measurably larger than a plain-row bootstrap's seBoot on the SAME
** genuinely-clustered data, mirroring tests/quaids_robust_test.e's
** identical closed-form finding.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_robust_bootstrap_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/quaidsrobust.src
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

seed = 204;
tobs = 1000;
nClustersTrue = 40;

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, clusterId } = _quaidsClusterSyntheticDGP(tobs, seed, nClustersTrue);

B = 10;
bootSeed = 42;

struct quaidsRobustBootOut rbOut;
rbOut = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, clusterId, B, bootSeed);

call check(rbOut.nRequested == B, "nRequested echoes B");
call check(rbOut.nCompleted + rbOut.nFailed <= rbOut.nAttempts, "completed+failed <= attempts");
call check(rbOut.nCompleted >= 1, "at least one bootstrap replication completed");
call check(rbOut.nClusters == nClustersTrue, "nClusters matches the true number of clusters");

call check(rows(rbOut.seBoot) == rows(rbOut.b) and cols(rbOut.seBoot) == cols(rbOut.b), "seBoot matches b's shape");
call check(sumc(sumc(rbOut.seBoot .== rbOut.seBoot)) == rows(rbOut.seBoot)*cols(rbOut.seBoot), "seBoot contains no NaN/missing values");
call check(minc(minc(rbOut.seBoot)) >= 0, "seBoot are all non-negative");
call check(rows(rbOut.bBoot) == rbOut.nCompleted and cols(rbOut.bBoot) == rows(rbOut.b)*cols(rbOut.b), "bBoot has one row per completed replication and one column per coefficient");

/* Reshape/cell-position regression guard -- written from day one. */
seBootCheck = reshape(stdc(rbOut.bBoot), cols(rbOut.b), rows(rbOut.b))';
call check(maxc(maxc(abs(rbOut.seBoot - seBootCheck))) == 0, "seBoot's cells are correctly positioned relative to bootOut.b (reshape/transpose regression guard)");

/* Exact echo of the base (unresampled) point estimate / closed-form SE. */
struct quaidsOut qOutBase;
qOutBase = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
struct quaidsRobustOut rOutBase;
rOutBase = quaidsRobustFit(qOutBase, w, prices, totexp, aCtl, clusterId);
call check(maxc(maxc(abs(rbOut.b - rOutBase.b))) == 0, "rbOut.b exactly echoes the base (unresampled) reduced-form point estimate");
call check(maxc(maxc(abs(rbOut.seRobust - rOutBase.se))) == 0, "rbOut.seRobust exactly echoes the base closed-form robust/cluster SE");

call printQuaidsRobustBootstrap(rbOut);
call check(1, "printQuaidsRobustBootstrap() runs without error");


/* --- Core plausibility check: a cluster-aware bootstrap's seBoot is
   measurably larger than a plain-row bootstrap's seBoot on the SAME
   genuinely-clustered data -- mirroring quaids_robust_test.e's identical
   closed-form finding. Same seed for both, so the only difference is
   whether whole clusters or plain rows are resampled. --- */

struct quaidsRobustBootOut rbOutNaive;
rbOutNaive = quaidsRobustBootstrapFit(w, intcpt, prices, totexp, instr, aCtl, 0, B, bootSeed);

call check(meanc(meanc(rbOut.seBoot)) > meanc(meanc(rbOutNaive.seBoot)), "cluster-aware bootstrap seBoot is measurably larger than a plain-row bootstrap's seBoot on genuinely clustered data");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "ROBUST/CLUSTER BOOTSTRAP TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "ROBUST/CLUSTER BOOTSTRAP TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
