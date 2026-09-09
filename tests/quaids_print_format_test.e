/*
** quaids_print_format_test.e
**
** Regression guard for a real bug found while building
** examples/00_real_data_quickstart.e (public-release Phase 4, PR-402):
** printQuaidsElas() (src/quaidselas.src) called `format /rd 7, 0;` (a
** 0-decimal display format) to control a header row's column width, but
** never reset the format before returning -- this is GAUSS's own GLOBAL
** print state, not scoped to the proc, so it silently leaked into every
** subsequent plain `print` statement in the CALLING script. A real fitted
** coefficient like -0.3465 printed as "0" afterward; income elasticities
** like 2.267/1.507/0.346/-0.079 printed as "2"/"2"/"0"/"0" -- values that
** look like small integers/booleans, not obviously wrong at a glance,
** which is exactly why no prior test caught it: every existing `check()`
** in this suite compares the underlying stored numeric VALUE (never
** affected by print format), never the printed console TEXT.
**
** The same `format /rd ...; ... (no reset)` pattern existed in 7 other
** printer procs across the codebase (quaidsshares.src, quaidsrobust.src
** x2, quaidsreplicate.src, quaidszerocorrect.src, quaidscurvature.src x2
** -- the last two need the optional `optmt` package to produce a
** quaidsCurvOut/quaidsCurvBootOut to print in the first place, so they
** are not covered here; they share the identical fix, applied the same
** way). All were fixed by adding `format /rd 16, 8;` (GAUSS's own
** observed default numeric format) immediately before each proc's own
** `endp;`. This file directly regression-tests that fix for two of them
** (printQuaidsElas, the original and most severe case since it used
** 0 decimals; printQuaidsShares, a second file/proc, to confirm the fix
** pattern generalizes) using the technique this repo's own
** quaids_pubtable_test.e already established: write real output to a
** file via GAUSS's own `output file=...` redirection, read it back as
** text, and check the printed digits directly -- the only way to test a
** print-FORMAT bug, since the underlying stored values were never wrong.
**
** This file also regression-tests a second, unrelated real bug found in
** the same session, by the same real-data example: printQuaids()
** (src/quaids.src) crashed with `error G0058: Index out of range` on any
** fit with `nint==0` (no extra intercept-shifter columns -- e.g. real
** published data like Blanciforti86 with no demographic variables),
** because `qOut.ivCor[1:qOut.nint, i]`/`qOut.homogCor[1:qOut.nint, i]`
** used `1:0` as a range, which GAUSS treats as invalid rather than
** empty. Fixed with an explicit `qOut.nint > 0` guard in both places.
**
** Run from the tests/ directory:
**   tgauss -b -x quaids_print_format_test.e
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsshares.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
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

/* Same read-back-a-text-file technique as tests/quaids_pubtable_test.e's
   own readWholeFile()/fileContains() helpers. */
proc (1) = readWholeFile(fname);
    local fh, sa, cls;
    fh = fopen(fname, "r");
    sa = fgetsa(fh, 100000);
    cls = close(fh);
    retp(sa);
endp;

proc (1) = fileContains(sa, needle);
    retp(maxc(strindx(sa, needle, 1)) > 0);
endp;

/* A distinctive marker value with digits in every decimal place, so
   rounding to fewer decimals (or to a bare integer) is unambiguous:
   a leaked 0-decimal format would print "3", a leaked 4-decimal format
   would print "3.1416" -- neither contains the full "14159265" mantissa
   checked for below. */
markerValue = 3.14159265;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(500, 204, 0, 1);

aCtl = quaidsControlCreate;
aCtl.linear = 1;
aCtl.maxiter = 1;
aCtl.homogenous = 1;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);
call check(qOut.converged == 1, "base fit converged");

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptPt = m_[1:1+nint];
pricesPt = m_[1+nint+1:1+nint+n];
totexpPt = m_[1+nint+n+1];

/* --- printQuaidsElas(): the original, most severe (0-decimal) leak --- */

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);

output file = print_format_probe_elas.txt reset;
call printQuaidsElas(elasOut);
print "MARKER:" markerValue;
output off;

sa = readWholeFile("print_format_probe_elas.txt");
call check(fileContains(sa, "14159265"), "printQuaidsElas() does not leak its 0-decimal format to a later print statement");

/* --- printQuaidsShares(): a second file/proc, confirming the fix generalizes --- */

sharesOut = quaidsSharesFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);

output file = print_format_probe_shares.txt reset;
call printQuaidsShares(sharesOut);
print "MARKER:" markerValue;
output off;

sa = readWholeFile("print_format_probe_shares.txt");
call check(fileContains(sa, "14159265"), "printQuaidsShares() does not leak its format to a later print statement");

/* --- printQuaids() with nint==0: the G0058 crash fix --- */

aCtlNoInt = quaidsControlCreate;
aCtlNoInt.linear = 1;
aCtlNoInt.maxiter = 1;
aCtlNoInt.homogenous = 1;

qOutNoInt = quaidsFit(w, 0, prices, totexp, instr, aCtlNoInt);
call check(qOutNoInt.nint == 0, "no-intercept-shifter fit has nint==0");

output file = print_format_probe_noint.txt reset;
call printQuaids(qOutNoInt);
output off;
call check(1, "printQuaids() with nint==0 completed without crashing (G0058-class index-range regression guard)");

print "";
if nfail == 0;
    print "PRINT FORMAT TEST: ALL " $+ ntos(ncheck) $+ " CHECKS PASSED";
else;
    print "PRINT FORMAT TEST: " $+ ntos(nfail) $+ " OF " $+ ntos(ncheck) $+ " CHECKS FAILED";
endif;
