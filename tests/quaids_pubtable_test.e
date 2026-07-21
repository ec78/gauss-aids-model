/*
** quaids_pubtable_test.e
**
** Milestone 6: validates the pubtable adapter (src/pubtable_quaids.src) --
** ptModelFromQuaids/ptFromQuaids (coefficient tables), ptModelFromQuaidsElas/
** ptFromQuaidsElas/ptTablesFromQuaidsElas (elasticity tables), and the
** ptFromQuaidsFamily dispatcher -- against a real quaidsFit()/
** quaidsElasFit() result.
**
** Checks exact numeric parity between ptModel.estimates/stdErrors and the
** qOut.bestB/qOut.bestV (or elasOut.er/ser) values they're built from --
** not just "it runs" -- plus row/column shape and row-name/title checks,
** and an end-to-end export smoke test (LaTeX/Markdown/CSV) confirming the
** exported files exist and contain the expected content.
**
** Requires the pubtable package installed (this machine has it at
** c:\gauss26\pkgs\pubtable). Run from the tests/ directory:
**   tgauss -b -x quaids_pubtable_test.e
*/

new;
library pubtable;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../src/pubtable_quaids.src;
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

/* Reads a whole text file back as a string array (one row per line), for
   substring checks against exported tables -- checked line-by-line with
   strindx() rather than joined into one string, since strjoin() does not
   collapse a string array with an empty quote-char argument. */
proc (1) = readWholeFile(fname);
    local fh, sa, cls;
    fh = fopen(fname, "r");
    sa = fgetsa(fh, 100000);
    cls = close(fh);
    retp(sa);
endp;

/* True if any line of a string array (as returned by readWholeFile())
   contains needle. */
proc (1) = fileContains(sa, needle);
    retp(maxc(strindx(sa, needle, 1)) > 0);
endp;


struct quaidsControl aCtl;
aCtl = quaidsControlCreate;
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

{ w, intcpt, prices, totexp, instr, trueParams } = _quaidsSyntheticDGP(3000, 204, 1, 1);
struct quaidsOut qOut;
qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

n = qOut.n;
k = rows(qOut.bestB);


/* --- ptModelFromQuaids: one equation's coefficient column. --- */

struct ptModel mdl1;
mdl1 = ptModelFromQuaids("Good 1", qOut, 1);
call check(maxc(abs(mdl1.estimates - qOut.bestB[., 1])) == 0,
    "ptModelFromQuaids: estimates match qOut.bestB column exactly (good 1)");
call check(maxc(abs(mdl1.stdErrors - sqrt(diag(qOut.bestV[1:k, 1:k])))) == 0,
    "ptModelFromQuaids: stdErrors match qOut.bestV block exactly (good 1)");
call check(rows(mdl1.termNames) == k, "ptModelFromQuaids: termNames has K rows");
call check(mdl1.termNames[1] $== "CONSTANT", "ptModelFromQuaids: first term name is CONSTANT");
call check(mdl1.termNames[k] $== qOut.unam[qOut.nu], "ptModelFromQuaids: last term name is the final IV-residual name");

struct ptModel mdl3;
idx3 = seqa((3-1)*k+1, 1, k);
mdl3 = ptModelFromQuaids("Good 3", qOut, 3);
call check(maxc(abs(mdl3.estimates - qOut.bestB[., 3])) == 0,
    "ptModelFromQuaids: estimates match qOut.bestB column exactly (good 3)");
call check(maxc(abs(mdl3.stdErrors - sqrt(diag(qOut.bestV[idx3, idx3])))) == 0,
    "ptModelFromQuaids: stdErrors match qOut.bestV block exactly (good 3)");
call check(maxc(abs(mdl3.estimates - mdl1.estimates)) > 1e-6,
    "ptModelFromQuaids: different goods give genuinely different columns");


/* --- ptFromQuaids: N-column comparison table. --- */

struct ptTable coefTbl;
coefTbl = ptFromQuaids(qOut);
call check(cols(coefTbl.body) == n, "ptFromQuaids: one column per good");
call check(rows(coefTbl.rowNames) == rows(coefTbl.body), "ptFromQuaids: rowNames match body row count");
call check(coefTbl.title $== qOut.model $+ " results", "ptFromQuaids: title reflects fitted model");


/* --- Elasticities at the sample mean (same construction as
   tests/quaids_elasticities_test.e). --- */

nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp~w~instr);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

struct quaidsElasOut elasOut;
elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);

struct ptModel elasMdl;
elasMdl = ptModelFromQuaidsElas("Income elasticities", elasOut);
call check(maxc(abs(elasMdl.estimates - elasOut.er)) == 0,
    "ptModelFromQuaidsElas: estimates match elasOut.er exactly");
call check(maxc(abs(elasMdl.stdErrors - elasOut.ser)) == 0,
    "ptModelFromQuaidsElas: stdErrors match elasOut.ser exactly");
call check(rows(elasMdl.termNames) == n, "ptModelFromQuaidsElas: one term name per good");

struct ptTable elasIncomeTbl;
elasIncomeTbl = ptFromQuaidsElas(elasOut);
call check(elasIncomeTbl.title $== "Income elasticities", "ptFromQuaidsElas: title is Income elasticities");
call check(cols(elasIncomeTbl.body) == 1, "ptFromQuaidsElas: single column");

struct ptTable elasTbls;
elasTbls = ptTablesFromQuaidsElas(elasOut);
call check(rows(elasTbls) == 3, "ptTablesFromQuaidsElas: returns 3 tables");
call check(elasTbls[1].title $== "Income elasticities", "ptTablesFromQuaidsElas[1]: income elasticities");
call check(strindx(elasTbls[2].title, "Uncompensated", 1) > 0, "ptTablesFromQuaidsElas[2]: uncompensated price elasticities");
call check(strindx(elasTbls[3].title, "Compensated", 1) > 0, "ptTablesFromQuaidsElas[3]: compensated price elasticities");
call check(rows(elasTbls[2].body) == 2*n and cols(elasTbls[2].body) == n,
    "ptTablesFromQuaidsElas[2]: n x n elasticity matrix rendered as 2n x n (value+SE rows)");
call check(rows(elasTbls[3].body) == 2*n and cols(elasTbls[3].body) == n,
    "ptTablesFromQuaidsElas[3]: n x n elasticity matrix rendered as 2n x n (value+SE rows)");

/* Numeric parity: the rendered uncompensated-elasticity value cell for
   (good 1, good 1) should round to elasOut.ep[1,1] at the table's own
   formatting precision (3 decimals). */
call check(abs(stof(elasTbls[2].body[1, 1]) - elasOut.ep[1, 1]) < 5e-4,
    "ptTablesFromQuaidsElas[2]: rendered (1,1) cell matches elasOut.ep[1,1] to formatting precision");


/* --- ptFromQuaidsFamily dispatcher. --- */

struct ptTable dispQ;
dispQ = ptFromQuaidsFamily(qOut);
call check(dispQ.title $== coefTbl.title, "ptFromQuaidsFamily(qOut): matches ptFromQuaids");

struct ptTable dispE;
dispE = ptFromQuaidsFamily(elasOut);
call check(dispE.title $== elasIncomeTbl.title, "ptFromQuaidsFamily(elasOut): matches ptFromQuaidsElas");


/* --- End-to-end export smoke test: LaTeX/Markdown/CSV. --- */

call ptExport(coefTbl, "pubtable_test_coef.tex");
call ptExport(coefTbl, "pubtable_test_coef.md");
call ptExport(coefTbl, "pubtable_test_coef.csv");

texText = readWholeFile("pubtable_test_coef.tex");
mdText = readWholeFile("pubtable_test_coef.md");
csvText = readWholeFile("pubtable_test_coef.csv");

call check(fileContains(texText, "\\begin{tabular}"), "export .tex: contains a LaTeX tabular environment");
call check(fileContains(texText, "GAMMA"), "export .tex: contains a GAMMA_ coefficient row");
call check(fileContains(mdText, "|"), "export .md: contains Markdown table pipes");
call check(fileContains(csvText, ","), "export .csv: contains comma-separated fields");
call check(fileContains(csvText, "CONSTANT"), "export .csv: contains the CONSTANT row label");


print;
print "-----------------------------------------------------------";
if nfail == 0;
    print ftos(ncheck, "PUBTABLE ADAPTER TEST: ALL %*.*lf CHECKS PASSED", 1, 0);
else;
    print ftos(nfail, "PUBTABLE ADAPTER TEST: %*.*lf CHECKS FAILED", 1, 0);;
    print ftos(ncheck, " (of %*.*lf total)", 1, 0);
endif;
print "-----------------------------------------------------------";
