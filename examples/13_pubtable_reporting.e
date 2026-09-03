/*
** 13_pubtable_reporting.e
**
** Publication-quality LaTeX/Markdown/CSV/RTF/HTML/XLSX table export via
** the optional pubtable adapter (src/pubtable_quaids.src). This adapter
** is NOT part of the installed quaids package (it has a hard compile-
** time dependency on pubtable's own struct types), so it is #included
** from the source tree directly, after both quaids and pubtable are
** loaded -- see docs/USAGE_GUIDE.md's "Reporting (pubtable)" section.
**
** Requires the pubtable package installed. Run from the examples/
** directory (needed for the ../src/ relative #include below):
**   tgauss -b -x 13_pubtable_reporting.e
*/

new;
library pubtable, quaids;
#include quaids.sdf
#include ../src/pubtable_quaids.src
#include example_data.src

{ w, intcpt, prices, totexp, instr } = quaidsExampleData(3000, 204);

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

/* ---------------------------------------------------------------------
** 1. A coefficient table straight from a fitted quaidsOut
** --------------------------------------------------------------------- */

coefTbl = ptFromQuaids(qOut);
call ptExport(coefTbl, "quaids_coefficients.tex");
call ptExport(coefTbl, "quaids_coefficients.md");
call ptExport(coefTbl, "quaids_coefficients.csv");

/* ---------------------------------------------------------------------
** 2. Elasticity tables at the sample mean
** --------------------------------------------------------------------- */

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptMean = m_[1:1+nint];
pricesMean = m_[1+nint+1:1+nint+n];
totexpMean = m_[1+nint+n+1];

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptMean, pricesMean, totexpMean, aCtl);
elasTbls = ptTablesFromQuaidsElas(elasOut);   // 3x1: income, uncompensated, compensated

call ptExport(elasTbls[1], "quaids_income_elasticities.md");
call ptExport(elasTbls[2], "quaids_uncompensated_elasticities.tex");
call ptExport(elasTbls[3], "quaids_compensated_elasticities.csv");

/* ---------------------------------------------------------------------
** 3. A full applied-workflow table bundle in one call
** --------------------------------------------------------------------- */

wfOut = quaidsWorkflowFit(w, intcpt, prices, totexp, instr, aCtl);
workflowTbls = ptTablesFromQuaidsWorkflow(wfOut);
call ptExportAll(workflowTbls, "quaids_workflow");

print "Exported quaids_coefficients.{tex,md,csv},";
print "quaids_income_elasticities.md, quaids_uncompensated_elasticities.tex,";
print "quaids_compensated_elasticities.csv, and quaids_workflow* table files";
print "to the current working directory.";
