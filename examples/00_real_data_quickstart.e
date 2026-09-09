/*
** 00_real_data_quickstart.e
**
** A complete, real-data walkthrough: load a real published dataset from
** a CSV file, run preflight diagnostics, fit the stable LA-AIDS
** baseline, interpret the output, compare against an independent
** reference, compute elasticities, and export a results table. Every
** OTHER example in this suite uses synthetic data chosen to isolate one
** feature at a time -- see example_data.src's own header for why those
** synthetic shares look unrealistic. This example needs no optional
** package (`optmt`/`pubtable`); it runs from a clean `library quaids;`
** installation.
**
** Data: Blanciforti, Green & King (1986) U.S. food-consumption data,
** 1947-1978 (32 annual observations), 4 food categories -- see
** ../tests/fixtures/published/SOURCE.md for the full citation, license
** note, and provenance. The same file backs
** tests/quaids_published_validation_test.e, which independently
** asserts (as part of the automated suite) the exact comparison step 5
** below only prints -- run it yourself with
** `tgauss -b -x quaids_published_validation_test.e` from tests/.
**
** Run from the examples/ directory:
**   tgauss -b -x 00_real_data_quickstart.e
*/

new;
library quaids;

/* ---------------------------------------------------------------------
** 1. Load real data from a CSV and map columns
**
** loadd() reads a CSV into a GAUSS dataframe; selecting columns by name
** avoids assembling matrices by hand and avoids share/price columns
** silently getting out of order. quaidsFit() requires LOG prices and
** LOG total expenditure, not raw levels -- the CSV's price and
** expenditure columns are levels, so ln() them before fitting. This is
** a common, easy-to-miss data-preparation step.
** --------------------------------------------------------------------- */

data = loadd("../tests/fixtures/published/blanciforti86_food32.csv");

goodNames = "Food1"$|"Food2"$|"Food3"$|"Food4";
shareVars = "wFood1"$|"wFood2"$|"wFood3"$|"wFood4";
priceVars = "pFood1"$|"pFood2"$|"pFood3"$|"pFood4";

w = data[., shareVars];
prices = ln(data[., priceVars]);      // CSV holds price LEVELS, not logs
totexp = ln(data[., "xFood"]);        // total food expenditure
instr = ln(data[., "xAgg"]);          // total expenditure across all 11
                                      // commodity groups in the original
                                      // source -- a strong instrument
                                      // for xFood (correlation ~0.97 in
                                      // logs, first-stage R^2 ~0.998)

print "Observations:" rows(w) "years, " cols(w) "food categories:" goodNames';

/* Real published/survey data is rounded to a few decimal places, so raw
   shares essentially never sum to floating-point-exact 1 the way a
   constructed accounting identity does -- this CSV's own rows sum to
   1 +/- up to 0.001, not 1 +/- 1e-8. quaidsPreflight()'s adding-up check
   is a hard failure at a much tighter tolerance (1e-6), by design (a
   real adding-up violation usually means a genuine data error, e.g. a
   missing or duplicated category). Row-normalizing is the standard fix
   for real, rounded data: dividing each share by its own row's sum
   forces exact adding-up while barely perturbing values that were
   already within 0.1% of summing to 1. */
w = w./sumc(w');

/* ---------------------------------------------------------------------
** 2. Preflight: screen the data before spending time fitting
** --------------------------------------------------------------------- */

aCtl = quaidsControlCreate();
aCtl.linear = 1;     // LA-AIDS: the stable baseline -- see README.md's
aCtl.maxiter = 1;    // "Model & Feature Support Tiers" table
aCtl = quaidsSetHomogeneity(aCtl, 1);

pOut = quaidsPreflight(w, 0, prices, totexp, instr, aCtl);
call printQuaidsPreflight(pOut);

if not pOut.ok;
    print "Preflight found a hard failure -- stopping before fitting.";
    end;
endif;

/* ---------------------------------------------------------------------
** 3. Fit
** --------------------------------------------------------------------- */

qOut = quaidsFit(w, 0, prices, totexp, instr, aCtl);

if not qOut.converged;
    print "quaidsFit() did not converge.";
    end;
endif;

call printQuaids(qOut);

/* ---------------------------------------------------------------------
** 4. Interpret: what does this fit actually establish?
**
** qOut.converged==1 here is close to automatic -- aCtl.maxiter=1 (LA-
** AIDS) never iterates, so there is no fixed-point convergence question
** to ask (see README.md's "What qOut.converged proves, and what it does
** not" note). What DOES matter for a real applied result: adding-up/
** homogeneity/symmetry hold by construction (checked below), and the
** symmetry-given-homogeneity test reports whether the DATA supports the
** symmetry restriction this fit imposed, not just the estimator.
** --------------------------------------------------------------------- */

alphaHat = qOut.bS[1, .];
gammaHat = qOut.bS[2:5, .];
betaHat = qOut.bS[6, .];

print "";
print "Adding-up (sum of alpha == 1):" sumc(alphaHat');
print "Homogeneity (gamma row sums == 0):" sumc(gammaHat');
print "Symmetry test given homogeneity: stat =" qOut.symStat ", p-value =" qOut.symPval;
print "  (a small p-value means the data rejects the symmetry restriction";
print "   this fit imposed -- see quaidsJointTest() for a standalone test";
print "   on an unconstrained fit)";

/* ---------------------------------------------------------------------
** 5. Verify your environment: compare against an independent reference
**
** These alpha/beta values are from R's micEconAids::aidsEst(...,
** instNames=...) -- 3SLS with the same log(xAgg) instrument, computed
** independently in R and checked into
** tests/quaids_published_validation_test.e as a real automated
** assertion (not just printed here). GAUSS uses a control-function/
** residual-inclusion IV approach; R uses 3SLS -- both are valid,
** asymptotically consistent estimators for the same model, so they are
** expected to agree closely but not bit-for-bit. If your own
** environment reproduces alpha/beta within the tolerance below, your
** GAUSS installation and this package are both working correctly.
** --------------------------------------------------------------------- */

let alphaR[1, 4] = -0.329388137 0.032896406 0.292524909 1.003966821;
let betaR[1, 4] = 0.373054464 0.100376128 -0.093012355 -0.380418237;
tol = 0.05;

maxDiff = maxc(maxc(abs((alphaHat|betaHat) - (alphaR|betaR))));

print "";
print "alpha (this fit):      " alphaHat;
print "alpha (R reference):   " alphaR;
print "beta  (this fit):      " betaHat;
print "beta  (R reference):   " betaR;
print "Max abs difference (expect <" tol "):" maxDiff;
print "Matches within tolerance?" maxDiff < tol;

/* ---------------------------------------------------------------------
** 6. Elasticities at the sample mean
** --------------------------------------------------------------------- */

n = qOut.n;
nint = qOut.nint;
m_ = meanc(qOut.intcptFull~prices~totexp);
intcptPt = m_[1:1+nint];
pricesPt = m_[1+nint+1:1+nint+n];
totexpPt = m_[1+nint+n+1];

elasOut = quaidsElasFit(qOut.bestB, qOut.bestV, intcptPt, pricesPt, totexpPt, aCtl);
call printQuaidsElas(elasOut);

/* ---------------------------------------------------------------------
** 7. Export a results table
**
** No optional package required: this writes a plain text summary via
** GAUSS's own output redirection. For publication-quality LaTeX/
** Markdown/CSV/HTML/XLSX export instead, see 13_pubtable_reporting.e
** (requires the optional `pubtable` package).
** --------------------------------------------------------------------- */

output file = blanciforti_results.txt reset;
print "Blanciforti86 food demand -- LA-AIDS results";
print "Budget categories:" goodNames';
print "alpha:" alphaHat;
print "beta :" betaHat;
print "gamma:";
print gammaHat;
print "Income elasticities:" elasOut.er';
output off;

print "";
print "Results written to examples/blanciforti_results.txt";
print "";
print "Next: see examples/README.md for the rest of the suite (synthetic";
print "data, chosen to isolate one feature at a time).";
