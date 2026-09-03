/*
** 01_basic_estimation.e
**
** The starting point for the example suite: build a small synthetic
** household budget-share dataset, choose a model, fit it, and read the
** printed report. See docs/USAGE_GUIDE.md's "Choosing An API" and
** "Choosing A Model" sections for the full picture; examples/README.md
** lists the rest of the suite (dataframe input, diagnostics, hypothesis
** tests, elasticities, welfare, zero-share correction, robust/replicate
** standard errors, curvature, survey weighting, workflows, reporting).
**
** Run from the examples/ directory:
**   tgauss -b -x 01_basic_estimation.e
*/

new;
library quaids;
#include example_data.src

/* ---------------------------------------------------------------------
** 1. The data
**
** quaidsExampleData() (examples/example_data.src) simulates a household
** budget survey: 5 spending categories (w, TxN budget shares that sum to
** 1 in every row), an intercept-shifter (intcpt, here a household-size
** variable), log prices for each category (prices), log total
** expenditure (totexp), and an instrument for total expenditure
** (instr, a log wage-income variable -- quaidsFit() always treats total
** expenditure as endogenous and instruments it, so this argument is
** never optional).
** --------------------------------------------------------------------- */

{ w, intcpt, prices, totexp, instr } = quaidsExampleData(3000, 204);
goodNames = quaidsExampleGoodNames();

print "Budget categories:" goodNames';
print "Observations:" rows(w);
print "Row sums of w (adding-up; should all equal 1):" minc(sumc(w')) maxc(sumc(w'));

/* A genuine caveat, not a bug: this synthetic DGP's INDIVIDUAL shares are
   not bounded to [0,1] the way real expenditure fractions are, even
   though they always sum to exactly 1 (checked above) -- see
   example_data.src's own header comment for why. Real survey shares
   don't look like this; this dataset is tuned for reliable estimator
   convergence across the whole example suite, not visual realism. */
print "Mean w (can be outside [0,1] -- see the comment above):" meanc(w)';

/* ---------------------------------------------------------------------
** 2. The control structure
**
** quaidsControlCreate() returns a struct with sensible defaults. The
** two switches that matter most:
**
**   aCtl.linear   1 = AIDS (Stone or iterated-translog price index)
**                 0 = QUAIDS (adds a quadratic log-expenditure term)
**   aCtl.maxiter  1 = one-step LA-AIDS (Stone price index, never iterates)
**                 >1 = iterate the (nonlinear) translog price index
**
** All three named models (LA-AIDS, iterated AIDS, QUAIDS) are really
** just this one estimator under different settings -- see
** docs/USAGE_GUIDE.md#choosing-a-model-la-aids-vs-iterated-aids-vs-quaids
** for the full table. This example fits QUAIDS, since the synthetic data
** genuinely has a quadratic term.
** --------------------------------------------------------------------- */

aCtl = quaidsControlCreate();
aCtl.linear = 0;         // QUAIDS
aCtl.maxiter = 100;      // iterate the translog price index to convergence
aCtl.homogenous = 1;     // impose homogeneity, and test/report symmetry
aCtl.err = .0001;        // relative parameter-change convergence tolerance

/* ---------------------------------------------------------------------
** 3. Fit and read the report
**
** quaidsFit() is silent -- it does 100% of the estimation with no
** console output and returns a quaidsOut struct. printQuaids() then
** reproduces the full estimation-stage report (IV first stage,
** iteration summary, homogeneity-constrained coefficients,
** overidentification test, symmetry test, symmetry-constrained
** coefficients) from that struct alone. Splitting fitting from printing
** this way is what lets a script check qOut.converged before deciding
** whether the result is even worth printing.
** --------------------------------------------------------------------- */

qOut = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

if not qOut.converged;
    print "quaidsFit() did not converge -- see docs/USAGE_GUIDE.md's";
    print "'Choosing A Model' section on aCtl.relax for a mitigation.";
    end;
endif;

call printQuaids(qOut);

print "";
print "Model actually fit:" qOut.model;
print "Iterations to converge:" qOut.iterations;
print "Symmetry rejected at the 5% level?" qOut.symPval < .05;

/* ---------------------------------------------------------------------
** 4. The one-call alternative
**
** quaids() is the original, backward-compatible entry point: it calls
** quaidsFit(), prints the same estimation report via printQuaids(), and
** additionally prints elasticities at four fixed points (mean/quartiles)
** plus descriptive statistics and the Slutzky diagnostic -- all in one
** call, at the cost of not being able to inspect qOut before it prints.
** Prefer quaidsFit()/printQuaids() (above) in scripts; use quaids() for
** a quick, fully-printed one-liner.
** --------------------------------------------------------------------- */

print "";
print "=== The one-call quaids() wrapper (same fit, fuller printed report) ===";
{ b1, v1, b2, v2 } = quaids(w, intcpt, prices, totexp, instr, aCtl);

print "";
print "Next: 02_dataframe_input.e (loading data as a named-column dataframe),";
print "03_preflight_diagnostics.e (screening data before fitting), and the";
print "rest of the suite -- see examples/README.md for the full list.";
