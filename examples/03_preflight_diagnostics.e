/*
** 03_preflight_diagnostics.e
**
** quaidsPreflight() is a silent, estimator-free screen you can run
** BEFORE quaidsFit(): dimensions, non-finite values, share adding-up,
** zero/negative shares, price/expenditure/instrument variation, IV
** strength, design invertibility, cluster counts, and a basic
** convergence-risk hint. It never estimates anything itself -- it is
** meant to catch common data problems cheaply, before spending time on
** a fit that is likely to fail or mislead. See docs/USAGE_GUIDE.md's
** "Choosing An API" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 03_preflight_diagnostics.e
*/

new;
library quaids;
#include example_data.src

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;

/* ---------------------------------------------------------------------
** 1. A pass with a real, honest warning
**
** quaidsExampleZeroData() (used again in 07_zero_share_correction.e)
** has genuine zero budget shares -- some households report no spending
** at all on some goods, a real corner solution, not an error.
** quaidsPreflight() flags this as a WARNING (pOut.ok stays 1): it is
** worth knowing about (it is exactly what motivates
** quaidsZeroFit()/07_zero_share_correction.e), but it does not by
** itself make the data unusable for a plain quaidsFit() call.
** --------------------------------------------------------------------- */

{ wZero, intcptZero, pricesZero, totexpZero, instrZero } = quaidsExampleZeroData(3000, 1);

print "=== Preflight on a dataset with genuine zero shares ===";
pOut = quaidsPreflight(wZero, intcptZero, pricesZero, totexpZero, instrZero, aCtl);
call printQuaidsPreflight(pOut);

print "pOut.ok:" pOut.ok " (warnings don't block estimation; hard errors do)";
print "zeroShareFrac:" pOut.zeroShareFrac;

/* ---------------------------------------------------------------------
** 2. A second, added warning -- one price held perfectly constant
**
** Real data occasionally has this problem too: a price series that
** never actually varies in your sample carries no identifying
** information for that good's own price effect.
** --------------------------------------------------------------------- */

pricesFlat = pricesZero;
pricesFlat[., 1] = meanc(pricesZero[., 1]) * ones(rows(pricesZero), 1);

print "";
print "=== Same data, with good 1's price also held constant ===";
pOutFlat = quaidsPreflight(wZero, intcptZero, pricesFlat, totexpZero, instrZero, aCtl);
print "nWarnings (was" pOut.nWarnings ", now):" pOutFlat.nWarnings;
print "lowPriceVariation flag:" pOutFlat.lowPriceVariation;
print "minPriceStd (good 1's price std dev is now exactly 0):" pOutFlat.minPriceStd;

/* ---------------------------------------------------------------------
** 3. A hard failure -- negative budget shares
**
** Unlike zero shares or low price variation (warnings), a NEGATIVE
** share cannot represent a real expenditure fraction at all --
** quaidsPreflight() reports pOut.ok == 0 for it, a hard stop.
**
** The standard example dataset (quaidsExampleData(), used in
** 01_basic_estimation.e and most of this suite) actually triggers this
** on its own, with no perturbation needed: see example_data.src's own
** header comment for why (the deterministic price/expenditure swing in
** this synthetic DGP family routinely pushes some cells below 0, a
** known, honestly-documented limitation of simulated data, not
** something that would happen with real survey shares). This is exactly
** the kind of signal quaidsPreflight() exists to catch in real data --
** a negative share there would mean investigating the source data, not
** proceeding to quaidsFit().
** --------------------------------------------------------------------- */

{ w, intcpt, prices, totexp, instr } = quaidsExampleData(3000, 204);

print "";
print "=== Preflight on the standard (01_basic_estimation.e) dataset ===";
pOutStd = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl);
print "pOut.ok:" pOutStd.ok " (0 == hard failure)";
print "negativeShareCount:" pOutStd.negativeShareCount " of" rows(w)*cols(w) "cells";
print "";
print "A script that gates on pOut.ok would stop here rather than call";
print "quaidsFit() -- 01_basic_estimation.e calls quaidsFit() directly on";
print "this same data anyway, since the negative cells are a known,";
print "documented property of this particular synthetic dataset, not a";
print "real data-quality problem to fix.";
