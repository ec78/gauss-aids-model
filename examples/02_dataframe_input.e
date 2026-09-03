/*
** 02_dataframe_input.e
**
** quaidsFull() selects columns from an already-loaded GAUSS dataframe by
** name instead of requiring the caller to assemble plain matrices by
** hand -- useful when your data already lives in a dataframe (e.g. from
** loadd("mydata.csv")). This example builds a small dataframe from the
** same synthetic data as 01_basic_estimation.e and confirms
** quaidsFull() reproduces quaidsFit()'s exact output on it. See
** docs/USAGE_GUIDE.md's "Choosing An API" section.
**
** Run from the examples/ directory:
**   tgauss -b -x 02_dataframe_input.e
*/

new;
library quaids;
#include example_data.src

{ w, intcpt, prices, totexp, instr } = quaidsExampleData(3000, 204);
goodNames = quaidsExampleGoodNames();

/* ---------------------------------------------------------------------
** 1. Assemble a dataframe with meaningful column names
**
** shareVars and priceVars are matched BY POSITION (shareVars[i] and
** priceVars[i] must be the same good) -- there is no name-matching
** magic linking "Food" to "Food_Price" automatically, so the two string
** arrays below are built in the same good order deliberately.
** --------------------------------------------------------------------- */

priceNames = goodNames $+ "_Price";

data = asDF(w[.,1], goodNames[1]);
i = 2;
do while i <= rows(goodNames);
    data = dfaddcol(data, goodNames[i], w[.,i]);
    i = i + 1;
endo;

i = 1;
do while i <= rows(goodNames);
    data = dfaddcol(data, priceNames[i], prices[.,i]);
    i = i + 1;
endo;

data = dfaddcol(data, "TotalExpenditure", totexp);
data = dfaddcol(data, "WageIncome", instr);
data = dfaddcol(data, "HouseholdSize", intcpt);

print "Dataframe columns:" getcolnames(data)';
print "Dataframe rows:" rows(data);

/* ---------------------------------------------------------------------
** 2. Fit via quaidsFull() -- column names, not column positions
** --------------------------------------------------------------------- */

aCtl = quaidsControlCreate();
aCtl.linear = 0;
aCtl.maxiter = 100;
aCtl.homogenous = 1;
aCtl.err = .0001;

qOutFrame = quaidsFull(data, goodNames, priceNames, "TotalExpenditure",
    "WageIncome", "HouseholdSize", aCtl);

/* ---------------------------------------------------------------------
** 3. Confirm it matches a direct quaidsFit() call on the same data
** --------------------------------------------------------------------- */

qOutMatrix = quaidsFit(w, intcpt, prices, totexp, instr, aCtl);

print "";
print "quaidsFull() converged:" qOutFrame.converged;
print "Max abs difference between quaidsFull() and quaidsFit() coefficients:";
print maxc(maxc(abs(qOutFrame.bestB - qOutMatrix.bestB)));
print "(should be exactly 0 -- quaidsFull() just selects columns and calls quaidsFit() internally)";

/* ---------------------------------------------------------------------
** 4. If your data has no extra demographic shifters, pass 0
** instead of a string array, matching quaidsFit()'s own intcpt == 0
** convention.
** --------------------------------------------------------------------- */

qOutNoExtra = quaidsFull(data, goodNames, priceNames, "TotalExpenditure",
    "WageIncome", 0, aCtl);
print "";
print "quaidsFull() with no extra intercept shifters -- nint:" qOutNoExtra.nint;
