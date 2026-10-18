/*
** p6_scale_fit.template.e -- TVP-AIDS P6 scale re-measurement template.
** run_p6.ps1 substitutes 5, 500, 5 and runs the
** generated copy from this directory (GAUSS26_CFG override + outer
** hard timeout). See RESULTS.md.
*/
new;
library cmlmt, tsmt, sslib;
#include ../../../src/quaids.sdf
#include ../../../src/quaidsutil.src
#include ../../../src/quaidsiv.src
#include ../../../src/quaidselas.src
#include ../../../src/quaidsslutzky.src
#include ../../../src/quaids.src
#include ../../../src/quaidstvp.src
#include ../../../src/quaidstvpkalman.src
#include ../../../src/quaidstvpmle.src
#include ../../../src/quaidstvpsmooth.src
#include ../../../src/quaidstvpelas.src
#include ../../../src/quaidstvpfit.src
#include ../../../examples/example_data.src

nSub = 5;
tobs = 500;
seed = 204;
capIters = 500;
capTime = 5;

#include p6_fit_body.inc

