new;
library optmt;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.sdf;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsutil.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsiv.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidselas.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaidsslutzky.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/src/quaids.src;
#include quaidsfit_anderson_prototype.src;
#include c:/Users/eclow/Documents/GitHub/gauss-aids-model/tests/quaidsfixtures.src;
#include direct_objective.src;
#include pilot_lib.src;

/*
** pilot_resume_test.e -- V6 (resume) fixture for pilot_resume_test.ps1.
** A tiny job: AIDS seed 1, starts 1-2, all three methods (6 cells) into
** pilot_resume_test/. The driver runs it, reruns it (must add 0 rows),
** deletes the last 2 rows, and reruns it (must recompute exactly those
** 2 cells, identical apart from wall time).
*/
pilotModels = 0;
pilotSeeds = 1;
pilotStarts = 1|2;
pilotMethods = 1|2|3;
pilotOutDir = "pilot_resume_test";
pilotVerbose = 0;
#include pilot_harness_core.src;
