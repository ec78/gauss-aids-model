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
** e2_run_m0.e -- E2 confirmation run (CORRECTNESS_VALIDATION_PROTOCOL.md):
** model 0, HELD-OUT seeds 31-200, all 16 starts, all 3 methods -> e2_out/.
** Resumable: rerunning skips every cell already in
** e2_out/cells_m0.csv. Run from this directory:
**   C:\gauss26\tgauss.exe -b -x e2_run_m0.e
*/
pilotModels = 0;
pilotSeeds = seqa(31, 1, 170);
pilotStarts = seqa(1, 1, 16);
pilotMethods = 1|2|3;
pilotOutDir = "e2_out";
pilotVerbose = 0;
#include pilot_harness_core.src;
print "e2_run_m0.e: run complete.";
