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
** pilot_run_m0.e -- Stage B of the matched pilot (PILOT_PROTOCOL.md):
** model 0, seeds 1-30, all 16 starts, all 3 methods -> pilot_out/.
** Resumable: rerunning skips every cell already in
** pilot_out/cells_m0.csv. Run from this directory:
**   C:\gauss26\tgauss.exe -b -x pilot_run_m0.e
*/
pilotModels = 0;
pilotSeeds = seqa(1, 1, 30);
pilotStarts = seqa(1, 1, 16);
pilotMethods = 1|2|3;
pilotOutDir = "pilot_out";
pilotVerbose = 0;
#include pilot_harness_core.src;
print "pilot_run_m0.e: run complete.";
