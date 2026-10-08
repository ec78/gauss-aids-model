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
** e2_determinism.e -- A1 determinism re-run for the E2 confirmation run
** (CORRECTNESS_VALIDATION_PROTOCOL.md). Re-solves a fixed, prespecified
** subset -- seeds 31, 115 and 200, both models, all starts and methods --
** into a SEPARATE directory e2_det_out/, to be compared with e2_out/ by
**   python pilot_analyze.py --compare e2_det_out e2_out
** (every column except wallTime must be bit-identical). Run only after
** e2_run_m0.e / e2_run_m1.e have finished. Run from this directory:
**   C:\gauss26\tgauss.exe -b -x e2_determinism.e
*/
pilotModels = 0|1;
pilotSeeds = 31|115|200;
pilotStarts = seqa(1, 1, 16);
pilotMethods = 1|2|3;
pilotOutDir = "e2_det_out";
pilotVerbose = 0;
#include pilot_harness_core.src;
print "e2_determinism.e: run complete.";
