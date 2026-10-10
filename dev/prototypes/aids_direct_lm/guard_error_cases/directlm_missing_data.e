/*
** Expected-failure guard test: the private AIDS direct-LM prototype
** (_quaidsDirectFit, dev/prototypes/aids_direct_lm/quaidsdirect.src) must reject missing values in the data.
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;
#include quaidsfixtures.src;

{ w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(400, 12345, 0, 1);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;
w[7, 2] = miss(0, 0);

dOut = _quaidsDirectFit(w, intcpt, prices, totexp, instr, aCtl);
