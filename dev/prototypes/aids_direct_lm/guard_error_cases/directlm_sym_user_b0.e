/*
** Expected-failure guard test (symmetry protocol I8): arm S of the private AIDS direct-LM prototype
** (_quaidsDirectSymFit, quaidsdirectsym.src; it inherits _quaidsDirectCheckInputs, so the message keeps the _quaidsDirectFit prefix) must reject a user starting value in aCtl.b0 (never validated).
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src
#include ../src/quaidsiv.src
#include ../src/quaidselas.src
#include ../src/quaidsslutzky.src
#include ../src/quaids.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirect.src;
#include ../dev/prototypes/aids_direct_lm/quaidsdirectsym.src;
#include quaidsfixtures.src;

{ w, intcpt, prices, totexp, instr, tp } = _quaidsSyntheticDGP(400, 12345, 0, 1);

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();
aCtl.linear = 1;
aCtl.homogenous = 1;
aCtl.b0 = ones(8, 5);

dOut = _quaidsDirectSymFit(w, intcpt, prices, totexp, instr, aCtl);
