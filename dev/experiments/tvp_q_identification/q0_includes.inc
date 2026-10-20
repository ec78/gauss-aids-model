/*
** q0_includes.inc -- source chain shared by every Q0 job (run from this
** directory). Phase C's conv_includes.inc list, plus Phase C's
** conv_lib.src (reused unchanged) and this study's q0_lib.src.
*/
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
#include ../../../tests/quaidsfixtures.src
#include ../../../examples/example_data.src
#include ../tvp_conv_study/conv_lib.src
#include q0_lib.src
