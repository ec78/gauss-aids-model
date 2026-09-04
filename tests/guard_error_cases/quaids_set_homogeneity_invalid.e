/*
** Expected-failure guard test: quaidsSetHomogeneity() must reject
** anything other than scalar 0 or 1 (here, a non-binary scalar).
*/

new;
#include ../src/quaids.sdf;
#include ../src/quaidsutil.src

aCtl = quaidsControlCreate();
aCtl = quaidsSetHomogeneity(aCtl, 2);
