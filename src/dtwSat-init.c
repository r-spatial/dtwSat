// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>
#include "f.h"

static R_FortranMethodDef fortranMethods[] = {
  {"bestmatches",  (DL_FUNC) &F77_SUB(bestmatches), 11},
  {"computecost",  (DL_FUNC) &F77_SUB(computecost),  7},
  {"g",            (DL_FUNC) &F77_SUB(g),            4},
  {"tracepath",    (DL_FUNC) &F77_SUB(tracepath),   11},
  {NULL, NULL, 0}
};

void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
  __attribute__ ((visibility ("default")))
#endif
R_init_dtwSat(DllInfo *dll) {
	R_registerRoutines(dll, NULL, NULL, fortranMethods, NULL);
	R_useDynamicSymbols(dll, FALSE);
}


