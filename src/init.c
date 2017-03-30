#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(bestmatches)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(computecost)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(g)(void *, void *, void *, void *);
extern void F77_NAME(tracepath)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"bestmatches", (DL_FUNC) &F77_NAME(bestmatches), 11},
  {"computecost", (DL_FUNC) &F77_NAME(computecost),  7},
  {"g",           (DL_FUNC) &F77_NAME(g),            4},
  {"tracepath",   (DL_FUNC) &F77_NAME(tracepath),   11},
  {NULL, NULL, 0}
};

void R_init_dtwSat(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);
}


