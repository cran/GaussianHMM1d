#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
extern void est_hmm_1d_new(void *, void *, void *, void *,void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
static const R_CMethodDef CEntries[] = {
  {"est_hmm_1d_new",    (DL_FUNC) &est_hmm_1d_new,    14},
  {NULL, NULL, 0}
};

void R_init_GaussianHMM1d(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
