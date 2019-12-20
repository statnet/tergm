#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP godfather_wrapper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCDyn_wrapper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCDynSArun_wrapper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"godfather_wrapper",    (DL_FUNC) &godfather_wrapper,     8},
    {"MCMCDyn_wrapper",      (DL_FUNC) &MCMCDyn_wrapper,      14},
    {"MCMCDynSArun_wrapper", (DL_FUNC) &MCMCDynSArun_wrapper, 19},
    {NULL, NULL, 0}
};

void R_init_tergm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
