#ifndef VCRPART_OLMM_H
#define VCRPART_OLMM_H
/* #define USE_FC_LEN_T /* added 2022-02-03, see https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings */

#include <R.h>
#include <Rdefines.h>

double olmm_gLink(double x, int link);
double olmm_GLink(double x, int link);

SEXP olmm_setPar(SEXP x, SEXP par);
SEXP olmm_update_marg(SEXP x, SEXP par);
SEXP olmm_update_u(SEXP x);
SEXP olmm_pred_marg(SEXP x, SEXP eta, SEXP W, SEXP n, SEXP pred);
SEXP olmm_pred_margNew(SEXP x, SEXP etaNew, SEXP WNew, SEXP subjectNew, 
		       SEXP nNew, SEXP pred);
#endif
