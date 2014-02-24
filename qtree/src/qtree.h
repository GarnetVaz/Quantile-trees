#ifndef _QTREE_
#define _QTREE_

#include<vector>
using namespace std;

#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP qtreeCPP(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
/* List qtreeCPP(NumericMatrix, NumericVector, double, int, int, double); */


#endif
