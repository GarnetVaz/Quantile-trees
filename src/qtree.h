#ifndef _QTREE_
#define _QTREE_

#include<Rcpp.h>
#include<armadillo>
#include<vector>

struct ourVector {
  arma::uvec li;
  arma::uvec ri;
  int i;
  double val;
  bool empty;
  double quantile;
  double sold;
};

#include <Rcpp.h>
RcppExport SEXP qtreeCPP(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);


#endif
