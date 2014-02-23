#ifndef _QTREE_
#define _QTREE_

#include<vector>
using namespace std;

struct nodeStruct {
  vector<unsigned int> li;
  vector<unsigned int> ri;
  int i;
  double val;
  bool empty;
  double quantile;
  double sold;
};


#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP qtreeCPP(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
/* List qtreeCPP(NumericMatrix, NumericVector, double, int, int, double); */


#endif
