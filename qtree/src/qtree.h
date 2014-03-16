#ifndef _QTREE_
#define _QTREE_



struct ourVector {
  arma::uvec li;
  arma::uvec ri;
  int i;
  double val;
  bool empty;
  double quantile;
  double sold;
};

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List qtreeCPP(const arma::mat &X, const arma::vec &Y,
		    double mindev, int mincut, int minsize, double tau)

#endif
