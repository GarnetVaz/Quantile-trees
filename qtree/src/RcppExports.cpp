// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// qtreeCPP
List qtreeCPP(NumericMatrix pred, NumericVector resp, double minDev, int minCut, int minSize, double tau);
RcppExport SEXP qtree_qtreeCPP(SEXP predSEXP, SEXP respSEXP, SEXP minDevSEXP, SEXP minCutSEXP, SEXP minSizeSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type pred(predSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type resp(respSEXP );
        Rcpp::traits::input_parameter< double >::type minDev(minDevSEXP );
        Rcpp::traits::input_parameter< int >::type minCut(minCutSEXP );
        Rcpp::traits::input_parameter< int >::type minSize(minSizeSEXP );
        Rcpp::traits::input_parameter< double >::type tau(tauSEXP );
        List __result = qtreeCPP(pred, resp, minDev, minCut, minSize, tau);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
