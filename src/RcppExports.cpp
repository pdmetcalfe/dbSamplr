// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gsProbs
Rcpp::NumericVector gsProbs(const double p, const Rcpp::IntegerVector& sizes, const Rcpp::IntegerVector& crits);
RcppExport SEXP dbSamplr_gsProbs(SEXP pSEXP, SEXP sizesSEXP, SEXP critsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type crits(critsSEXP);
    rcpp_result_gen = Rcpp::wrap(gsProbs(p, sizes, crits));
    return rcpp_result_gen;
END_RCPP
}
