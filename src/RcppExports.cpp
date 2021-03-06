// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// PBasis
arma::mat PBasis(const arma::mat& data, const arma::umat& basis, const bool& intercept, const std::size_t& n_terms, const arma::uvec& reccur_limit);
RcppExport SEXP OptionPrice_PBasis(SEXP dataSEXP, SEXP basisSEXP, SEXP interceptSEXP, SEXP n_termsSEXP, SEXP reccur_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type n_terms(n_termsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type reccur_limit(reccur_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(PBasis(data, basis, intercept, n_terms, reccur_limit));
    return rcpp_result_gen;
END_RCPP
}
// LBasis
arma::mat LBasis(const arma::mat& data, const arma::umat& basis, const bool& intercept, const std::size_t& n_terms, const arma::uvec& reccur_limit);
RcppExport SEXP OptionPrice_LBasis(SEXP dataSEXP, SEXP basisSEXP, SEXP interceptSEXP, SEXP n_termsSEXP, SEXP reccur_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type n_terms(n_termsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type reccur_limit(reccur_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(LBasis(data, basis, intercept, n_terms, reccur_limit));
    return rcpp_result_gen;
END_RCPP
}
// BermudaPutPolicy
arma::umat BermudaPutPolicy(const arma::cube& path, const arma::mat& expected, const double& strike, const double& discount, const arma::umat& basis, const std::string& basis_type);
RcppExport SEXP OptionPrice_BermudaPutPolicy(SEXP pathSEXP, SEXP expectedSEXP, SEXP strikeSEXP, SEXP discountSEXP, SEXP basisSEXP, SEXP basis_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type path(pathSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type expected(expectedSEXP);
    Rcpp::traits::input_parameter< const double& >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< const double& >::type discount(discountSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type basis_type(basis_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(BermudaPutPolicy(path, expected, strike, discount, basis, basis_type));
    return rcpp_result_gen;
END_RCPP
}
// BermudaPutAddDual
arma::mat BermudaPutAddDual(const arma::cube& path, Rcpp::NumericVector subsim_, const arma::mat& expected, const double& strike, const double& discount, const arma::umat& basis, const std::string& basis_type);
RcppExport SEXP OptionPrice_BermudaPutAddDual(SEXP pathSEXP, SEXP subsim_SEXP, SEXP expectedSEXP, SEXP strikeSEXP, SEXP discountSEXP, SEXP basisSEXP, SEXP basis_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type path(pathSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type subsim_(subsim_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type expected(expectedSEXP);
    Rcpp::traits::input_parameter< const double& >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< const double& >::type discount(discountSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type basis_type(basis_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(BermudaPutAddDual(path, subsim_, expected, strike, discount, basis, basis_type));
    return rcpp_result_gen;
END_RCPP
}
// BermudaPutBounds
Rcpp::List BermudaPutBounds(const arma::cube& path, Rcpp::NumericVector subsim_, const arma::mat& expected, const double& strike, const double& discount, const arma::umat& basis, const std::string& basis_type);
RcppExport SEXP OptionPrice_BermudaPutBounds(SEXP pathSEXP, SEXP subsim_SEXP, SEXP expectedSEXP, SEXP strikeSEXP, SEXP discountSEXP, SEXP basisSEXP, SEXP basis_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type path(pathSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type subsim_(subsim_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type expected(expectedSEXP);
    Rcpp::traits::input_parameter< const double& >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< const double& >::type discount(discountSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type basis_type(basis_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(BermudaPutBounds(path, subsim_, expected, strike, discount, basis, basis_type));
    return rcpp_result_gen;
END_RCPP
}
// BermudaPutLSM
Rcpp::List BermudaPutLSM(const arma::cube& path, const double& strike, const double& discount, const arma::umat& basis, const bool& intercept, const std::string& basis_type);
RcppExport SEXP OptionPrice_BermudaPutLSM(SEXP pathSEXP, SEXP strikeSEXP, SEXP discountSEXP, SEXP basisSEXP, SEXP interceptSEXP, SEXP basis_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type path(pathSEXP);
    Rcpp::traits::input_parameter< const double& >::type strike(strikeSEXP);
    Rcpp::traits::input_parameter< const double& >::type discount(discountSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type basis(basisSEXP);
    Rcpp::traits::input_parameter< const bool& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type basis_type(basis_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(BermudaPutLSM(path, strike, discount, basis, intercept, basis_type));
    return rcpp_result_gen;
END_RCPP
}
