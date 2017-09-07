// Copyright 2017 <jeremyyee@outlook.com.au>
// Lower and upper bounds for the true value
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "inst/include/basis.h"

// Extracting the prescribed policy on a set of sample paths
//[[Rcpp::export]]
arma::umat BermudaPutPolicy(const arma::cube& path,
                            const arma::mat& expected,
                            const double& strike,
                            const double& discount,
                            const arma::umat& basis,
                            const std::string& basis_type) {
  // Extract parameters
  const std::size_t n_dec = path.n_slices;
  const std::size_t n_path = path.n_rows;
  const std::size_t n_dim = path.n_cols;
  const std::size_t n_terms = expected.n_rows;
  bool intercept = true;
  if (n_terms == arma::accu(basis)) {
    intercept = false;
  }
  // Extract information about regression basis
  arma::mat reg_basis(n_path, n_terms);
  arma::uvec reccur_limit(basis.n_rows);
  reccur_limit = ReccurLimit(basis);
  // Extract the prescribed policy
  arma::umat policy(n_path, n_dec - 1, arma::fill::zeros);
  arma::vec fitted(n_path);
  arma::mat states(n_path, n_dim);
  for (std::size_t tt = 0; tt < n_dec - 1; tt++) {
    states = path.slice(tt);
    if (basis_type == "power") {
      reg_basis = PBasis(states, basis, intercept, n_terms, reccur_limit);
    } else if (basis_type == "laguerre") {
      reg_basis = LBasis(states, basis, intercept, n_terms, reccur_limit);
    }
    fitted = reg_basis * expected.col(tt);  // fitted continuation values
    for (std::size_t pp = 0; pp < n_path; pp++) {
      if ((strike - states(pp, 0)) > std::max(fitted(pp), 0.)) {
        policy(pp, tt) = 1;  // 1 = exercise, 0 = continue
      }
    }
  }
  return policy;
}

// Compute the additive duals
//[[Rcpp::export]]
arma::mat BermudaPutAddDual(const arma::cube& path,
                            Rcpp::NumericVector subsim_,
                            const arma::mat& expected,
                            const double& strike,
                            const double& discount,
                            const arma::umat& basis,
                            const std::string& basis_type) {
  // Extract parameters
  const std::size_t n_dec = path.n_slices;
  const std::size_t n_path = path.n_rows;
  const std::size_t n_dim = path.n_cols;
  const std::size_t n_terms = expected.n_rows;
  bool intercept = true;
  if (n_terms == arma::accu(basis)) {
    intercept = false;
  }
  arma::uvec reccur_limit(basis.n_rows);
  reccur_limit = ReccurLimit(basis);
  const arma::ivec s_dims = subsim_.attr("dim");
  const std::size_t n_subsim = s_dims(0);
  arma::cube subsim(subsim_.begin(), n_subsim, n_dim, n_path * (n_dec - 1), false);
  // Additive duals
  arma::mat add_dual(n_path, n_dec - 1, arma::fill::zeros);
  arma::mat path_basis(n_path, n_terms);
  arma::mat subsim_basis(n_subsim, n_terms);
  arma::vec fitted(n_path);
  arma::vec subsim_fitted(n_subsim);
  arma::vec temp1(n_subsim);
  arma::mat states(n_subsim, n_path);
  std::size_t tt, pp;
  arma::mat compare1(n_subsim, 3, arma::fill::zeros);
  arma::mat compare2(n_path, 3, arma::fill::zeros);
  for (tt = 0; tt < n_dec - 2; tt++) {
    // Find the average of the subsimulation paths
    for (pp = 0; pp < n_path; pp++) {
      states = subsim.slice(n_path * tt + pp);
      if (basis_type == "power") {
        subsim_basis = PBasis(states, basis, intercept, n_terms, reccur_limit);
      } else if (basis_type == "laguerre") {
        subsim_basis = LBasis(states, basis, intercept, n_terms, reccur_limit);
      }
      // Fitted expected value function for subsims
      compare1.col(0) = strike - states.col(0);
      compare1.col(2) = subsim_basis * expected.col(tt + 1);
      // Reward for subsims
      temp1 = arma::max(compare1, 1);
      add_dual(pp, tt) += arma::sum(temp1);
    }
    add_dual.col(tt) = (1.0 / n_subsim) * add_dual.col(tt);
    // Find the realised values. Reg basis for paths.
    if (basis_type == "power") {
      path_basis = PBasis(path.slice(tt + 1), basis, intercept, n_terms, reccur_limit);
    } else if (basis_type == "laguerre") {
      path_basis = LBasis(path.slice(tt + 1), basis, intercept, n_terms, reccur_limit);
    }
    // Fitted expected value function for realised paths
    compare2.col(0) = strike - path.slice(tt + 1).col(0);
    compare2.col(2) = path_basis * expected.col(tt + 1);
    add_dual.col(tt) -= arma::max(compare2, 1);
  }
  // Find the duals for the last time epoch
  // Find the average of the subsimulation paths
  tt = n_dec - 2;
  for (pp = 0; pp < n_path; pp++) {
    // Average for each path
    compare1.col(0) = strike - subsim.slice(n_path * tt + pp).col(0);
    temp1 = arma::max(compare1.cols(0, 1), 1);
    add_dual(pp, tt) += arma::sum(temp1);
  }
  add_dual.col(tt) = (1.0 / n_subsim) * add_dual.col(tt);
  // Find the realised value
  compare2.col(0) = strike - path.slice(tt + 1).col(0);
  add_dual.col(tt) -= arma::max(compare2.cols(0, 1), 1);  
  return add_dual;
}

// Lower and upper bounds for the true value
//[[Rcpp::export]]
Rcpp::List BermudaPutBounds(const arma::cube& path,
                            Rcpp::NumericVector subsim_,
                            const arma::mat& expected,
                            const double& strike,
                            const double& discount,
                            const arma::umat& basis,
                            const std::string& basis_type) {
  // Extract parameters
  const std::size_t n_dec = path.n_slices;
  const std::size_t n_path = path.n_rows;
  const std::size_t n_dim = path.n_cols;
  const std::size_t n_terms = expected.n_rows;
  bool intercept = true;
  if (n_terms == arma::accu(basis)) {
    intercept = false;
  }
  const arma::umat path_policy =
      BermudaPutPolicy(path, expected, strike, discount, basis, basis_type);
  arma::mat mart =
      BermudaPutAddDual(path, subsim_, expected, strike, discount, basis, basis_type);
  // Initialise with last time
  arma::mat compare(n_path, 3, arma::fill::zeros);
  arma::mat states(n_path, n_dim);
  states = path.slice(n_dec - 1);
  arma::vec primals(n_path);
  compare.col(0) = strike - states.col(0);
  primals = arma::max(compare.cols(0, 1), 1);
  arma::vec duals = primals;
  // Perform the backward induction
  arma::uword policy;
  arma::vec exercise(n_path);
  arma::uword next;
  for (int tt = (n_dec - 2); tt >= 0; tt--) {
    primals = discount * primals;
    duals = discount * duals;
    states = path.slice(tt);
    compare.col(0) = strike - states.col(0);
    exercise = arma::max(compare.cols(0, 1), 1);
    // Primal values
    for (std::size_t ii = 0; ii < n_path; ii++) {
      if (path_policy(ii, tt) == 1) {
        primals(ii) = exercise(ii);
      } else {
        primals(ii) = mart(ii, tt) + primals(ii);
      }
    }
    // Dual values
    compare.col(2) = mart.col(tt) + duals;
    duals = arma::max(compare, 1);
  }
  return Rcpp::List::create(Rcpp::Named("primal") = primals,
                            Rcpp::Named("dual") = duals);
}
