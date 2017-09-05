// Copyright 2017 <jeremyyee@outlook.com.au>
// Pricing Bermuda put option using LSM
////////////////////////////////////////////////////////////////////////////////

#include "inst/include/basis.h"

// Regression coefficients from singular value decomposition
arma::vec SVDCoeff(const arma::mat& xreg,
                   const arma::vec& yreg) {
  arma::mat u;
  arma::vec s;
  arma::mat v;
  arma::svd_econ(u, s, v, xreg, "both", "dc");
  arma::vec d = (u.t()) * yreg;
  arma::vec temp(xreg.n_cols, arma::fill::zeros);
  int r = rank(xreg);
  temp(arma::span(0, r - 1)) =
      d(arma::span(0, r - 1)) / s(arma::span(0, r - 1));
  arma::vec svd_coeff = v * temp;
  return svd_coeff;
}

// Pricing Bermuda put using LSM
//[[Rcpp::export]]
arma::vec BermudaPutLSM(const arma::cube& path,
                        const double& strike,
                        const double& discount,
                        const arma::umat& basis,
                        const bool& intercept,
                        const std::string& basis_type) {
  // Extract parameters
  const std::size_t n_dec = path.n_slices;
  const std::size_t n_path = path.n_rows;
  const std::size_t n_dim = path.n_cols;
  // Extract information about regression basis
  std::size_t n_terms = arma::accu(basis);  // Number of features in basis
  if (intercept) { n_terms++; }
  arma::uvec reccur(basis.n_rows);
  reccur = ReccurLimit(basis);
  // Perform the Bellman recursion starting at last time epoch
  arma::mat states = path.slice(n_dec - 1);
  arma::mat path_values(n_path, 2, arma::fill::zeros);  // (exercise, value)
  path_values.col(0) = strike - states.col(0);
  path_values.col(1) = arma::max(path_values, 1);
  arma::vec expected(n_terms);  // Regression fit
  arma::mat reg_basis(n_path, n_terms);
  std::size_t n_money;  // number of in_money paths
  arma::uvec sorted(n_path);  // sort index
  arma::uvec cols_value = arma::linspace<arma::uvec>(0, 1, 2);  // col indices
  arma::uvec cols_path = arma::linspace<arma::uvec>(0, n_dim - 1, n_dim);
  arma::mat path_money(n_path, 2);  // Store in_money paths
  arma::mat value_money(n_path, 2);  // Store in_money values, (states, value)
  arma::vec fitted(n_path);  // Fitted continuation value
  // Perform Backward induction
  for (int tt = (n_dec - 2); tt >= 0; tt--) {
    path_values.col(1) = discount * path_values.col(1);  // Discount one step
    // Exercise value
    states = path.slice(tt);
    path_values.col(0) = strike - states.col(0);
    // Select paths that are in the money
    n_money = 0;
    for (std::size_t pp = 0; pp < n_path; pp++) {
      if (path_values(pp, 0) > 0) {
        n_money++;
      }
    }
    if (n_money == 0) {
      continue;
    }
    n_money--;  // Change to C++ indexing
    // Sort by sorted
    sorted = arma::sort_index(path_values.col(0), "descend");
    path_money.submat(0, 0, n_money, n_dim - 1) =
        states.submat(sorted.subvec(0, n_money), cols_path);
    value_money.submat(0, 0, n_money, 1) =
        path_values.submat(sorted.subvec(0, n_money), cols_value);
    // Regression basis
    if (basis_type == "power") {
      reg_basis.submat(0, 0, n_money, n_terms - 1) =
          PBasis(path_money.submat(0, 0, n_money, n_dim - 1), basis, intercept, n_terms, reccur);
    } else if (basis_type == "laguerre") {
      reg_basis.submat(0, 0, n_money, n_terms - 1) =
          LBasis(path_money.submat(0, 0, n_money, n_dim - 1), basis, intercept, n_terms, reccur);
    }
    expected = SVDCoeff(reg_basis.submat(0, 0, n_money, n_terms - 1),
                        value_money.col(1).subvec(0, n_money));
    // The fitted continuation value
    fitted.subvec(0, n_money) = reg_basis.rows(0, n_money) * expected;
    for (std::size_t ww = 0; ww < n_money; ww++) {
      if ((strike - path_money(ww, 0)) > fitted(ww)) {
        path_values(sorted(ww), 1) = strike - path_money(ww, 0);
      }
    }
  }
  return path_values.col(1);
}
