#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Find r-squareds per combination for linear regression
//' 
//' @description
//' `KDRsolveCombos_cpp` finds the linear regression r-squared for each of many 
//'   combinations of variables.  It can handle multiple dependent variables
//'   at once.
//' 
//' `KDRsolveCombos1_cpp` does the same but is somewhat faster when there is only one 
//'   dependent variable.
//' 
//' @details [KDRsolveCombos()] is a friendlier wrapper for these functions; it 
//'   adds adjusted r-squareds if desired and automatically chooses which 
//'   of these to call.
//' 
//' @param combo An integer matrix of combinations to solve for, one column
//'   per combination and as many rows as the size of the combinations (all of
//'   which must be the same size on any one call).
//' @param scpMatrix A sum-of-cross-products matrix, from which the constant
//'   term in the regression has aleady been swept out.  (Easily accomplished by
//'   pre-centering all variables before creating the SCP matrix.)
//' @param depvar For KDRsolveCombos_cpp, an integer vector of the row/column
//'   positions in `scpMatrix` of the dependent variables 
//'   (numbered `1..N` as in R).
//'   For KDRsolveCombos1_cpp, a single integer indicating 
//'   the only dependent variable.
//' @returns For `KDRsolveCombos_cpp`, a matrix of r-squareds found, with 
//'   one column per dependent variable and one row per combination in combo.  
//' 
//' For `KDRsolveCombos1_cpp`, a vector of r-squareds found, one value for each 
//'   combination in combo.
//'
//' @export
// [[Rcpp::export]]
arma::mat KDRsolveCombos_cpp(const arma::imat& combo,
                       const arma::mat& scpMatrix,
                       const arma::uvec& depvar) {
  
  // depvar is assumed to be 1-based coming from R
  arma::uvec dv = depvar - 1;
  
  const arma::uword ncomb = combo.n_cols;     // number of combos to figure for
  const arma::uword k = combo.n_rows;         // size of all combos
  const arma::uword ld = dv.n_elem;           // number of dependent variables
       // KDRsolveCombos1_cpp is up to 5-10% faster if just one dependent
  
  arma::mat rsq(ncomb, ld, arma::fill::zeros);  // Our final result:
                                              // one r-squared per combo per DV
  
  arma::vec vstart(ld);                       // starting r-squared, per DV
  for (arma::uword j = 0; j < ld; ++j)
    vstart(j) = scpMatrix(dv(j), dv(j));
  
  for (arma::uword icol = 0; icol < ncomb; ++icol) {   // the big loop
    
    // combo indices (convert 1-based → 0-based)
    arma::uvec ccomb = arma::conv_to<arma::uvec>::from(combo.col(icol) - 1);
    
    // Extract needed submatrices for this combo
    arma::mat SCP_xx = scpMatrix.submat(ccomb, ccomb);
    arma::mat SCP_xy = scpMatrix.submat(ccomb, dv);
    arma::mat SCP_yx = scpMatrix.submat(dv, ccomb);
    
    // Solve SCP_xx * B = SCP_xy
    arma::mat B = arma::solve(SCP_xx, SCP_xy, arma::solve_opts::fast);
    
    // Compute diagonals of (SCP_yx * B) and r-squareds from them.
    for (arma::uword j = 0; j < ld; ++j) {
      double acc = 0.0;
      for (arma::uword i = 0; i < k; ++i)
        acc += SCP_yx(j, i) * B(i, j);
      rsq(icol, j) = acc / vstart(j);
    }
  }
  
  return rsq;
}

//' @rdname KDRsolveCombos_cpp
//' @export
// [[Rcpp::export]]
arma::vec KDRsolveCombos1_cpp(const arma::imat& combo,
                              const arma::mat& scpMatrix,
                              unsigned int depvar) {
  
  // This is the only-one-dependent variable version of KDRsolveCombos_cpp.
  
  // depvar is 1-based from R
  const arma::uword dv = depvar - 1;
  
  const arma::uword ncomb = combo.n_cols;
  const arma::uword k     = combo.n_rows;
  
  arma::vec rsq(ncomb, arma::fill::zeros);   // final r-squared, one per combo
  
  // s_yy
  const double vstart = scpMatrix(dv, dv);
  
  arma::uvec ccomb(k);
  arma::mat SCP_xx(k, k);
  arma::vec s_xy(k);
  arma::rowvec s_yx(k);
  arma::vec b(k);
  
  for (arma::uword icol = 0; icol < ncomb; ++icol) {
    
    // load combo (1-based → 0-based)
    arma::uvec ccomb = arma::conv_to<arma::uvec>::from(combo.col(icol) - 1);
    
    // extract submatrices cleanly
    SCP_xx = scpMatrix.submat(ccomb, ccomb);
    
    s_xy = scpMatrix.submat(ccomb, arma::uvec{dv});
    s_yx = scpMatrix.submat(arma::uvec{dv}, ccomb);
    
    // solve SCP_xx * b = s_xy
    arma::solve(b, SCP_xx, s_xy, arma::solve_opts::fast );
    
    // dot product
    rsq(icol) = arma::dot(s_yx, b) / vstart;
  }
  
  return rsq;
}
