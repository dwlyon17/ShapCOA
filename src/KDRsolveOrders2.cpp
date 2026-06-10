#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Compute r-squareds for many different orders of variables
//' 
//' This function computes r-squared for each step, in order, in a list
//'   of variables, for many different orderings of variables.  It can 
//'   handle multiple dependent variables.
//'
//' This function combines an orderings-based driver with a regression-based
//'   valuation function in C++ and 
//'   is the fastest option in the `ShapCOA` package to get 
//'   results for orderings-based linear key driver regressions.
//'   
//' It does not directly produce adjusted r-squareds or Shapley Values, 
//'   but the r-squareds it returns are used by its usual caller to get those.
//'   
//' It can handle size-limited Shapley Values in that the length of each
//'   ordering need not include all variables (so, length is the
//'   `maxitems` of some other functions) and the caller can enforce
//'   the `minitems` idea.  
//'   
//' @param orders An integer matrix of orderings, one per column, 
//'   often from [getCOA()].
//' @param scpMatrix Cross-products matrix, without a constant 
//'   (either built with all pre-centered variables, or with the constant
//'   already swept out and removed from the matrix).  
//'   The dependent variable(s) must be the last row/column(s) in `scpMatrix`.
//'    See [cpBuild()].
//' @param ndep Integer number of dependent variables.  If zero, it is 
//'   calculated as the number of columns in `scpMatrix` minus the 
//'   number of rows in `orders` (i.e., assuming all items in all orders).
//' @returns An array of R-squared values, with first two dimensions 
//'   matching those of the `orders` parameter and a third dimension 
//'   equal to `ndep`.
//' @export
// [[Rcpp::export]]
arma::cube KDRsolveOrders_cpp ( const arma::umat& orders,
                               const arma::mat& scpMatrix, 
                               int ndep=0 )  {
  int p = scpMatrix.n_rows;    // total number of variables
  int maxlen = orders.n_rows;  // number of X's from length of each ordering
  if ( ndep == 0 )  ndep = p - maxlen;  
                               // scpMatrix is ONLY X's & Y's, in that order
  int norders = orders.n_cols; // Number of orderings to do
  
  
  // =================== Setup what we can before the combo loop
  // Dependent variable indices (last ndep)
  arma::uvec depvar(ndep);            // row/col of dependents in original SCP
  arma::uvec depnew(ndep);            // row/col of deps in re-arranged CP
  for (int k = 0; k < ndep; k++) {
    depvar(k) = p - ndep + k;
    depnew(k) = maxlen + k;           // different if shortened orderings
  }
  
  // Build sweep order (0-based).  
  int pwork = maxlen + ndep;          // Size of working CP matrix
           // (smaller if orderings are not full-length)
  arma::uvec ordersu(pwork);          // First maxlenl will vary by ordering
  ordersu.tail(ndep) = depvar;        // Dependent vars always unswept and last
  
  // Save original Corrected SSQs for base in r-squared computations.
  arma::vec cssorig(ndep);               // starting SSQ for r-sq denominator
  for (int k = 0; k < ndep; k++)
    cssorig(k) = scpMatrix(depvar(k), depvar(k));  // from original positions
  
  // Output matrix: One r-sq per ordering step per dependent variable
  arma::cube rsq(maxlen, norders, ndep, fill::zeros);
  
  // =================== Now ready to loop over orderings
  
  for ( int iorder = 0; iorder < norders; iorder++ )  {
  
    // First, reorder cross-product matrix to match sweep order.  
    ordersu.head(maxlen) = orders.col(iorder) - 1;  // leaves depvars at end
    arma::mat CP = scpMatrix(ordersu, ordersu);     // subset/re-order CP
  
    // ===== Sweep loop
    for (int ipiv = 0; ipiv < maxlen; ipiv++) {
    
      double pivot = CP(ipiv, ipiv);
    
      // Partial Gauss–Jordan sweep
      for (int i = ipiv + 1; i < pwork; i++) {
        double factor = CP(i, ipiv) / pivot;
        CP(i, span(ipiv + 1, pwork - 1)) -=
          factor * CP(ipiv, span(ipiv + 1, pwork - 1));
        CP(i, ipiv) = -factor;
      }
    
      CP(ipiv, span(ipiv + 1, pwork - 1)) /= pivot;
    
      // ===== R-squared bookkeeping, still inside sweep loop
      for (int k = 0; k < ndep; k++) {
        // R-squared 
        rsq(ipiv,iorder,k) =  1.0 - CP(depnew(k), depnew(k)) / cssorig(k);
      }      // end of loop over dependents
    }        // end of sweep loop
  }          // of loop over orders
  
  return rsq ;  
}
