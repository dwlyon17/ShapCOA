#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Valuation function for various forms of TURF and TURF on MaxDiff analyses
//' 
//' This function computes the the total value of combinations of various
//'   forms of TURF analysis for multiple orderings of items.
//'   It is most useful for TURF on MaxDiff data.  
//'   It is tightly coupled with [SVsForMaxDiff()], which handles the final
//'     aggregation into Shapley Values, including differencing and 
//'     unrotating the orderings.  
//'     It relies on its caller to check parameters and center or exponentiate
//'     utilities as needed.
//'     Even though it is exported and available,
//'     it is not *intended* for use in other situations.
//'     
//' @param orders An integer matrix of orderings, one per column,
//'   typically from [getCOA()].
//' @param utilities The exponentiated MaxDiff utilities, in a numeric 
//'   matrix with one column per item and one row per respondent.  
//'   For standard `0/1` TURF, the `0` and `1`.
//' @param weights Vector of respondent weights, one per row of `utilities`.  
//'   Weights may be zero to subset data but should never be negative.
//'   The parameter is required; use a vector of 1's for unweighted data.
//' @param maxdiff Logical, whether `utilities` contains exponentiated MaxDiff
//'   utilities (the default).
//'   If not, `0-1` reach data is assumed but is not checked.
//' @param xform How to transform total exponentiated utilities 
//'   scores for a combination into the final score for the combination.  
//'   Ignored if `maxdiff` is `FALSE`.  
//'   See [SVsForMaxDiff()] for details.
//' @param threshold For MaxDiff, the numeric threshold to use for 
//'   post-thresholding of MaxDiff scores.  
//'   If `maxdiff` is `FALSE`, the "depth" to use in evaluating a standard
//'   `0/1` TURF on "reach".
//' @inheritParams SVsForMaxDiff
//' @returns  A `double` matrix of combination values, with the same shape and
//'   size as the `orders` matrix.
//' @export
// [[Rcpp::export]]
arma::mat SVsMaxDiffCalc_cpp ( const arma::umat& orders, arma::mat& utilities, 
                               const arma::vec& weights, 
                               const bool maxdiff=true,  const int xform=1,
                               double threshold=0.0,  const int tasksize=0 )  {
  int nitems = utilities.n_cols;
  int nresp = utilities.n_rows;   
  int norders = orders.n_cols;      
  int ordlen = orders.n_rows;       // less than nitems if size-limited
  
  double wtsum = sum(weights);         // Used to get weighted means
  
                            // mat for results MUST be outside the order loop
  arma::mat scores(ordlen, norders);  // scores by position by ordering
  scores.fill(arma::datum::nan);      
  // following are used within the loop, but pre-allocated for efficiency
  arma::uvec orderk(ordlen); 
  arma::vec sums_sofar(nresp);                // running total utilities
  arma::mat rawsums(nresp,ordlen);            // total utils per R per step
  
  for ( int iorder = 0; iorder < norders; iorder++ )  {
    orderk = orders.col(iorder) - 1;          // temp, curr order, 0-based
    sums_sofar.zeros();
    for ( int iord = 0; iord < ordlen; iord++ )  {  // THIS is the "quick loop"
      sums_sofar += utilities.col(orderk(iord));    //   over resps is implicit
      rawsums.col(iord) = sums_sofar;         // total score per R thru orders
    }
    // Now back in the orders loop only.  Convert utility sums to final
    //   "values" of the different combinations along the way.
    if ( maxdiff )  {
      switch (xform) {                 // for MaxDiff, transform sum(exps)
        case 0: break;                        // MNL-style, nothing to change
        case 1: rawsums /= ( rawsums + tasksize - 1 );  
                break;                 // Sawtooth non-anchor default
        case 2: rawsums /= ( rawsums + 1 ); 
                break;                 // Sawtooth anchored
      }
    }

    if ( threshold > 0 )                // apply any post-thresholding 
         scores.col(iorder) = ( conv_to<mat>::from(rawsums >= threshold).t() 
                                 * weights ) / wtsum;
    else scores.col(iorder) = ( rawsums.t() * weights ) / wtsum;
                                        // weight and get wtd means
  }                        // end of the orderings loop
  
  return scores;           // let caller put it altogether!
                           // these are not differenced, not unrotated
}
