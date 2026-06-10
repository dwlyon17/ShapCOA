#' Evaluate a standard 0-1 TURF on a single combination
#'
#' This is the simplest-minded possible TURF valuation function, intended
#'   mostly for testing of various drivers.
#'
#' Inputs are not checked in any way; the caller is just trusted.
#' 
#' Drivers should be told `multi=FALSE` and `nonull=TRUE` when applicable.
#' 
#' @param comb The combination to evaluate, a vector of item numbers.
#' @param tdata The 0/1 data, a matrix with one column per item 
#'   and one row per respondent.
#' @param depth The depth required for an overall reach 
#'   (i.e., the number of indiviudal item reaches required for an 
#'   overall reach).
#' @param wts The respondent weights, a vector with one element per row of
#'   `tdata`, or a constant 1 for unweighted data.
# @returns The percentage of respondents reached at the given depth, a scalar.
#' @export
turf1 <- function ( comb, tdata, depth=1, wts=1 )  {
  # scoring/value function for TURF with specifiable depth, with weighting
  
  hits <- rowSums ( tdata[,comb,drop=FALSE] )
  if ( length(wts) > 1 )  { turf <- sum ( wts * ( hits >= depth ) ) / sum(wts)
  } else turf <- mean ( hits >= depth ) 
}
