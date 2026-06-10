# SVs for TURF on MaxDiff functions

#' Approximate Shapley Values for TURF Analysis on MaxDiff data
#' 
#' This function computes Shapley Values using COAs or other orderings of
#'   items, for various forms of TURF on MaxDiff data
#'   (and for standard 0-1 TURFs).
#'   It uses the valuation function [SVsMaxDiffCalc_cpp()].  
#'   
#' @details For standard 0-1 TURFs and many variants of it, including depth>1 and many
#'   measures other than standard reach, this function is slow and unnecessary.
#   Instead, use fast exact SV calculations from [SVsForTURF()].
#'   
#' Three major versions of scoring combinations of MaxDiff items are supported,
#'   specified by the `xform` parameter:
#'   multinomial logit (MNL), the default Sawtooth individual probability score
#'   and the (default for anchored MaxDiffs) 
#'   Sawtooth anchored probability score.
#'
#' Two ways to measure the TURF "reach" are directly supported: the total
#'   weighted probability scores (WPS), and whether the WPS clear a specified
#'   threshold ("post-thresholding".  
#'   
#' A third measurement option is standard 0-1 TURF resulting 
#'   from "pre-thresholding," also known as "item thresholding", 
#'   where a threshold is applied in advance to each individual item, 
#'   based on either utilities or some kind of score.  
#   However, this results in a standard 0-1 TURF, 
#   for which Shapley Values area better handled by [SVsForTURF()].
#'   If this is done here, the `utilities` parameter must contain only `0's` 
#'   and `1's`. In other words, the pre-thresholding must be done by the caller;
#'   it is not handled by this function.
#'        
#' Weighted data is supported, as are size-limited Shapley Values.  
#'
#' The `orders` need not be full-length (i.e., need not have every item
#'   occuring in every ordering), so this works for "size-limited Shapley 
#'   Values".  The number of items is deduced from the maximum value appearing
#'   in `orders`.  If `maxitems` is specified but `orders` are full-length,
#'   only the first `maxitems` rows are used.
#' 
#' `tasksize` only applies if `xform` is 1, since the default weighted 
#'    probability calculation uses the task size.  
#'    It is required in that case.
#'
#' ## MaxDiff Scores
#' `xform` controls how utilities are centered or anchored and how
#'   they are transformed into scores for the TURF analysis:  
#'  
#'   * {`xform=0   `}  MNL-based scoring.  Utilities are not centered. 
#'            After exponentiating, individual scores are re-percentaged to 
#'            add to `100`.  The score for a combination of items is the sum
#'            of the individual item scores. 
#'   * {`xform=1   `}  The default individual utility probability 
#'            scoring.  Utilities are first zero-centered.  The scores are
#'            `E / ( E + a - 1 )`, where `E` is the exponentiated utility for 
#'              an item and `a` is the number of items 
#'              in the original tasks, given by `tasksize`.
#'            For the score for a combination of items, `S` is the sum of the 
#'              exponentiated utilities for the items in the combination and
#'              then calculated as `S / ( S + a - 1 )`. 
#'   * {`xform=2   `}  Anchored probability scores.  The utilities
#'            are assumed to have already been zero-anchored and are not 
#'            further centered.  Scores for items and combinations are
#'            the same as for `xform=1` except that the final transformation
#'            is `S / ( S + 1)` and the task size is not used. 
#' 
#' Note that the scores used for the TURF are *not* the same as the 
#'   "rescaled weighted probability scores" Sawtooth Software 
#'   uses for *aggregate* reporting.  
#'   The default rescaled WPS uses the `xform=1` calculation for the items
#'   detailed above, then rescales the results to add to 100 before averaging
#'   over respondents.  The anchored WPS uses the same calculation 
#'   at the respondent level as the `xform=1` case 
#'   (so the tasksize does matter now) and then multiplies the result by 
#'   `100 / (1/a)` so that items with utility equal to the anchor get a final
#'   score of `100` and the maxium score is `100a`.  
#'   
#' @param orders An integer matrix of orderings, one per column, often a COA 
#'   from [getCOA()].  
#'   If `NULL`, a COA will be generated via [getCOA()], using the 
#'   `ncoa` parameter.`  
#' @param utilities MaxDiff logit *utilities* (**not** transformed *scores*).  
#'   Need not be zero-centered or zero-anchored.  Alternatively, for 
#'   standard 0-1 TURF, the `0/1` scores.
#' @param xform How to transform MaxDiff utilities into scores for analysis.  
#'   See details below.  In brief, 0 means MNL, 
#'   1 means default individual utility probability score,
#'   2 means anchored probability score, and
#'   NULL means use sums of unexponentiated utilities (i.e., 0/1 data).
#' @param anchor If numeric, the column number of the anchor item 
#'   in an anchored MaxDiff;
#'   it will be used for centering and then removed.
#'   If `FALSE`, `utilities` are not from anchored data.  
#'   If `TRUE`, `utilities` are anchored but have *already* been zero-anchored
#'     and the anchor column has already been removed.
#' @param tasksize The number of items per task in the original MaxDiff tasks.
#'   This is the "a" in the published Sawtooth formulas; 
#'     it is used in computing the weighted probability scores when `xform` 
#'     is 1 (unanchored MaxDiff).
#'   Required if `xform` is `1`, ignored otherwise.
#' @param threshold If TRUE or numeric, and `xform` is not `NULL`, 
#'   weighted probability scores are subjected to a *post*-threshold, 
#'   either of the specified non-zero numeric value, 
#'     or of a standard default of `0.5` for MNL (`xform=0`),
#'     `0.9` for unanchored WPS (`xform=1`) 
#'     or `0.5` for anchored Maxdiff (`xform=2`).
#'   If `FALSE`, TURF calculations are on the weighted probability scores 
#'     themselves, without thresholding.  
#'   If `xform` is NULL, this is the "depth" parameter for a standard 
#'     `0/1` TURF.
#' @param weights Vector of respondent weights, one per row of `df`.  
#'   Unweighted data if `NULL`.  Weights may be zero to subset data but should
#'   never be negative.
# Next is for minitems and maxitems
#' @inheritParams SVsByOrders
#' @param ncoas How many COAs to "stack" for extra precision, 
#'   applicable only if `order` is `NULL` and the COA is generated by 
#'   this function.
#' @param ... Available to pass `random` or `silent` parameters 
#'   onward to [getCOA()].  Irrelevant `orders` is supplied.
#' @returns A vector of approximate Shapley Values, one value per column
#'   in `utilities`.
#' @export
SVsForMaxDiff <- function ( orders=NULL, utilities, xform=NULL, 
                            anchor=FALSE, tasksize=0, 
                            threshold=TRUE, weights=NULL, 
                            minitems=0, maxitems=0, ncoas=1, ... )  {
  # Note: orders is typically left NULL for auto COA generation, but must be
  #       first to work with future parallel options.
  
  # Fix up defaulted xform
  if ( is.null(xform) )  {       # assume 0/1 unless evidence of anchoring
    if ( is.logical(anchor) && anchor )  { xform <- 2
    } else if ( is.numeric(anchor) && anchor > 0 )  { xform <- 2
    } else xform <- -1    # easier internally, for 0/1 TURF
  }
    
  # ==== Deal with utilities.  Check, center, exponentiate, rescale ...
  if ( min(utilities) >= 0.0 && xform != -1  )  {
      stop ( "SVsForMaxdiff: utilities must be raw utilities, ",
             "not exponentiated or scored" )
  }
  
  # Anchoring first
  if ( is.numeric(anchor) )  {
    if ( anchor > ncol(utilities) )  {
      stop ( "SvsForMaxDiff: anchor must be a column number in utilities.\n" )
    }
    anchname <- colnames(utilities)[anchor]
    if ( is.null(anchname) )  anchname <- paste0 ( "Column ", anchor )
    cat ( "SVsForMaxDiff: Zero-anchoring on", anchname, ".\n" ) 
    utilities <- ( utilities - utilities[,anchor] ) [,-anchor]
                                  # zero-anchor and toss the anchor
    if ( is.null(xform) )  xform <- 2    # default xform
    anchor <- TRUE                # convert numeric to logical now
  }
  if ( anchor && xform != 2 )  {
    stop ( "SVsForMaxDiff: anchor is TRUE but xform is not 2." )
  }
  
  # Standard unanchored
  if ( xform == 1 )  {
    if ( is.null(tasksize) || tasksize == 0 )  {
      stop ( "SVsForMaxDiff: tasksize= must be specified ",
             "when xform=1 for default MaxDiffs" );
    }
    utilities <- utilities - rowMeans(utilities)    # zero-center
  } 
    
  if ( xform == -1 )  {              # not MaxDiff data!
    if ( !isTRUE(all.equal ( range(utilities), c(0,1) )) )  {
      cat ( "SVsForMaxDiffs WARNING:", "utilities for standard 0/1 TURF",
            "(xform=NULL) are NOT only 0's and 1's.\n" )
    }
  } else utilities = exp(utilities)  # exponentiate now-centered MaxDiff utils
  
  # Simple MNL-- rescale to (0,1)
  if ( xform == 0 )  utilities = utilities / rowSums(utilities)

  # ==== Generate the COA needed, or check any orderings caller gave us.
  #      Apply the maximum items, if it matters.
  nitems <- ncol(utilities)         # # items (less any original anchor)
  if ( is.null(maxitems) )  maxitems <- 0
  if ( maxitems == 0 )  maxitems <- nitems
  if ( is.null(orders) )  {
    orders <- getCOA ( nitems, ncoas, ... )  
                                        # Get (possibly stacked) COA we need
  } else {                              # If actually given orderings, not usual
    if ( nitems < nrow(orders) )  {     # Nonsense!
      stop ( "SVsForMaxDiff: Orderings given are longer than ",
             "number of columns of utilities" )
    }                                   # implied maxitems
    if ( nitems != nrow(orders) )  {    # sizes don't match
      if ( nrow(orders) != maxitems )  {
        stop ( "SVsForMaxDiff: Length of each ordering is neither ",
               "number of items nor maxitems." )
      } else {                          # legit implied maxitems
        cat ( "SVsForMaxDiff: Orderings are not full for all variables.",
              " Size-limited SVs will be done.\n" )
        maxitems <- nrow(orders)
      }
    }
  } 
  if ( maxitems != nitems )  orders <- orders[1:maxitems,,drop=FALSE]
  
  # ==== Manage the weights
  if ( is.null(weights) )   weights <- rep(1,nrow(utilities)) 
  if ( !checkWeights ( weights, nrow(utilities), 
                       "SVsOnMaxDiff" ) )  return (NULL)
  
  # ==== Default the post-threshold
  if ( is.null(threshold) )  threshold <- TRUE    # default if omitted
  if ( threshold == 0 )  threshhold <- TRUE       # default, should do it
  if ( is.logical(threshold) )  {
    if ( threshold )  {
      if ( xform == -1 )  {        # a standard 0/1 TURF
        threshold <- 1
        cat ( "SVsForMaxDiffs: defaulting to threshold or 'depth' of 1",
              "for standard 0/1 TURF." )
      } else  threshold <- c ( 0.5, 0.9, 0.5 ) [xform+1]
                                 # standard numeric defaults
    } else threshold <- 0        # no threshold wanted
  }                              # else numeric, leave the value alone
  
  # Compute the total values, result same shape as orders.
  values <- SVsMaxDiffCalc_cpp ( orders, utilities, weights, 
                                 maxdiff=(xform!=-1), 
                                 xform=xform, threshold=threshold, 
                                 tasksize=tasksize )
  
  # Turn totals into increments and apply minimum items, if any
  values <- increments ( values )          # We get back 2-dimensional only
  if ( minitems > 0 )  {                   # drop items below the minimum
    drops <- - (1:minitems)
    values <- values[drops,,drop=FALSE]   
    orders <- orders[drops,,drop=FALSE]    # in case min=max, don't lose a dim
  }
  
  # Add up totals, get means for results.  
  # First get counts for denominator.  
  dim(values) <- c ( dim(values), 1 )     # Make it 3-D for totalup calls
  cdummy <- array ( 1, dim(values) )
  counts <- totalup ( cdummy, orders )[,1]
  svs <- totalup ( values, orders ) / counts
  svs[,1]
}


