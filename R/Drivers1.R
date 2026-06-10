# Drivers for arbitrary valuation functions. 

# Includes production-usable orders-based and combinations-based drivers.
# Includes an orders-based driver for the simplest-to-write possible 
#   valuation function.
# Includes an orders-based driver specifically for linear KDRs, 
#   implementing adjusted r-squareds and multiple dependent variables.
# All allow minimum and maximum numbers of items for 
#   size-limited Shapley Values.
# Also some internal utility functions for the orderings-based drivers.

#' Exact, combinations-based, Shapley Values on an arbitrary valuation function
#' 
#' This function calculates exact Shapley Values using the combinations 
#'   approach, for any user-supplied valuation function.  
#' 
#' This is a general driver for the combinations approach.  That approch is 
#'   typically not computationally feasible for more than 25 to 30 items, but 
#'   it does yield exact, not approximate, Shapley Values.
#' 
#' For linear Key Driver Regressions,
#'   [KDRsolveCombos()] is a natural counterpart to this function. 
#'    
#' ## The Valuation Function:
#' The first argument to the value function is the combination(s) it must
#'   compute for.  Further arguments are provided through the `...` mechanism.
#'    
#' If `multi` is 'FALSE`, the value function will be called once for each 
#'   combination, with its first parameter being the combination *vector*.  
#'   It would usually return a single numeric result.
#' 
#' The value function can return multiple measures, if desired, as a vector.
#'   The number of measures need not be specified in parameters; 
#'   it is determined by a preliminary dummy call to the function.
#'   
#' If `multi` is TRUE, the value function will be called with an entire
#'   *matrix* of combinations, one per column.  On any one call, all 
#'   combinations will be of the same size.
#'   The value returned is usually a vector with one value per combination.
#'   one row per position in the orderings and one column per ordering 
#'   In the case of multiple measures, the return value is an array with
#'   one row per combination and one column per measure.
#'   
#' Some valuation functions cannot handle the `NULL` combination and many do
#'   not need to.  But logistic regression, for example, might need to evaluate 
#'   a NULL combo to establish a non-zero log-likelihood starting point for 
#'   other variables to move from.  The `nonull` parameter specifies whether
#'   to call the value function with a null combination, as opposed to just
#'   assuming its value is zero.
#' 
#' @param vfunc The valuation function (see below).
#' @param nitems The number of items or variables in the problem.
#' @param maxitems The maximum combination size to consider, defaulting to
#'   all items.  If this is 5, for example, value increments from size 4
#'   to 5 are included, but from size 5 to 6 are not.  
#'   If 0 or `NULL`, the value of `nitems` is used.
#' @param minitems The minimum combination size to consider, defaulting to 
#'   zero.  If this is 2, for example, value increments from combinations of
#'   size 2 to 3 are included in the calculation, but increments from size 0 
#'   to 1 and size 1 to 2 are not.
#' @param multi Whether `vfunc` can be called with a matrix of many combinations
#'   of a given size, or if `FALSE`, must be called for each combination
#'   individually.  Solutions are a bit faster with value functions written to
#'   allow this.
#'   If `vfunc` is [KDRsolveCombos()], `multi` must be `TRUE`.
#' @param nonull Whether `vfunc` can/should be called for the `NULL` combinaton,
#'   or a value of zero should simply be assumed in that case.  
#'   If `vfunc` is [KDRsolveCombos()], `nonull` must be `TRUE`.
#' @param silent Whether to print progress/information messages to console 
#'   as the function works.
#' @param ... Additional parameters to pass to the valuation function `vfunc`.
#' @returns A matrix of exact Shapley Values on the `vfunc` value increments, 
#'   one column per result type provided
#'   by `vfunc`, and one row per item.  If only one result is provided by 
#'   `vfunc`, a vector with one value per item.
#' @importFrom utils combn
#' @export
SVsByCombos <- function ( vfunc=KDRsolveCombos, nitems, 
                          maxitems=nitems, minitems=0, 
                          multi=TRUE, nonull=TRUE, silent=FALSE, ... )  {
  # Timing tests show no significant speed loss from adding max, min, nonull
  #   and multi features to this.
  # We work with only two combination sizes at once to prevent memory issues
  #   from trying to do all at once.  
  # Basic approach borrowed from relaimpo, but more memory-efficient.
  
  if ( is.null(maxitems) || maxitems == 0 )  maxitems <- nitems
  if ( maxitems > nitems || minitems >= maxitems )  {
    cat ( "SVsByCombos: For", nitems, "items, maximum of", maxitems,
          "and minimum of", minitems, "are invalid.\n" )
    return ( NULL ) 
  }
  
  # Deduce number of scores returned from vfunc (typically meaning number 
  #   of different dependent variables, but could be alternate scores 
  #   like r-squared and adjusted r-squared as well).
  # Call with one combination and see what happens.
  if ( multi )  { combo <- matrix(1:2,2,1)     # one combo as a matrix
  } else combo <- 1:2                          # or as a vector
  ndv <- NCOL ( vfunc ( combo, ... ) )         # see what comes back!
  if ( !silent )  cat ( "Value function returns", ndv, "values per item.\n" )
  
  contrib <- matrix ( 0, nitems, ndv )
  
  c1 <- combn ( nitems, 0 )        
  if ( minitems == 0 && nonull )  {
    r1 <- matrix ( 0, 1, ndv )
    cstart <- 1                    # don't do null in main loop
  } else cstart <- minitems
  
  for ( csize in cstart:maxitems )  {
    c2 <- combn ( nitems, csize )      # combinations of next size up
    if ( !silent )  cat ( "Working in combinations of", csize, "\n" )
    if ( multi )  { r2 <- vfunc ( c2, ... )
    } else {
      r2 <- matrix ( NA, ncol(c2), ndv )
      #      r2 <- rep ( NA, ncol(c2) )         # pre-allocate results
      for ( combo in 1:ncol(c2) )  r2[combo,] <- vfunc ( c2[,combo], ... )
    }
    if ( csize > minitems )  {       # if have results for 2 sizes already
      cc1 <- 1:ncol(c1)              # will use repeatedly in setdiff ...
      for ( j in 1:nitems )  {
        wj <- which ( colSums(c2==j) > 0 )    # larger ones WITH j
        noj <- setdiff ( cc1, which ( colSums(c1==j) > 0 ) )
        # smaller ones without j
        contrib[j,] <- contrib[j,] + colMeans ( r2[wj,,drop=FALSE] ) - 
          colMeans ( r1[noj,,drop=FALSE] ) 
      }
    }
    # Larger will now be the smaller on next pass of loop
    c1 <- c2
    r1 <- r2
  }
  
  if ( ndv == 1 )  contrib <- contrib[,1]      # Usually return vector 
  contrib <- contrib / ( maxitems - minitems ) # convert totals to averages
  contrib                                      # That's the result
}

#' Turn ordered result value measures into incremental values
#'
#' This internal function "diffs" a set of total combo values to obtain
#'   the increment added by each item.  It presumes a starting value of
#'   zero, and tolerates even single-value results 
#'   (i.e., for a single combination size).
#'   
#' @param result A  vector, matrix or 3-D array (one result vector 
#'   per column) of results to be differenced.  
#'   Note that the result should be in order by steps taken, 
#'   not by original variable numbers. 
#' @returns A vector or matrix of same shape and size as the parameter, 
#'   of increments across the positions.
#' @keywords internal
increments <- function ( result )  {
  if ( NROW(result) == 1 )  return ( result )
  switch ( length(dim(result)), 
         rbind ( result[1], diff(result) ),                  # vector
         rbind ( result[1,], diff(result) ),                 # matrix
         { for ( i in 1:dim(result)[3] )  {                  # 3-D
            result[,,i] <- rbind ( result[1,,i], 
                                   diff(result[,,i]) )
           }
           result 
         } )
}

#' Total SV item contributions by variable, undoing the ordering they were in
#' 
#' This internal function totals contributions by each variable across all
#'   orderings, effectively "unscrambling" the ordering of resultsof results.
#'
#' @param pieces Numeric array of results, one ordering per column, 
#'   one row per step in `orders` and a third dimension for multiple
#'   measures obtained from the orderings.
#' @param orders The integer matrix of orderings the `pieces` follow,
#'   with the same shape and size as the first two dimensions  of `pieces`.
#' @returns A vector or matrix (one fewer dimension than `pieces`, of totals
#'   for each variable or item (as many as the maximum value in orders)
#'   and as many columns as the third dimension of `pieces`. 
#' @keywords internal
totalup <- function ( pieces, orders )  {
  ndep <- if ( is.na(dim(pieces)[3]) )  1  else dim(pieces)[3]
  nitems <- max(orders)
  totals <- matrix ( 0, nitems, ndep )
  for ( i in 1:ncol(orders) )  {
    totals[orders[,i],] <- totals[orders[,i],] + pieces[,i,]
  }
  totals
}

#' Approximate orderings-based Shapley Values on an arbitrary valuation 
#' function using COAs or other orderings
#' 
#' This function calculates approximate Shapley Values using the orderings 
#'   approach, for an appropriate user-supplied valuation function.
#'   
#' @details The `orders` need not be full-length (i.e., need not have every item
#'   occuring in every ordering), so this works for "size-limited Shapley 
#'   Values".  The number of items is deduced from the maximum value appearing
#'   in `orders`.  If `maxitems` is specified but `orders` are full-length,
#'   only the first `maxitems` rows are used.
#'  
#' ## The Valuation Function:
#' The first argument to the value function is the ordering(s) it must
#'   compute for.  Further arguments are provided through the `...` mechanism.
#'    
#' If `multi` is 'FALSE`, the value function will be called once for each 
#'   ordering, with its first parameter being the ordering vector.  
#'   It would usually return a vector of results, for each combination size
#'   as it steps through the ordering.  It must *NOT* return increments
#'   as opposed to total values; the differencing is done by this function.  
#'   It must *NOT* return values in item-number order, but must keep the
#'   ordering order.  In other words, if item 5 is the first in an ordering,
#'   the first value returned by `vfunc` must be the contribution of item
#'   5 as the first one in.
#' 
#' The value function can return multiple measures, if desired, as a matrix
#'   with one column per measure (e.g., raw and adjusted r-squared, or 
#'   results for two different dependent variables).  The number of measures
#'   need not be specified in parameters; it is determined by a preliminary
#'   dummy call to the function.
#'   
#' If `multi` is TRUE, the value function will be called with the entire
#'   orders matrix (as well as one test call with a single ordering as a 
#'   single-column matrix).  The value returned is usually a matrix with
#'   one row per position in the orderings and one column per ordering 
#'   (in other words, the same shape as the `orders` matrix).  In the case
#'   of multiple measures, the return value has a third dimension, which is
#'   the number of measures.
#'   
#' @param orders An integer matrix of orderings, one per column, often a COA 
#'   from [getCOA()].  
#' @param vfunc The valuation function (see below).
#' @param silent If TRUE, does not comment on number of values returned
#'   by the value function `vfunc`.
#' @param minitems Minimum combination size to consider, defaulting to zero.
#'   If this is 2, for example, contributions from the first two items
#'   in each ordering are excluded. 
#' @param maxitems Maximum combination size to consider, treated as unlimited
#'   if zero or `NULL`. 
#' @param multi Whether `vfunc` can be called with a matrix of many orderings
#'   at once, or if `FALSE`, must be called for each ordering
#'   individually.  Solutions are (usually trivially) faster with 
#'   value functions written to allow this.
#' @param ... Additional parameters to pass to the valuation function `vfunc`.
#' @returns A vector or matrix, depending on which `vfunc` returns, with the
#'   approximate Shapley Values for each item/variable in the orderings and one
#'   column, if a matrix, for each measure returned by `vfunc`.
#' @export
SVsByOrders <- function ( orders, vfunc=KDRsolveOrders_cpp, silent=FALSE, 
                          minitems=0, maxitems=0, multi=FALSE, ... )  {
  # ==== Defaults and checking ...
  maxlen <- nrow(orders)      # length of each ordering
  nitems <- max(orders)       # number of items
  if ( is.null(maxitems) )  maxitems <- 0
  if ( maxitems == 0 )  maxitems <- nitems
  if ( maxitems > 0 )  {             # Cut down if not already done
    if ( maxitems != maxlen )  {     # Not already cut down
      if ( maxlen != nitems )   {    # Not full either?
        stop ( "SVsByOrders: Length of each ordering is neither ",
               "number of variables nor maxitems." )
      } else {
        orders <- orders[1:maxitems,,drop=FALSE]   # cut it down from full
        maxlen <- nrow(orders)
      }
    }
  }

  # ==== Deduce number of distinct scores returned from func (typically
  #      meaning number of different dependent variables, but could be alternate
  #      scores like r-squared and adjusted r-squared as well).
  # Call with one ordering and see what happens.
  if ( multi )  {
    test1 <- vfunc ( orders[,1,drop=FALSE], ... )
    nres <- NCOL(test1)               # # of distinct results
    if ( length(dim(test1)) == 3 )  { nres <- dim(test1)[3]
    } else nres <- 1
    nvres <- dim(test1)[1]            # # results per ordering per measure
  } else {
    test1 <- vfunc ( orders[,1], ... )
    nres <- NCOL(test1)               # # of distinct results
    nvres <- length(test1) / nres     # # of results per ordering per measure
  }
  if ( !silent )  cat ( "SVsByOrders: Value function returns", 
                        nres, "values per item.\n" )
  if ( nvres != dim(orders)[1] )  {
    stop ( "SVsByOrders: Value function returns wrong number of values", 
           "per ordering.  Expected", dim(orders[1]), "but got", nvres, "." )
  }

  # ==== Now the basic work of the value function.
  if ( !multi )  {              
    contribs <- array ( NA, c( ncol(orders), nitems, nres ) )
                                # create array to fill in
    for ( i in 1:ncol(orders) )  contribs[,i,] <- vfunc ( orders[,i], ... )
                                # step through all orders
  } else {                      # one call does it all!
    contribs <- vfunc ( orders, ... )
    if ( !identical ( dim(contribs)[1:2], dim(orders) ) )  {
      stop ( "SVsByOrders: Result from vfunc with multi does not match ",
             "dimensions of the order parameter.\n",
             "orders dim is ", paste ( dim(orders), collapse=" " ), 
             "vfunc result is ", paste ( dim(contribs), collapse=" " ) )
    }
    if ( length(dim(contribs)) == 2 )  dim(contribs) <- c(dim(contribs),1)
                                    # add 3rd dimension for convenience
  }

  # Turn totals into increments and apply minimum items, if any
  contribs <- increments ( contribs )
  if ( minitems > 0 )  {                   # drop items below the minimum
    drops <- - (1:minitems)
    contribs <- contribs[drops,,,drop=FALSE]   
    orders <- orders[drops,,drop=FALSE]    # in case min=max, don't lose a dim
  }

  # Add up totals, get means for results.
  # First get counts for denominators.  Same for all measures.
  cdummy <- array ( 1, c ( dim(contribs)[1:2], 1 ) )
  counts <- totalup ( cdummy, orders )[,1]
  svs <- totalup ( contribs, orders ) / counts
  if ( dim(svs)[2] == 1 )  svs[,1]  else svs
}

#' Simplest possible approximate orderings-based Shapley Values on an 
#' arbitrary combination-based value function using COAs or other orderings
#' 
#' This function evaluates a user-supplied value function over multiple 
#'   orderings of items or variables (likely COAs) and averages the results to
#'   obtain approximate Shapley Values.  
#'   
#' This function will work for any value function that can evaluate 
#'   a single combination, the simplest possible requirement for
#'   a value function.
#'   But, it is **not efficient**  in that such a value function can 
#'   never make any use of any efficiencies  
#'   in successive "steps into" an ordering.
#'   
#' @details The `orders` need not be full-length (i.e., need not have every item
#'   occuring in every ordering), so this works for "size-limited Shapley 
#'   Values".  The number of items is deduced from the maximum value appearing
#'   in `orders`.  Unlike some other drivers, this one does *not* shorten 
#'   orders according to a `maxitems` parameter.
#'
#' Minimum combination sizes are supported by the `minitems` parameter.
#'   
#' ## The Valuation Function:
#' The first argument to the value function is a vector giving the combination 
#'   it must compute for.  
#'   Further arguments to it are provided through the `...` mechanism.
#'   
#' Multiple measures per combination are **not** supported by this driver;
#'   the value function must return a single numeric value.
#'    
#' @param vfunc The valuation function (see below).
#' @param silent Whether to comment to console 
#'   when implied maximum combination size is not same as number of items.
# Next is to inherit orders and minitems
#' @inheritParams SVsByOrders 
#' @param ...  Additional parameters passed to the value function.
#' @returns A numeric vector of the average value increments, which is an 
#'   approximation of the Shapley Values if each ordering covers all
#'   variables and the orderings are some type of random sample (e.g., a COA).
#' @export
SVsByOrders1 <- function ( orders, vfunc, minitems=0, silent=FALSE, ... )  { 
  # Very inefficient and generally unnecessary approach, which is
  # why it is tagged internal.  Would work for user, but better to use
  # an orderings solver.
  
  nitems <- max(orders)      # max instead of nrow handles size-limited SVs
  maxlen <- nrow(orders)     # length of each ordering
  if ( nitems != maxlen && !silent )  {
    cat ( "SVsByOrders1: doing size-limited Shapley Values.\n" )
  }
  contribs <- matrix ( NA, nitems, ncol(orders) )
  
  for ( i in 1:ncol(orders) )  {   # step through all orders
    val <- rep(NA,nitems)
    order1 <- orders[,i,drop=FALSE]
    for ( j in 1:maxlen )  val[j] <- vfunc ( order1[1:j,,drop=FALSE], ... )
    valinc <- diff ( c ( 0, val ) )   # make it incremental
    if ( minitems > 0 )  valinc[1:minitems] <- NA
    if ( maxlen != nitems )  {    # if orderings are not full ...
      order1 <- unique ( c ( order1, 1:nitems ) )
    }                    # be sure all vars get pulled somewhere even if NAs
    contribs[,i] <- valinc [order(order1)]    # put back into variable order
  }    # We save results rather than totaling on the spot to allow for 
  #   shortened orderings, and use the `NA` logic below as well.
  asvs <- rowSums ( contribs, na.rm=TRUE ) / rowSums ( !is.na(contribs) )
  asvs
}

#' Approximate Shapley Values in Key Driver Regressions using COAs 
#'   or other orderings
#' 
#' This function computes SVs for r-squared and/or adjusted r-squared in 
#'   linear key driver regressions.  
#'   It can compute for multiple dependent variables 
#'   and can handle size-limited Shapley Values.
#'
#' This is a wrapper around [KDRsolveOrders_cpp()] that 
#'   adds the options of adjusted r-squared and size-limited SVs to 
#'   the results.
#'   It is the fastest option in the `ShapCOA` package to get 
#'   COA-based approximate Shapley Values on KDRs.
#'   
#' The `orders` need not be full-length (i.e., need not have every item
#'   occuring in every ordering), so this works for "size-limited Shapley 
#'   Values".  The number of items is deduced from the maximum value appearing
#'   in `orders`.  If `maxitems` is specified but `orders` are full-length,
#'   only the first `maxitems` rows are used.
#'   
#' In the (unlikely) event that `orders` gives all possible orderings, 
#'   the SVs computed are exact.
#'   
#' @param scpMatrix Cross-products matrix, *without* a constant 
#'   (either built with all pre-centered variables, or with the constant
#'   already swept out and removed from the matrix).  
#'   The dependent variable(s) must be the last row/column(s) in `scpMatrix`.
#'   See [cpBuild()].
#'   If `adjusted == TRUE`, it must contain 
#'   an `ESS` attribute specifying the total number of respondents or 
#'   effective sample size used to create it.
#' @param ndep Integer number of dependent variables.  If zero, it is 
#'   calculated as the number of columns in `scpMatrix` minus the 
#'   number of rows in `orders` (i.e., assuming all items in all orders).
# Next is to inherit adjusted and both
#' @inheritParams adjustOrNot   
# Next is to inherit orders, minitems and maxitems
#' @inheritParams SVsByOrders   
#' @returns A vector or matrix, depending on which `vfunc` returns, with the
#'   approximate Shapley Values for each item/variable in the orderings and one
#'   column, if a matrix, for each measure returned by `vfunc`.
#' @returns A vector (typically) or matrix (if `both`, or if `ndep>1`) 
#'   of R-squared (or adjusted) increments. 
#'   If a matrix, columns correspond to the dependent variables, and 
#'   if `both`, further columns for adjusted r-squareds, in the same order. 
#' @export
SVsForKDRsByOrders <- function ( orders, scpMatrix, ndep=0, adjusted=NULL, 
                                 both=FALSE, minitems=0, maxitems=0 )  {
  # Set defaults and figure out what we are doing
  if ( is.null(maxitems) )  maxitems <- 0
  p <- nrow(scpMatrix)        # total number of variables, including dependents
  maxlen <- nrow(orders)      # length of each ordering
  nitems <- max(orders)       # number of items
  if ( ndep == 0 )  ndep <- p - maxlen    # deduce # of dependents
  if ( maxitems > 0 )  {
    if ( maxitems != maxlen )  {     # Not already cut down
      if ( maxlen != nitems )   {    # Not full either?
        stop ( "SVsForKDRsByOrders: Length of each ordering is neither ",
               "number of variables nor maxitems." )
      } else {
        orders <- orders[1:maxitems,,drop=FALSE]   # cut it down from full
        maxlen <- nrow(orders)
      }
    }
  }
  
  maxvars <- nrow(scpMatrix) - ndep
  tt <- adjustOrNot ( adjusted, both, scpMatrix, maxvars )
  doadj <- tt$doadj                             # handle adjusted/both stuff
  doraw <- tt$doraw
  totdf <- tt$totdf
  
  nres <- if ( doadj && both )  1  else  2;     # number of results/dependent
  kadj <- if ( doadj && both )  0  else  ndep;  # offset btwn raw and adj
  nres <- nres * ndep;             # total result columns, nres per dependent
  
  # ===== Get the real work of regressions done
  rsq <- KDRsolveOrders_cpp ( orders, scpMatrix, ndep )
  
  # Calculate adjusteds if needed
  if ( doadj )  {
    adjrsq <- array(NA,dim(rsq))
    for ( kdf in 1:NROW(rsq) )  {  
      adjrsq[kdf,,] <- 1 - ( 1 - rsq[kdf,,] ) * ( totdf / ( totdf - kdf ) )
    }
  }
  
  # Turn totals into increments
  for ( i in 1:ndep )  {
    if ( !doadj || both )  rsq[,,i] <- increments ( rsq[,,i,drop=FALSE] )   
    if ( doadj )  adjrsq[,,i] <- increments ( adjrsq[,,i,drop=FALSE] )
  }
  
  # Apply minimum items, if any
  if ( minitems > 0 )  {                   # drop items below the minimum
    drops <- - (1:minitems)
    rsq <- rsq[drops,,,drop=FALSE]         # must do even if not needed due to
    #   counts logic below!
    if ( doadj )  adjrsq <- adjrsq[drops,,,drop=FALSE]   
    orders <- orders[drops,,drop=FALSE]    # in case min=max, don't lose a dim
  }
  
  # Add up totals; get means, for straight r-squared.
  # First get counts for denominator.  Same for r-sq and adjusted, same all
  #   dependent variables.
  cdummy <- array ( 1, c ( dim(rsq)[1:2], 1 ) )
  counts <- totalup ( cdummy, orders )[,1]
  if ( doadj )  svadj <- totalup ( adjrsq, orders ) / counts
  if ( !doadj || both )  svs <- totalup ( rsq, orders ) /  counts
  # SVs for straight r-squared
  
  if ( both )  { svs <- cbind ( svs, svadj )
  } else if ( doadj )   svs <- svadj 
  if ( dim(svs)[2] == 1 )  svs[,1]  else svs
}

