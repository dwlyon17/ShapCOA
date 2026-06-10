# Generate orders of entry, not with COAs, but randomly, cyclically, 
#   Williams-ly or for all possible.

#' Generate random orderings
#' 
#' This function generates random orderings, as many as desired, for as many 
#' variables as desired.  
#' 
#' This does not prevent generation of duplicates, which should be an 
#' issue only for large numbers of orders on few variables.  It is intended
#' for use in researching COA performance vs. random, not recommended for
#' "live" data use.

#' @param nvars Number of variables for which to generate orderings, which will 
#'   be the length of each ordering.
#' @param norders Number of orderings to generate.
#' @returns Array of dimension `nvars` by `norders`, with one ordering 
#'   per column.
#' @importFrom stats runif
#' @export 
genRandOrders <- function ( nvars, norders )  {
  # Generate norders purely random orders for nvars variables
  orders <- matrix ( runif(nvars*norders), nvars, norders )
  for ( i in 1:norders )  orders[,i] <- order ( orders[,i] )
  orders
}

#' Generate cyclic orderings
#' 
#' This function generates cyclic orderings, as many as desired, for as many 
#' variables as desired.  
#' 
#' This generates only full cycles, so the number of orderings may slightly 
#' exceed the number requested by `norders` when `norders/nvars` is not an 
#' integer.  
#' 
#' When multiple cycles are requested, each is random, with no attempt to 
#' prevent duplicates.  This is intended for use in research and play, not 
#' recommended for "live" data use.

#' @inheritParams genRandOrders
#' @param norders Number of orderings (minimum) to generate.
#' @param ncycles Number of full cycles to generate; overrides `norders` if
#'   not `NULL`.
#' @param random If `TRUE`, first set of cycles should be random, not `1:n`.
#' @returns Array of dimension `nvars` by (at least) `norders`, 
#'   with one ordering per column.
#' @importFrom stats runif
#' @export
genCyclic <- function ( nvars, norders=1, ncycles=NULL, random=FALSE )  {
  if ( is.null(ncycles) )  ncycles <- ceiling ( norders / nvars )
  cyc1 <- matrix ( 1:nvars, nvars, nvars, byrow=TRUE )
  cycle <- c(nvars,1:(nvars-1))              # mapping to cycle last to first
  for ( i in 2:nvars )  {
    cyc1[i,] <- cycle[cyc1[i-1,]]
  }                                          # now have one cyclic design
  cycles <- cyc1
  if ( ncycles > 1 || random )  {
    for ( i in 2:(ncycles+random) )  {
      remap <- order ( runif ( nvars ) )     # new ordering for remap
      cycles <- rbind ( cycles, matrix ( remap[cyc1], nvars, nvars ) )
    }
    if ( random )  cycles <- cycles[-(1:nvars),]   # toss first non-random one
  }
  t(cycles)                                  # transpose variables into rows
}

#' Generate Williams designs of orderings
#' 
#' This function generates Williams orderings, as many as desired, for as 
#' many variables as desired.  
#' 
#' Williams designs are a cyclic set of orderings plus its "fold-over" or
#' reverse orderings.  Each design has 2 x nvars orderings.
#' 
#' The function generates only full cycles, so the number of orderings may 
#' slightly exceed the number requested by norders when norders/nvars/2 
#' is not an integer.  
#' 
#' When  multiple cycles are requested, each is random, with no 
#' attempt to prevent duplicates.  This is intended for use in research 
#' and play, not recommended for "live" data use.
#' 
#' @inheritParams genRandOrders
#' @inheritParams genCyclic 
#' @inheritDotParams genCyclic
#' @returns Array of dimension `nvars` by (at least) `norders`, with one 
#'   ordering per column.  
#'   If `ncycles` is specified, there will be `nvars * 2 * ncycles` columns.
#' @export
genWilliams <- function ( nvars, norders=NULL, ... )  {
  # Williams idea and name discovered and suggested by Keith Chrzan.
  if ( !is.null(norders) )  norders <- ceiling(norders/2)     
                                    # since we double, need only 1/2 to start
  cycl <- genCyclic ( nvars, norders=norders, ... )
  cbind ( cycl, cycl[nvars:1,] )    # reverse order of rows ...
}

#' Generate all possible orderings
#' 
#' This function generates all possible orderings of as many variables 
#' as specified.  It is really just a wrapper for [combinat::permn()], as
#' "orderings" are permutations.
#' 
#' The function warns if called with `nvars>10` but attempts it anyway.  
#' Memory problems are likely in that case.  
#' 
#' While useful in computing exact Shapley Values, using all orderings is 
#' much less efficient than a combinations-based approach, and not 
#' feasible in problems with more than about 12 items or variables.

#' @inheritParams genRandOrders
#' @returns Array of dimension `nvars` by `nvars!`, with one ordering 
#'   per column.
#' @importFrom combinat permn
#' @export
genAllOrders <- function ( nvars )  {
  # Just call package function, convert to array and transpose as wanted.
  if ( nvars > 10 )  {
    cat ( "Warning: listing all permutations for nvars > 10", 
          " is not reasonable.\n" )
    cat ( "Trying anyway ...\n" )
  }
  matrix ( unlist(combinat::permn(nvars)), nrow=nvars ) 
}
