# Utility functions internal to ShapCOA, mostly for checking of common 
#   parameters

#' Check whether respondent weights are OK
#' 
#' This function creates an error if weights are not of right length 
#'   or any are negative.
#'   
#' @param weights The weights vector to be checked.
#' @param nresp Number of respondents and expected number of weights.
#' @param caller Name of function calling us, for use in messages/errors.
#' @returns `TRUE` if weights look OK, `FALSE` if any problem found.
#' @keywords internal
checkWeights <- function ( weights, nresp, caller="??" )  {
  if ( length(weights) != nresp )  {     # wrong length?
    stop ( paste0 ( caller, ": Data has ", nresp, " rows but only ", 
                   length(weights), " weights were given.\n" ) )
  }
  if ( any(weights<0) )  {
    stop ( paste0 ( caller, ": Negative weights are not allowed.  ",
           "See position ", which(weights<0)[1], 
           " in weights, for example.\n" ) )
  }
  TRUE
}

#' Check validity of, and default, `minitems` and `maxitems` parameters to 
#' drivers
#' 
#' This function checks the legitimacy of `minitems` and `maxitems` parameters
#'   provided to SV-calculation drivers.  It also defaults them as 
#'   appropriate and returns their final values.
#'   The intention is to make handling of these parameters consistent 
#'     across drivers, while avoiding duplication of code.
#'   Invalid values produce a "stop" with a message.
#'   
#' Ideally, `nitems` will be provided by the caller.  
#' 
#' An alternative is to provide a data matrix, SCP matrix or orderings (COA)
#'   matrix from which the total number of items can be deduced.  
#'   Specifically:
#'   Usually either a raw data array with as many columns as items
#'   (assumed if there are more rows than columns),
#'   or a sum-of-crossproducts matrix with as many rows and columns as items
#'     plus one (assumed if there are as many rows as columns; 
#'     note that this will not work if there are multiple dependent variables
#'     in a key drivers regression), 
#'   or an orderings array (such as a COA) with either as many rows as items
#'     or as many rows as the maxitems value desired
#'     (assumed if it has more columns than rows).
#'   
#' @param minitems Minimum combination size to use for 
#'   size-limited Shapley Values.
#' @param maxitems Minimum combination size to use for 
#'   size-limited Shapley Values.
#' @param orders An orderings matrix, data matrix or SCP matrix from which
#'   the number of items can be deduced.
#' @param nitems The number of items or variables in the problem or dataset;
#'   must be specified by name if `orders` is omitted.
#' @param caller Name of function calling us, for use in messages/errors.
#' @returns A 2-item numeric vector with the final values of 
#'   minitems and maxitems.
#' @keywords internal 
ckMinMax <- function ( minitems=NULL, maxitems=NULL, 
                       orders=NULL, nitems=NULL, caller="??" )  {
  # Easy one first!
  if ( is.null(minitems) )  minitems <- 0    # Nothing dropped on low end
  
  # See what we can deduce from the orderings.  maxlen is implied # of items
  maxlen <- NA                               # we know nothing yet
  if ( !is.null(orders) )  {                 # orders implies things ...
    orddim <- dim(orders)
    if ( !is.null(orddim) )  {
      if ( orddim[1] == orddim[2] )  { maxlen <- orddim[1] - 1   # SCP matrix
      } else if ( orddim[1] > orddim[2] )  { maxlen <- orddim[2] # data matrix
      } else {                             # orderings matrix
        maxlen <- orddim[1]
        if ( is.null(nitems) )  nitems <- max(orders)
      }
    } else {
      stop ( caller, ": cannot recognize 'orders' parameter", call.=FALSE )
    }
    if ( !is.null(nitems) && nitems < maxlen )  {     # Nonsense!
      stop ( caller, ": Orderings given are longer than ",
             "number of columns of utilities" )    # implied maxitems
    }   
    if ( is.null(nitems) && !is.na(maxlen) )  nitems <- maxlen
  }
  
  # Default the max.  Usually all, but short orderings imply a maximum too.
  if ( !is.null(nitems) )  {                 # we know how many in total
    if ( is.null(maxitems) )  maxitems <- 0  # first-round default
    if ( maxitems == 0 )  {                  # now 
      maxitems <- min ( maxlen, nitems, na.rm=TRUE )  
                                             # implied max if short orderings
      if ( maxitems != nitems )  {           # if not the obvious answer
        message ( caller, ": maxitems=", maxitems, " implied by dimension ",
                "of orders parameter.\n",
                " Size-limited SVs will be done.\n" )  # Make sure user knows
      }     # Note this only happens if we determined maxitems here
    }                                        
    if ( !is.null(orders) )  {               # check max vs. orderings
      ok <- maxlen==nitems || maxlen==maxitems  # orders either all or pre-cut
      if ( !ok )  {
        stop ( caller, ": Length of each ordering ", maxlen, " is neither ",
               "number of items ", nitems, " nor maxitems ", maxitems,
               call.=FALSE )
      }
    }
    if ( maxitems > nitems )  {
      stop ( caller, ": maxitems ", maxitems, " exceeds total items ", nitems,
             call.=FALSE )
    }
    if ( minitems >= nitems )  {
      stop ( caller, ": minitems ", minitems, " must be strictly less than ",
            "total items ", nitems, call.=FALSE )
    }
  }
  
    
  # Stuff we can check even if nitems is unknown ...
  if ( minitems > maxitems )  {
    stop ( caller, ": minitems ", minitems, " exceeds maxitems ", maxitems,
           call.=FALSE )
  }
  if ( minitems < 0 || maxitems < 1 || minitems%%1!=0 || maxitems%%1!=0 )  {
    stop ( caller, ": negative or non-integer numbers of items in ",
           "minitems ", minitems, ", maxitems ", maxitems, call.=FALSE )
  }
    
  # Send fixed up values back to caller
  return ( c(minitems,maxitems) )   
}
