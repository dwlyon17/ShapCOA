# KDR solvers and general-purpose drivers for combination-approach Shapley
# Values, plus utilities to build CP matrices

#' Reconcile KDR adjusted and both flags, set totdf
#'
#' Internal function for setting up the KDR adjusted and both flags, and
#'   total degrees of freedom for adjusted r-squared.
#'
#' @param adjusted If numeric and not zero, requests the adjusted r-squared
#'   and gives the total d.f. to use for it.  
#'   If TRUE, requests adjusted r-squared using d.f. from the `ESS` attribute
#'   of the `scpMatrix` argument.  
#' @param both If `TRUE`, compute raw r-squared in addition to adjusted.
#' @param scpMatrix The SCP matrix to be used in the KDRs. 
#'   If `adjusted == TRUE`, it must contain 
#'   an `ESS` attribute specifying the total number of respondents or 
#'   effective sample size used to create it.
#' @param maxvars Total number of independent variables/items 
#'   (which the total d.f. must exceed).
#' @param caller Name of function calling us, for use in messages/errors.
#' @returns A list with three elements: doadj (logical, whether to compute
#'   adjusted r-squared), doraw (logical, whether to do straight r-squared),
#'   and totdf (total d.f. to use for adjusted r-squared, with 1 already
#'   subtracted for the constant term).
#' @keywords internal
adjustOrNot <- function ( adjusted=FALSE, both=FALSE, 
                          scpMatrix=NULL, maxvars=0, caller="??" )  {
  # Rigamarole for adjusted or not, and different ways of conveying
  #   the total degrees of freedom to use for adjusting.  
  if ( is.null(adjusted) )  adjusted <- both            # both=>adjusted
  if ( !adjusted && both )  {
    stop ( caller, "adjusted FALSE, both TRUE, will not cut it!" )
  }
  # adjusted now logical or numeric (not NULL)
  totdf <- 0                        # assume not adjusting
  if ( is.logical(adjusted) )  {
    if ( adjusted )  totdf <- attr(scpMatrix,"ESS")   
                    # total DF is effective sample size from the CP matrix
  } else {
    totdf <- adjusted        # if not logical, caller gave us the total DF
  }
  if ( adjusted )  {              # check for reasonable total d.f.
    totdf <- totdf - 1            # for convenience to callers
    if ( maxvars > 0 && totdf <= maxvars ) {    # This is just wrong!!
      stop ( caller, "adjusted must be > number of Xs",
             "(residual df for adjustment)" )
    }
  }
  return ( list ( doadj=adjusted, doraw=(!adjusted||both), totdf=totdf ) )
}

#' Construct a cross-product matrix for use in various key driver routines
#' 
#' This function builds the cross-product matrix needed for key driver 
#'   regressions.  
#'   It sets up for multiple dependent variables (always the
#'   last rows/columns of the result) and can either include a constant 
#'   row/column (always as the first row/column) 
#'   or pre-center all variables to avoid the need for one.
#'   
#' If `df` contains any non-numeric columns, they are ignored.  If there
#'   are constant columns, only the first is used.
#' 
#' The matrix produced has an `ESS` attribute giving Kish's effective sample 
#'   size (square of sum of weights, over sum of squares of weights) that can 
#'   be used as the total degrees of freedom when computing adjusted r-squared.
#'   
#' @param df A matrix or data frame of the X's and Y's for the result matrix.  
#' @param depvar Numeric vector of dependent variable columns in `df` OR
#'   character vector of dependent variable names.  If `NULL`, the
#'   final column of `df` is assumed to be the only dependent variable.
#' @param precenter Logical, whether to pre-center all variables and omit any
#'   constant row/column from the result.
#' @param silent Logical, whether to avoid sending commentary to console as 
#'   the CP matrix is constructed.
#' @param weights Vector of respondent weights, one per row of `df`,
#'   or character name of the weights column in `df`.  
#'   Unweighted data if `NULL`.  
#'   Weights may be zero to subset data but may never be negative.
#' @param missing If "none", any missing data in `df` is an error.
#'   If "listwise", any row/observation/respondent with any missing data
#'     is dropped from the CP matrix.
#'   If "pairwise", missing data is excluded only from the specific pairs
#'     of variables that include missings.  This can sometimes result in a
#'     cross-products or covariance matrix that is not positive definite
#'     and is thus not valid for regression, 
#'     as diagnosable by a negative eigenvalue for the SCP matrix,
#'     in which case an error is thrown.
#' @returns A square, symmetric cross-products matrix.  If `precenter` is 
#'   `FALSE`, the first row/column is for a constant variable, whether or not 
#'   one was found in `df`.  The dependent variables are always the last 
#'   row/columns.  Returns `NULL` if `weights` or `depvar` are invalid.
#' @importFrom stats complete.cases
#' @importFrom matrixStats colVars
#' @importFrom configural wt_cov
#' @export
cpBuild <- function ( df, depvar=NULL, precenter=TRUE, 
                      weights=NULL, silent=FALSE, 
                      missing=c("none","listwise","pairwise") )  {
  # Much inefficiency here, converting matrix to data.frame and back, but
  #   we are usually just called once and this allows flexibility in how
  #   dependents and weights are specified, and in which way df is structured
  #   in the first place.
  df <- data.frame ( df )      # In case matrix, adds Xn column names if needed
  
  # Convert numeric column parameters to column names 
  weights1 <- weights          # Remember original for error messages
  dep1     <- depvar           # Remember original for error messages
  if ( length(weights) == 1 && is.numeric(weights) )  {  
    weights <- colnames(df)[weights]      # wts as col # to name
  }
  if ( is.numeric(depvar) )  depvar <- colnames(df)[depvar] 
                                          # dependents as col #s to names
  
  # ==== Drop any non-numerics from the input
  drop <- rep(FALSE,ncol(df))
  for ( i in 1:ncol(df) )  drop[i] <- !is.numeric(df[[i]] )  
  if ( any(drop) && !silent )  {
    cat ( paste ( "cpBuild: Ignoring", colnames(df)[drop], "in input.\n" ) )
  }
  df <- df[,!drop]
  
  dfm <- as.matrix ( df )             # Matrix, not DF, for computing
  
  # ==== Handle the weights.  Can be given as numerics, or as a column name
  #      or number in the data frame.
  
  if ( is.null(weights) )  { weights <- rep(1,nrow(dfm)) # unweighted
  } else if ( length(weights) == 1 )  {    # weights are in the data frame
    wtpos <- match ( weights, colnames(dfm) )  # find it
    if ( is.na(wtpos) )  {                 # invalid/unfound weights
      cat ( "cpBuild: Invalid weight column", weights1, "\n" )
      browser()
      return (NULL)
    }
    weights <- dfm[,wtpos]                  # extract weights
    dfm <- dfm[,-wtpos]                     # remove from df
  }                              # not single, not character, must be vector
  
  if ( !checkWeights ( weights, nrow(dfm), "cpBuild" ) )  return (NULL)
                                            # ckWts will actually stop if bad
  effsize <- sum(weights)^2 / sum(weights^2)   # Kish's effective sample size
                # will be attached to result for use in adjusting r-squareds
  
  # ==== Find any constant column(s) in the df.  
  
  cons <- which ( matrixStats::colVars(dfm) == 0 )  # constant columns
  if ( length(cons) > 0 )  {
    consval <- colMeans(dfm[,cons,drop=FALSE])   # Hopefully, 1.0
    if ( any(consval != 1.0 ) )  {
      cat ( "cpBuild:  Warning: Some columns are constant (and not 1.0):",
            colnames(dfm)[cons[which(consval!=1.0)]], "\n" )
      cat ( "Building CP matrix anyway, but this is not good.\n" )
    }
  }
  
  cons <- which ( matrixStats::colVars(dfm) == 0 & colMeans(dfm) == 1.0 )
  if ( length(cons) > 0 )  {
    if ( length(cons) > 1 )  {
      warning ( "Multiple constant columns in data frame.  Using first." )
    }
    cname <- colnames(dfm)[cons[1]]
    if ( !precenter )  {                 # if we want a constant column
      if ( !silent )  cat(paste("cpBuild: Constant column is", cname, "\n" ))
      dfm <- cbind ( dfm[,cons[1]], dfm[,-cons] )   # make one first, drop all
      colnames(dfm)[1] <- cname
    }
  } else if ( !precenter )  {           # want constant, but none found
    if ( !silent )  cat ( "cpBuild: Adding constant column to input.\n" )
    dfm <- cbind ( Constant=1, dfm )    # add if missing
  }
  
  # ==== Handle dependent variable(s).  
  #      Move all y's to the very end.  (Can be more than one)
  
  if ( !is.null(depvar) )  {
    ys <- match ( depvar, colnames(dfm) )  
    if ( anyNA(ys) )  {
      cat ( "cpBuild: Dependent variable(s)", dep1[which(is.na(ys))], 
            "not in the data frame.\n" )
      return (NULL)
    }
  } else  {
    ys <- ncol(dfm)             # default: original last, before 1
    if ( !silent )  cat ( "cpBuild:", colnames(dfm)[ys], 
                          "taken as dependent variable.\n" )
  }
  ynames <- colnames(dfm)[ys]            # save the Y names
  dfm <- cbind ( dfm[,-ys,drop=FALSE], dfm[,ys,drop=FALSE] )   
                                         # move all y's to very end
#  colnames(dfm)[(ncol(dfm)+1-length(ys)):ncol(dfm)] <- ynames 
  
  # ==== Final wrap-up and actual CP matrix.
  
  if ( ncol(dfm) <= length(ys) + !precenter )  {   # Must be at least one X!
    stop ( paste0 ( "cpBuild: No independent variables ",
                    "(just dependent, weights and/or constants) found.\n" ) )
  }
  
  # Now we have dfm in the column order:  1 (maybe), Xs, Ys.
  
  # Finally ready for the easy part, except for missingness:  
  missing <- match.arg ( missing )        # clean up the parameter
  nonmiss <- complete.cases ( dfm )       # how many complete cases
  if ( any(!nonmiss) )  {                 # any missing data?
    if ( missing == "listwise" )  {
      message ( "cpBuild: Of ", nrow(dfm), " cases, ", sum(!nonmiss),
                " have missing data and are being dropped." )
      if ( sum(nonmiss) < ncol(df)+2 )  {
        stop ( "cpBuild: Too many cases dropped for this many variables." )
      }
      dfm <- dfm[nonmiss,]                # drop cases with missings
      weights <- weights[nonmiss]
      missing <- "none"                   # with cases dropped, just normal
    } else if ( missing == "none" )  {
      stop ( "cpBuild: missing='none' but there are ", 
             nrow(dfm)-sum(nonmiss), " cases with missing data." )
    }
  } else missing <- "none"                # if none missing, none is OK!
  
  switch ( missing, 
           none = {                       # none or fixed-up listwise
             if ( precenter )  {          # Pre-center if desired
               colmeans <- colSums(dfm*weights) / sum(weights)
               dfm <- t ( t(dfm) - colmeans )
             }
             dfm <- dfm * sqrt(weights)   # CPs will re-square the root!
             cp <- crossprod ( dfm )      # Just take the built-in crossproduct
             effsize <- sum(weights)^2 / sum(weights^2)   
                                          # Kish's effective sample size
             attr(cp,"ESS") <- effsize    # for us in adjusting r-squareds
           },
           
           pairwise={
             wtpres <- matrix ( weights, length(weights), ncol(dfm) )
                                          # weights per column/variable
             wtpres[is.na(dfm)] <- NA     # missing data -> missing weight
             effsize <- colSums ( wtpres, na.rm=TRUE ) ^ 2 /
               colSums ( wtpres^2, na.rm=TRUE )     # Kish's ESS per column
             if ( any(effsize<=0) || any(is.na(effsize) ) )  {
               stop ( "cpBuild: one variable is entirely missing!" )
             }
             effsizef <- exp ( mean ( log ( effsize ) ) )  # geo mean of values
             Ngood <- outer ( 1:ncol(dfm), 1:ncol(dfm), 
                              Vectorize ( function(i,j)
                                        sum(complete.cases(dfm[, c(i, j)])) ) )
             message ( "cpBuild: Using missing=pairwise, pairs have a maximum",
                       " of ", nrow(dfm)-min(Ngood), " missing cases\n",
                       "         ESS ranges from ", round(min(effsize),1), 
                       " to ", round(max(effsize),1), ", using ",
                       round(effsizef,1), " for any r-sq adjustment." )
             if ( precenter )  {
               colmeans <- colSums ( dfm*wtpres, na.rm=TRUE ) / 
                 colSums ( wtpres, na.rm=TRUE )  # wtd means for non-miss data
               dfm <- t ( t(dfm) - colmeans )     # NAs stay NA
             }
             cp <- wt_cov ( dfm, , weights, use="pairwise", unbiased=FALSE )
             eigvals <- eigen ( cp, symmetric=TRUE, only.values=TRUE )$values
             print ( min(eigvals) )
             if ( any(eigvals <= 0 ) )  {
               stop ( "cpBuild: missing=pairwise produced a covariance ",
                      "matrix with a negative eigenvalue(s),\n",
                      "         indicating incompatibility among ",
                      "pairwise estimates." ) 
             }
             attr(cp,"ESS") <- effsizef    # for us in adjusting r-squareds
           },
           
           { stop ( "cpBuild: invalid 'missing' parameter: ", missing ) 
           }
         )
  cp  
}

#' "Sweep out" and drop the constant row/column of a cross-products matrix
#' 
#' This function removes the constant from a cross-products matrix, preparing
#' it for use in various key driver regression situations.  Since the 
#' constant does not affect r-squared explained, its presence is simply a
#' nuisance when computing r-squareds and increments.
#' 
#' An easy alternative to calling this function is to use `precenter=TRUE` 
#'   with [cpBuild()], which eliminates any constant from the matrix.
#' 
#' If matrix `cp` has an `ESS` attribute giving the effective sample size 
#'   (see [cpBuild()]), it is preserved in the result. 
#' @param cp A cross-products matrix (square, symmetric).
#' @param constant The row/column number of the constant, often 1.
#' @returns A matrix with one less row and column than `cp`, with the 
#'   Gauss-Jordan sweep (or Beaton's `swp`) via [swp()] for the constant 
#'   having been done.  
#'   The result is still square, but no longer symmetric.
#' @export
cpConstantOut <- function ( cp, constant=1 )  {
  cpnew <- swp ( cp, constant ) [-constant,-constant]
  attr(cpnew,"ESS") <- attr(cp,"ESS")
  invisible(cpnew)
}

#' Do a Gauss-Jordan "sweep" (aka Beaton's `swp` operator) on a single 
#'   row/column of a cross-products matrix
#'
#' This function performs a single sweep or "swp" operation on a matrix.  
#'
#' This is a building block for full regressions, and for stepping through
#'   variables one at a time.  It is more for testing and play than for any
#'   real production use, as repeated calls over many orderings are very
#'   inefficient and various integrated C++ functions handle this better.
#'   
#' @param CP The cross-products matrix to be swept.  It may well have been
#'   partially swept before.
#' @param k The number of the row/column to be swept.
#' @returns The swept cross-products matrix, with same dimensions as `CP`.
#' @keywords internal
swp <- function ( CP, k )  {
  # For references, see matlib::swp and fastmatrix::sweep.operator help.
  pivot <- CP[k,k]
  newrow <- CP[k,] / pivot
  oldcol <- CP[,k]
  CP <- CP - tcrossprod ( oldcol, newrow ) 
  CP[k,] <-   newrow 
  CP[,k] <- - oldcol / pivot
  CP[k,k] <- 1 / pivot
  CP
}

#' Solve a linear KDR for one combination of variables and one dependent 
#'   variable
#' 
#' This function is a very simple (and not very efficient) illustration of 
#'   a valuation function.  It finds the linear regression r-squared for a
#'   single dependent variable regressed on a single subset/combination of 
#'   variables.  It is intended for use with [SVsByCombos()].
#' 
#' For practical use, see [KDRsolveCombos()].  Inspect the source of this function
#'   to see just how simple a valuation function can be.
#' 
#' It is more efficient to handle many combinations at once, so as to avoid
#'   repeated function calls.  If there are multiple dependent variables to be
#'   done, it is more efficient in a linear regression to handle them in 
#'   parallel as well.  But those efficiency issues are not large.  However,
#'   putting the key code in C++, as [KDRsolveCombos()] effectively does, 
#'   gives a huge efficiency gain.  
#'   For larger problems, the valuation function needs to be very 
#'   efficient, especially when used with the [SVsByCombos()] driver that gives
#'   exact Shapley Values.
#' 
#' @param combo An vector of variable numbers (rows/columns of `scpMatrix`)
#'   specifying the variables in the combination to be evaluated.
#' @param scpMatrix A sum-of-cross-products matrix, from which either 
#'   a) the constant term in the regression has already been swept out 
#'      (easily accomplished by pre-centering all variables before creating
#'      the SCP matrix) or b) the constant variable is included in `combo`.
#' @param depvar The row/column number in `scpMatrix` of the dependent variable.
#' @returns The (unadjusted) r-squared for the specified regression.
#' @export

KDRsolve1Combo <- function ( combo, scpMatrix, depvar )  {
  # A la relaimpo but vastly cut down!
  rsq <- scpMatrix[depvar,combo] %*% 
          solve ( scpMatrix[combo,combo], scpMatrix[combo,depvar] )  /
            scpMatrix[depvar,depvar]
  c(rsq)      # Get rid of (1,1) dimensioning
}

#' Solve slowly for linear r-squared for multiple combinations of variables
#' 
#' This function finds r-squareds of regressions on multiple combinations 
#' of one or more independent variables.  All combinations must be of the 
#' same size.  It is intended for use with [SVsByCombos()].
#' [KDRsolveCombos()] is faster, as it does the same things in C++.
#' 
#' @param combo An matrix of variable numbers (rows/columns of `scpMatrix`)
#'   specifying the variables in the combinations to be evaluated, 
#'   one column per combination,
#'   and with entries giving which rows/columns of the SCP matrix (next 
#'   parameter) the desired variables are in.  
#' @param depvar A vector of the row/column numbers in `scpMatrix` of the 
#'   dependent variables.
#' @inheritParams KDRsolveCombos1_cpp
#' @return If there are multiple dependent variables (`length(depvar)>1`), 
#'   an array of r-squared values, one row per column in `combo` and one column
#'   per dependent variable.  For a single dependent variable, a vector
#'   of r-squared values, one per combination in `combo`.
#' @export
KDRsolveComboR <- function ( combo, scpMatrix, depvar )  {  
  # A la relaimpo but vastly cut down!
  # This is around 10% (or less) faster than looped calls to 
  #   individual solutions.
  
  vstart <- scpMatrix[depvar,depvar]
  ncomb <- ncol(combo)            # number of combos to be done
  ndep <- length(depvar)          # number of dependent vars
  dvsel <- cbind(1:ndep,1:ndep)   # diag selector for dvs
  mdv <- ndep > 1                 # logical -- need to select result diagonal?
  rsq <- array ( 0, c ( ncomb, ndep ) )
  for ( icol in 1:ncomb )  {      # for each combination to be done
    ccomb <- combo[,icol]         # the combo for this iteration
    rsq[icol,] <- ( scpMatrix[depvar,ccomb] %*% 
        solve ( scpMatrix[ccomb,ccomb], scpMatrix[ccomb,depvar] )  / 
          vstart )  [dvsel]       # tested [dvsel] vs diag of, this is faster
  }
  if ( ndep == 1 )  rsq[,1]  else  rsq
}

#' (Adjusted) R-squared for linear regression on a set of combinations
#'
#' This is a valuation function for straight and/or adjusted r-squared for 
#'   linear regression on one or multiple dependent variables.  
#'   
#' This function is designed to 
#'   be called from [SVsByCombos()], and uses two 
#'   C++ functions to do the heavy lifting efficiently.
#'
#' Multiple combinations, all of the same size, are handled in a single call,
#'   so `multi=TRUE` should be specified to [SVsByCombos()] when it
#'   is told to use this function.
#'   
#' @param adjusted If not `NULL`, computes adjusted r-squared instead of 
#'   regular/straight r-squared.  
#'   Can be a numeric value giving the total degrees of freedom (DF) 
#'   to use for the adjustment, or can be `TRUE`,
#'   meaning to use the `ESS` (effective sample size) *attribute* of 
#'   `scpMatrix`, created by [cpBuild()].
#' @param both If `TRUE`, computes both straight and adjusted r-squareds.
#' @inheritParams KDRsolveCombos_cpp
#' @returns For a single dependent variable and `both=FALSE`, a vector of
#'   r-squared values, one per column of `combo`.  Otherwise, a matrix with as 
#'   many rows as `combo` has columns, and a column for each dependent variable,
#'   followed by, if `both` is `TRUE`, a column of adjusted r-squareds for
#'   each dependent variable.
#' @export
KDRsolveCombos <- function ( combo, scpMatrix, depvar, 
                             adjusted=NULL, both=FALSE )  {
  # Choose the better C++ routine; avoid overhead of multiple DVs if we can
  vfunc <- if (length(depvar)==1) KDRsolveCombos1_cpp  else KDRsolveCombos_cpp
  
  # Make sure depvar is numeric, not character
  if ( is.character(depvar) )  {
    depvarc <- depvar
    depvar <- match ( depvarc, colnames(scpMatrix) )
    if ( any(is.na(depvar)) )  {
      cat ( "KDRsolveCombos: Could not find dependent variable(s)", 
            depvarc[which(is.na(depvar))], "\n" )
      return (NULL)
    }
  }
  rsqs <- vfunc ( combo, scpMatrix, depvar )    # Basic computation
  
  maxvars <- nrow(scpMatrix) - length(depvar)
  tt <- adjustOrNot ( adjusted, both, scpMatrix, maxvars )
  doadj <- tt$doadj                             # handle adjusted/both stuff
  doraw <- tt$doraw
  totdf <- tt$totdf
  
  if ( doadj )  {
    adjrsqs <- 1.0 - ( 1.0 - rsqs ) * ( totdf / ( totdf - nrow(combo) ) )
    if ( both )  { rsqs <- cbind ( rsqs, adjrsqs )  
    } else         rsqs <- adjrsqs
  }
    rsqs
}
