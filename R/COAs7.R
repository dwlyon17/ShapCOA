# Make, and check, Galois field OAs (orthogonal arrays) and COAs (Component
#   Order of Addtiion) designs, including multipled COAs.

#' Check an orthogonal array to be sure it is OK (strength 2 only)
#' 
#' This function checks for equal pairwise precedence counts and for 
#' uncorrelated columns.
#' 
#' This is not meant to be called directly by a user (see [assessOA()] instead).
#' @param oades The orthogonal array to be checked, in 0-coded nrow>ncol format.
#' @param silent If `TRUE`, a message is printed only if the check fails,
#'   not if it succeeds.
#' @returns `TRUE` if the OA is valid, `FALSE` if not.
#' @importFrom stats cor
#' @keywords internal
ckOA <- function ( oades, silent=FALSE )  {
  # Fixed 10/8/25 to check equal pairwise counts, not just correlations.
  if ( !silent )  {
    cat ( "Checking OA of size ", dim(oades), min(oades), 
          "to", max(oades), "\n" )
  }
  ok = TRUE
  if ( NCOL(oades) <= 1 )  return (FALSE)
  for ( i in 1:(ncol(oades)-1) ) {
    tt <- length ( unique ( table ( oades[,i], useNA="ifany" ) ) )
    if ( tt > 1 )  stop ( "Column", i, "is not equi-distributed (balanced)." )
  }
  
  # Now check columns for correlations.  Report any (first 20) found.
  # It's a necessary but not sufficient condition to have zero correlations,
  # so this is not a necessary check.  But it's fast and spots problems
  # quickly.
  ccor <- cor ( oades )
  ccor[lower.tri(ccor,diag=TRUE)] <- NA     # drop diags and duplicates
  nz <- which ( ccor!=0, arr.ind=TRUE )
  if ( length(nz) > 0 )  {
    nzc <- cbind ( nz, round(ccor[nz],3) )
    colnames(nzc) <- c("col1","col2","corl")
    cat ( nrow(nzc), "columns are correlated" )
    if ( nrow(nzc) > 20 )  {
      cat ( ".  First 20 are:" )
      nzc <- nzc[1:50,]
    } else cat ( ":" )
    print ( t(nzc) )
    ok = FALSE
  }
  if ( ok && !silent )  cat ( "All column pairs are uncorrelated\n" )
  
  # Now check that pairwise values are equi-frequent between all possible pairs
  # of columns in the design.  This is the critical, defining, check.  If it
  # passes, then zero correlations are guaranteed.
  nvars <- ncol(oades)
  pcounts <- matrix(NA,nvars,nvars)
  mult <- 100000
  target <- ( nvars - 1 ) ^ 2
  nbad <- 0
  for ( i in 1:(nvars-1) )  {
    for ( j in (i+1):nvars )  {   # for all pairs of variables, 0:ncol-1
      if ( length ( unique ( oades[,i] + mult*oades[,j] ) ) != target )  {
        nbad <- nbad + 1
        if ( nbad < 10 )  { 
          cat ( "Columns", i, "and", j, "are not paired across all values.\n" )
        }
        ok <- FALSE
      }
    }
  }
  if ( nbad > 0 )  {
    cat ( nbad, "pairs of columns are not equally paired.\n" )
  } else if ( ok && !silent )  {
    cat ( "The OA checks out with all pairs of values equally",
          "present for all pairs of columns.\n" )
  }

  ok
}

#' Check a Component Order of Addition design for key properties
#' 
#' This function checks a COA design for equi-positioning and 
#' pairwise equi-precedence. 
#' 
#' This is not meant to be called directly by a user (see [assessCOA()] 
#' instead).
#' 
#' Various diagnostic messages are sent to the console in the event of failure.
#' 
#' COAs that have been size-limited (the orderings all shortened) or formed by
#'   dropping more than one item from a "true" COA will not pass all checks, 
#'   but are still acceptable for most use cases.
#' 
#' It is possible a design could pass
#'   both tests without being an actual COA.  In particular, "Williams"
#'   designs will pass.  However, they do not resemble COAs at all, and
#'   this function is only intended as a safety check on COA generation,
#' not to catch all possible frauds!
#' @param coades The COA array to be checked, in 0-coded nrow>ncol format
#' @inheritParams ckOA 
#' @returns `TRUE` if the COA passes all checks, `FALSE` if not.
#' @keywords internal
ckCOA <- function ( coades, silent=FALSE )  {
  nvars <- max(coades) + 1
  if ( !silent )  cat ( "Checking COA design of size ", dim(coades), "for", 
                        nvars,  "variables.\n" )
  ok <- TRUE
  
  # Be sure right values are in there.
  if ( !all.equal ( sort(unique(c(coades))), 0:(ncol(coades)-1) ) )  {
    cat ( "Design entries are not 0 to #vars minus 1 as expected.\n" )
    cat ( "Actual range is ", range(coades), ".\n" )
    ok <- FALSE
  }
  # Equipositioning test is easy:  constant sums for rows and columns.
  if ( length ( unique ( rowSums(coades) ) ) != 1 ) {
    nvars <- max(coades) + 1
    cat ( "Each item does not appear once and only once per ordering.\n" )
    if ( nvars > ncol(coades) )  {
      cat("(not surprising, since orderings are shorter than # of variables.\n")
    }
  }
  if ( length ( unique ( colSums(coades) ) ) != 1  )  {
    cat ( "Items do not appear with equal frequency in each position.\n" )
    ok <- FALSE
  }
  
  posn <- array ( NA, dim(coades) )  # position of each var in each row
  nvars <- ncol(coades)              # vars numbered 0 to ncol-1
  for ( k in 1:nrow(coades) )  posn[k,coades[k,]+1] <- 1:nvars
  
  if ( !any(is.na(posn)) )  {        # won't work if any NAs in here
    precounts <- matrix(NA,nvars,nvars)
    for ( i in 1:(nvars-1) )  {
      for ( j in (i+1):nvars )  {   # for all pairs of variables, 0:ncol-1
        precounts[i,j] <- sum ( posn[,i] < posn[,j] )
      }
    }
    upres <- unique ( c(precounts) )    # should have same value for all
    if ( length(upres) != 2 || upres[!is.na(upres)] != nrow(coades)/2 )  {
      cat ( "Not all pairs appear equally often before & after each other.\n" )
      cat ( "Pairwise counts matrix follows.\n" )
      ok2 <- ok <- FALSE
      print ( precounts )
    } else ok2 <- TRUE
    if ( ok2 & !silent )  {
      cat ( "The COA checks out with equal pairwise precedence.\n" )
    }
  } else cat ( "Not possible to check pairwise precedence.\n" )  # due to NAs
  ok 
}

#' Convert orthogonal array or COA design to internal "0-based nrow>ncol" format
#' 
#' This converts an orderings or combinations design 
#' * to use 0-based numbering of items, 
#' * to eliminate any constant 0 or 1 entry in all first positions, and 
#' * to have runs/orders/combinations in the rows and items in the columns
#'   (meaning the number of rows almost always exceed the number of columns).
#'   
#' The design is examined to see which of the output criteria is does not 
#'   already meet.  An already-compliant design is returned unchanged.
#'
#' Aside from development and testing, it is not clear what the use case for 
#'   this would be!
#'
#' The possibility of a leading 0 or 1 arises from designs used for logistic
#'   regression, where the constant must always be entered first and 
#'   explicitly.
#' 
#' @param des A design or orthogonal array (numeric matrix) to be converted. 
#' @returns The design with 0-based coding, with any constant first column 
#'   of ones removed, and with runs or orders in the rows, and items in the 
#'   columns (nrow>ncol format).
#' @seealso [userCode()] to do the reverse.
#'
#' @export
designCode <- function ( des )  {
  if ( min(des) == 1 )  des <- des - 1    # back to 0-based if necessary
  if ( ncol(des) > nrow(des) )  des <- t(des)  # transpose if needed
  if ( all ( des[1,] == 1 ) )  des <- des[-1,] - 1L
                                          # if constant forced in a la logistic
  des
}

#' Convert orthogonal array or COA design to 1-based ncol>nrow format
#' 
#' This converts an orderings or combinations design to a format equivalent
#' to what standard R function `combn` produces, and the format the functions
#' using COAs expect.  That format:
#' * uses 1-based numbering of items, 
#' * with runs/orders/combinations in the columns and items in the rows
#'   (meaning the number of columns almost always exceeds the number of rows).
#'   
#' The design is examined to see which of the output criteria is does not 
#' already meet.  An already-compliant design is returned unchanged.
#'
#' Aside from development and testing, it is not clear what the use case for 
#' this would be!
#' 
#' @param des A design or orthogonal array (numeric matrix) to be converted. 
#' @returns The design with 0-based coding, 
#'   and with runs or orders in the rows, and items in the 
#'   columns (nrow>ncol format).
#' @seealso [designCode()] to do the reverse.
#' @export
userCode <- function ( des )  {
  if ( min(des) == 0 )  des <- des + 1L        # back to 0-based if necessary
  if ( nrow(des) > ncol(des) )  des <- t(des)  # transpose if needed
  des
}

#' Assess an orthogonal array for validity at strength 2
#' 
#' This function is a friendly wrapper for [ckOA()].  It accepts an OA in 
#' either 0-based nrow>ncol format or 1-based ncol>nrow format.
#' @param oades The orthogonal array to be assessed.
#' @inheritDotParams ckOA silent
#' @returns `TRUE` if the OA is valid, `FALSE` if not.  
#' @seealso `ckOA` for details of the checks that are made.
#' @export
assessOA <- function ( oades, ... )  {
  oades <- designCode ( oades )
  ckOA ( oades, ... )
}

#' Assess a COA design for the usual COA properties
#' 
#' This function is a friendly wrapper for [ckCOA()].  It accepts a 
#' purported COA in 
#' either 0-based nrow>ncol format or 1-based ncol>nrow format.
#' @param coades The design to be assessed.
#' @inheritDotParams ckCOA silent
#' @returns `TRUE` if the design passes all checks, `FALSE` if not.  
#' @seealso `ckCOA` for details of the checks that are made.
#' @export
assessCOA <- function ( coades, ... )  {
  # "Assess" a COA design, which is not necessarily a pure COA.
  coades <- designCode ( coades )
  ckCOA ( coades, ... )
}

#' Find smallest prime or power of a prime of at least a given size
#'
#' This function finds the smallest size of an orthogonal array (OA) or 
#' Component Order of Addition (COA) design that is large enough
#' for a given number of items or variables.  A valid size must be a prime
#' number or a power of a prime number.
#'
#' This function can be used to pre-determine the size of OA or COA that will
#' be needed for a given problem, or to generate a list of sizes for 
#' research, play or presentation purposes.
#' @param nvars The number of items/variables needed.
#' @param all if `TRUE`, returns all primes/powers up to `nvars`.
#' @returns For `all` `FALSE`, an integer prime or power.  For `all` `TRUE`, 
#'   an integer vector of primes and powers, all less than `nvars`.`
#' @importFrom numbers Primes
#' @export
primeSize <- function ( nvars, all=FALSE )  {
  upper <- ceiling ( nvars * 1.5 )         # highest we will search
  primes <- numbers::Primes ( 2, upper )  # all the primes
  bases <- primes[primes<=upper]    # ones for ^2 and higher
  powers <- NULL
  for ( b in bases )  {
    maxp <- floor ( log(upper) / log(b) )   # maximum power to generate
    if ( maxp > 1 )  powers <- c ( powers, b^(2:maxp) )
  }
  allp <- sort ( c ( primes, powers ) )    # combine primes and powers
  if ( all )  allp[allp<=nvars]   else   min ( allp[allp>=nvars] ) 
                               # smallest one greater unless want all of them
}

#' Generate smallest orthogonal array for a given number of items/variables
#' 
#' This function generates the smallest Galois orthogonal array large enough 
#'   to handle a given number of items or variables.
#' 
#' @param nocheck Whether to skip running [ckOA()] on the generated OA.  
#'   Running the check can be time-consuming and is only needed in early 
#'   debugging.
#' @param silent Logical, passed on to [ckOA()] if `nocheck` is FALSE
#' @param random Logical, item numbers in the OA are randomly relabeled if 
#'   `TRUE`.  This relabeling does *not* make the resulting OA or any COAs 
#'   derived from it different in any practically useful way, and is *not* a
#'   substitute for [permuteCOA()].
#' @inheritParams primeSize
#' @returns An orthogonal array, with items numbered 0 to `nvars-1`, 
#'   `nvars` or slightly more (depending on next-largest prime or power) 
#'   columns and as many rows as the square of the number of columns.  
#'   This is the "0-based, nrow>ncol" format for an OA.
#' @importFrom lhs createBusht
#' @export
getOA <- function ( nvars, nocheck=FALSE, silent=FALSE, random=FALSE )  {
  nvar2 <- primeSize ( nvars )   # adjust size upward to one that works
  oades <- lhs::createBusht ( nvar2, nvar2+1, strength=2, bRandom=random )
                             # will be for levels 0:(nvars-1), with nvars+1
                             # columns, strength=2 -> pairwise independent
  # Non-silence and checking were mostly for debugging.  Checking is by far
  #   the slowest part of the process.
  if ( !silent )  {
    cat ( "Requested design for", nvars, if (nvars==nvar2) "and"  else "but", 
          "generated for", nvar2, "\n" )
  }
  if ( !nocheck )  ckOA ( oades )        # Be safe/sure
  invisible(oades)
}

#' Generate a single "true" Component Orthogonal Array (COA) 
#' 
#' This function generates a singleton pure COA, either from scratch or from
#' a supplied orthogonal array.
#' See Yang, Sun and Xu 2021 (Technometrics May 2021, 63:2).
#' 
#' The COA generated is exactly according to the Yang, Sun and Xu algorithm.
#' It may be larger than requested if the requested size is not a prime or
#' power of a prime, does *not* have columns permuted to alleviate
#' triple imbalance and is not checked by [ckCOA()].
#' 
#' This function is for internal use, only, by callers who handle the drops, 
#' permutations and checking.
#' @param oades Either the number of items/variables for which a COA is desired 
#'   OR an orthogonal array to be used, in the 0-based nrow>ncol format.
#' @inheritDotParams getOA nocheck silent random 
#' @returns A Component Order of Addition (COA) design, in the 0-based 
#'   nrow>ncol format, whose size matches that of the supplied OA in `oades` or
#'   is at least as large as the number specified by `oades`.
#' @keywords internal
makeCOA <- function ( oades, ... )  {
  if ( length(oades) == 1 )  oades <- getOA ( oades, ... )  
                                      # Get an OA of given size
  nc <- ncol(oades)
  size <- sqrt ( nrow(oades) )
  
  if ( size+1 != nc )  {
    stop ( "Dimensions of OA for makeCOA are ", nrow(oades), " ", nc, 
           ", not as expected." )
  }
  if ( min(oades) != 0 || max(oades) != size-1 )  {
    stop ( "Level numbers of OA for makeCOA are ", min(oades), " to ", 
           max(oades), ", not as expected." )
  }  
      
  # Enough checking already, here's the actual algorithm.
  oades <- oades [ order ( oades[,1], oades[,2] ), ]    # sort OA by first two
  maps <- oades[1:size,-1]       # first n rows control mapping of rest of cols
  coades <- oades[-(1:size),-1]  # remainder will be the COA after mapping
                                 # note that control col gets dropped
  
  invperm <- rep(NA,size)        # size the inverse permutation vector to start
  for ( j in 1:(nc-1) )  {       # all remaining columns in maps & oades 
    invperm[maps[,j]+1] <- 0:(size-1)     # how map top of column to 0:n-1
    coades[,j] <- invperm [coades[,j]+1]  # note +1 since start w/ 0:n-1
  }                                       # remap the column 
  coades
}

#' Drop excess items from a COA
#' 
#' This function removes extra items from a COA when the generated COA 
#' (typically from `makeCOA`) has more items than were desired.
#' 
#' This is not meant for direct calls by the user; it uses the 0-based 
#' nrow>ncol format and will not work on COAs delivered to the user by 
#' `getCOA`.
#'
#' When COAs are permuted to make multiple COAs, the permutation should be 
#'   done *before* any dropping of items, not after.  
#' 
#' @param coades The COA from which items are to be dropped, in 
#'   0-coded nrow>ncol format.
#' @param nvars The number of items/variables to keep.
#' @returns A new, reduced, Component Order of Addition (COA) design, 
#'   in the 0-based nrow>ncol format, with items numbered nvars or higher 
#'   (0-based) removed.  This result is not a "true" COA; it will have
#'   pairwise equi-precedence, but will not have equi-positioning if more 
#'   than one item was removed.  The result has the same number of rows as
#'   the original, but only `nvars` columns.
#' @export
dropCOA <- function ( coades, nvars )  {
  # Line below is for runs in rows, posns in columns, 0-based numbering
  coades2 <- matrix ( t(coades)[t(coades)<nvars], nrow(coades), nvars, 
                      byrow=TRUE )
  
  # Line below would be for 1-based ncol>nrow format
  #  coades2 <- matrix ( coades[coades<=nvars], nvars, ncol(coades) )
}

#' Permute columns of a COA to generate a new, non-overlapping, COA
#' 
#' This function permutes a given COA to create a new one, none of whose
#' orderings are in the original.  It can produce as many permuted versions
#' as desired, with no overlaps among any of them.  Each permuted version has
#' the same equi-positioning and equi-precedence properties as the original.
#' 
#' The new permutatoins are randomly generated on each call and thus are *not*
#' reproducible across runs unless the caller manages the random `set.seed`.
#'
#' If a single new permutation is wanted, use `ncoa=2` and drop the first 
#' COA returned (which is the original `coades`).
#'
#' When COAs are permuted to make multiple COAs, the permutation should be 
#'   done *before* any dropping of items, not after.  This function complains
#'   if asked to permute a COA that appears to have had items dropped.
#' 
#' @param coades The COA from which permutations are to be generated, in 
#'   0-coded nrow>ncol format.
#' @param ncoa The number of permutations to produce.  If less than 2, 
#' the original `coades` is returned as the result.
#' @returns A list, each element of which is a COA, with the first element in
#'   the list being the original `coades`.  
#'   The COAs are in the 0-based nrow>ncol format.  
#' @importFrom combinat permn
#' @importFrom stats runif
#' @export
permuteCOA <- function ( coades, ncoa=2 )  {
  if ( ncoa < 2 )  return ( list(coades) ) # No point if just want one (or none)
  
  ncols <- ncol ( coades )                
  ncheck <- primeSize ( ncols )         
  if ( ncols != ncheck )  {            # See if ncols is prime or power
    cat ( "Warning: permuteCOA is asked to permute a COA for ", ncols, 
          "items.\nThis appears to have had items dropped from a larger COA.\n",
          "Permutation should always be done BEFORE dropping items.\n",
          "Proceeding anyway.  Use caution!\n" )
  }
  # can be >nvars if nvars not prime^n[-1]
  if ( ncols < 16 )  {                  # beyond this, numbers are HUGE!
    maxperm <- factorial ( ncols - 2 )    # max permutations we can do
    if ( ncoa > maxperm )  {
      stop ( paste ( "permuteCOA: can make at most", factorial(ncols-2), 
                     "different COAs of size", ncols, "but", ncoa, 
                     "were requested" ), call.=FALSE )
    }
  } else maxperm <- 1e100               # No overflow or infinite, just big
  
  # All permutations keep the same first 2 columns.  That way, everything
  #   generated spans the space of all COAs without duplication of ANY
  #   individual orderings.
  if ( ncoa / maxperm > 0.0001 )  {  # any risk of random dups?  Then do all!
    permut <- cbind ( 1, 2,          # first two columns stay same
      matrix ( unlist(combinat::permn(3:ncols)), ncol=ncols-2, byrow=TRUE ) ) 
    permut <- permut[order(runif(ncoa-1))+1,,drop=FALSE]  
                                     # drop first (identity) and
                                     #   shuffle the order too
  } else {                           # just do random ones, toss any dups
    permut <- matrix ( NA, 2*ncoa, ncols )
    for ( i in 1:(ncoa*2) )  {       # times 2 to for any dups (overkill)
      permut[i,] <- c ( 1, 2, 
                        order ( runif ( ncols-2 ) ) + 2 )
    }
    permut <- permut[!duplicated(permut),,drop=FALSE] 
                                     # drop any dups; not likely
    permut <- permut[1:(ncoa-1),,drop=FALSE]  # keep just as many as we need
  }
    
  coalist <- list ( coades )              # start list of multiple COAs
  for ( i in 2:ncoa )  coalist <- c ( coalist, list ( coades[,permut[i-1,]] ) )
  return ( coalist )                    # permute as decided ...
}

#' Generate COA or COAs in final form for use elsewhere
#' 
#' This is the function intended for end-users to generate a COA, or a group
#' of non-overlapping ones, or a stacked COA made from the non-overlapping
#' ones.
#' 
#' If `nvars` is not prime or a power of a prime, this function reduces the 
#' base COA to the right size.  

#' Note that the set of COAs generated by `ncoa>1` is *never* reproducible,
#' unless the caller manages `set.seed` to make it so.  A single COA, however,
#' is automatically reproducible from call to call unless `random` is `TRUE`.
#' 
#' @param nvars The number of items/variables needed in the COA(s).  Cannot
#'   exceed 1000.
#' @param ncoa Number of COAs to generate; see `sep` below.
#' @param sep Ignored when `ncoa=1`.  If `TRUE`, return a list of separate COAs.  
#'   If `FALSE`, return a single "stacked" COA combining the `ncoa` 
#'   separate ones.
#' @param random If `TRUE`, pick a random initial permutation to form the 
#'   COA, meaning the result will NOT reproduce on future runs. 
#' @param nopermute If `TRUE`, the original algorithmic COA's columns are
#'   *not* permuted, resulting in a COA with very bad triplet imbalance.  
#'   This should never be used except for research or play.
#' @inheritDotParams getOA nocheck silent 
#' @returns A COA in 1-based ncol>nrow format as expected by most analysis 
#'   drivers.  If `ncoa=1`, it will have `nvars` rows and at least 
#'   `nvars * (nvars-1)` columns, more if `nvars` is not a prime or power 
#'   of a prime.  If `ncoa>1`, the result will be correspondingly larger if
#'   `sep` is `FALSE`, or a list of such COAs if `sep` is `TRUE`.
#' @importFrom stats runif
#' @export
getCOA <- function ( nvars, ncoa=1, sep=FALSE, random=FALSE, 
                     nopermute=FALSE, ... )  {
  # We call makeCOA, use permuteCOA if we need multiple variants, use dropCOA 
  # to trim the results if needed, and then add 1 and transpose to get into the 
  # layout used by analytic routines.
  
  coades <- makeCOA ( nvars, random=random, ... )        # An initial COA
  ncols <- ncol ( coades )           # can be >nvars if nvars not prime^n[-1]
  
  # A COA straight from the created OA has bad triplet-imbalance, EVEN IF
  #   random=TRUE is passed to the OA creation routine.
  # We fix that by randomly permuting the columns.  But we use a predefined
  #   permutation order to guarantee reproducibility, unless random=TRUE
  #   is given to us.
  if ( !nopermute )  {                   # Skip this next for testing only
    if ( random )  { porder <- order ( runif (ncols) )
    } else porder <- order ( permorder[1:ncols] )
    coades <- coades[,porder]                # The actual shuffling
  }
    
  # Now, if asked for multiple COAs, we further permute columns, always
  #   randomly, making a list of separate COAs.
  if ( ncoa <=1 )  { coaout <- list(coades) # make list for consistency with >1
  } else {
    coaout <- permuteCOA ( coades, ncoa )
  } 
  
  # If we generated for a larger number than wanted, drop the excess items.
  for ( i in 1:ncoa )  {
    if ( ncols != nvars )  coaout[[i]] <- dropCOA ( coaout[[i]], nvars )
                                        # if generated larger, cut it down
                                        # note this is AFTER permuting variants
    coaout[[i]] <- t ( coaout[[i]] + 1L )  # flip & use 1...n coding for final
  }
  
  # And finally, "stack" the separate COAs, if desired.
  if ( sep )  {  return ( coaout )       # return multiple separate ones as list
  } else {
    coastack <- do.call ( cbind, coaout )   # return all as one multiple design
    return ( coastack )
  }
}
