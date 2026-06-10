# Test Drivers 1.R

# ==== Get data to test on, select random vars to test
data_path <- system.file ( "GSSdata.RDS", package="ShapCOA" )
testd <- readRDS ( data_path )
testd <- data.frame ( as.matrix ( testd ) )

size1 <- floor ( runif(1) * 4 ) + 5         # how many variables to test on, 5-8
      # == If range 5:8 is changed, do NOT include any non-primes requiring
      #      drops of more than 1 item from COAs, else max/min won't match
ivars <- sample ( 16, size1 )               # pick variables to test
cat ( "Testing KDRs on variables:", ivars, "\n" )
testdd <- testd[,c(ivars,17)]
scpMatrix <- cpBuild ( testdd )

ordall <- genAllOrders ( size1 )

#     relaimpo exact reference
rel1 <- relaimpo::calc.relimp ( DepVar~., testdd, type="lmg" )$lmg
names(rel1) <- NULL           # Comparees don't keep names

#     my combinations version, exact
comb1 <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, 
                       scpMatrix=scpMatrix, depvar=ncol(scpMatrix) )
test_that("Combos driver matches RelaImpo", {
  expect_equal ( comb1, rel1 )
})

# Next few are duplicated by dtest calls below, but these are the simplest
#   and easiest-to-read tests

ord1 <- SVsByOrders ( ordall, , scpMatrix=scpMatrix, multi=TRUE )     # general driver
test_that("Orders driver matches Combos", {
  expect_equal ( ord1, comb1 )
})

ord2 <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix )   # KDR-specific
test_that("New KDR orders driver matches Combos", {
  expect_equal ( ord2, comb1 )
})

ord3 <- SVsByOrders1 ( ordall, KDRsolve1Combo, scpMatrix=scpMatrix,
                       depvar=nrow(scpMatrix) )  
                                    # 1 combo at a time
test_that("Dumb orders driver matches Combos", {
  expect_equal ( ord3, comb1 )
})


# dtest: run as many drivers as apply, and compare for a given
#        set of parameters
dtest <- function ( depvar=NULL, ndep=NULL, 
                    minitems=NULL, maxitems=NULL, adjusted=NULL, both=NULL )  {
  arglist <- as.list ( match.call() )[-1] 
  arglist$depvar <- NULL; arglist$ndep <- NULL
  callargs <- paste ( names(arglist), arglist,  sep="=", collapse=" " )
  # shamelessly inherits size1, scpMatrix and ordall
  # ... should be used to pass:
  #     minitems, maxitems -- size limited ranges
  #     adjusted, both -- to get adjusted r-squareds
  
  args <- list(nitems=size1,multi=TRUE,silent=TRUE,scpMatrix=scpMatrix)
  if ( !is.null(depvar) )    args <- c(args,list(depvar=rev(depvar)))
#  if ( !is.null(ndep) )      args <- c(args,ndep=ndep)
  if ( !is.null(minitems) )   args <- c(args,minitems=minitems)
  if ( !is.null(maxitems) )   args <- c(args,maxitems=maxitems)
  if ( !is.null(adjusted) )  args <- c(args,adjusted=adjusted)
  if ( !is.null(both) )      args <- c(args,both=both)
  comb1 <- do.call ( SVsByCombos, args )
  #    General orderings driver
  ok <- is.null(adjusted) && is.null(both) 
  if ( ok )  {
    args <- list(orders=ordall, multi=TRUE, scpMatrix=scpMatrix)
    if ( !is.null(minitems) )   args <- c(args,minitems=minitems)
    if ( !is.null(maxitems) )   args <- c(args,maxitems=maxitems)
    if ( !is.null(ndep) )      args <- c(args,ndep=ndep)
    ord1 <- do.call ( SVsByOrders, args )
    test_that(paste("General orders driver vs combo",callargs), {
      expect_equal ( ord1, comb1 )
    })
  }
  #    KDR-specific driver
  args <- list(orders=ordall, scpMatrix=scpMatrix)
  if ( !is.null(minitems) )   args <- c(args,minitems=minitems)
  if ( !is.null(maxitems) )   args <- c(args,maxitems=maxitems)
  if ( !is.null(ndep) )      args <- c(args,ndep=ndep)
  if ( !is.null(adjusted) )  args <- c(args,adjusted=adjusted)
  if ( !is.null(both) )      args <- c(args,both=both)
  ord2 <- do.call ( SVsForKDRsByOrders, args )
  test_that(paste("KDR orders driver vs combo",callargs), {
    expect_equal ( ord2, comb1 )
  })
  
  #    Dumb 1-combo at a time driver
  ok <- length(depvar) == 1 && is.null(adjusted) && is.null(both) 
  if ( ok )  {
    if ( !is.null(maxitems) && maxitems > 0 )  {
           ordall2 <- ordall[1:maxitems,,drop=FALSE]
    } else ordall2 <- ordall
    args <- list(orders=ordall2, vfunc=KDRsolve1Combo, scpMatrix=scpMatrix)
    if ( !is.null(depvar) )    args <- c(args,depvar=depvar)
    if ( !is.null(minitems) )   args <- c(args,minitems=minitems)
    ord3 <- do.call ( SVsByOrders1, args )
    test_that(paste("Dumb orders driver vs combo",callargs), {
      expect_equal ( ord3, comb1 )
    })
  }
  invisible(comb1)             # return one for further reconciliations
}

# Simple, plus adjusteds
dtest ( depvar=ncol(scpMatrix), ndep=1 )
dtest ( depvar=ncol(scpMatrix), ndep=1, adjusted=TRUE )
dtest ( depvar=ncol(scpMatrix), ndep=1, adjusted=TRUE, both=TRUE )
dtest ( depvar=ncol(scpMatrix), ndep=1, both=TRUE )
# Two dependents
dtest ( depvar=ncol(scpMatrix)-(0:1), ndep=2 )
dtest ( depvar=ncol(scpMatrix)-(0:1), ndep=2, adjusted=TRUE )
dtest ( depvar=ncol(scpMatrix)-(0:1), ndep=2, both=TRUE )
# Maximum sizes
dtest ( depvar=ncol(scpMatrix), ndep=1, maxitems=size1-3 )
dtest ( depvar=ncol(scpMatrix), ndep=1, maxitems=size1-3, adjusted=TRUE )
dtest ( depvar=ncol(scpMatrix)-(0:1), ndep=2, maxitems=size1-3, both=TRUE )
# Minimum sizes
dtest ( depvar=ncol(scpMatrix), ndep=1, minitems=size1-3 )
dtest ( depvar=ncol(scpMatrix), ndep=1, minitems=size1-3, adjusted=TRUE )
dtest ( depvar=ncol(scpMatrix)-(0:1), ndep=2, minitems=size1-3, both=TRUE )

# Now collect mins and maxs (returned for the combo version), and reconcile.
for ( ipass in 1:2 )  {
  ndep <- ipass
  depvars <- ncol(scpMatrix) - ( if(ipass==2) 0:1 else 0 )
  both <- ipass == 2
  comb0 <- dtest ( depvar=depvars, ndep=ndep, both=both )    # Total to save
  for ( sizecut in 1:(size1-1) )  {
    combmn <- dtest ( depvar=depvars, ndep=ndep, minitems=sizecut, both=both )
    combmx <- dtest ( depvar=depvars, ndep=ndep, maxitems=sizecut, both=both )
    test_that ( "", {
      expect_equal ( ( sizecut*combmx + (size1-sizecut) * combmn ) / size1,
                     comb0 )
    })
  }
}



# ==== Repeat with adjusted r-squareds
comb1a <- SVsByCombos ( KDRsolveCombos, size1, multi=TRUE, silent=TRUE,
                       scpMatrix=scpMatrix, depvar=ncol(scpMatrix),
                       adjusted=TRUE )

ord2a <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix, adjusted=TRUE ) 
test_that("New KDR orders driver ADJ rsq matches Combos", {
  expect_equal ( ord2a, comb1a )
})

ord2b <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix, 
                              adjusted=TRUE, both=TRUE )
comb2b <- cbind(comb1,comb1a)
attr(comb2b,"dimnames") <- NULL
test_that("Orders and combos match raw/adj", {
  expect_equal ( ord2b, comb2b )
})


# ==== Test maxitems and minitems, for combos and orders
#      Be sure corresponding min/max weight together to overall.
#      Be sure combo versions match orders versions.
for ( ipass in 1:3 )  {
  adjusted <- ipass == 2 
  both <- ipass == 3 
  if ( ipass == 3 )  adjusted <- NULL
#  cat ( "\nTESTING PASS", ipass, "  adjusted", adjusted, "  both", both, "\n" )
  comb1 <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, 
                         scpMatrix=scpMatrix, depvar=ncol(scpMatrix),
                         adjusted=adjusted, both=both)
  testsizes <- sort ( c ( 1, size1-1, sample ( 2:(size1-2), 2 ) ) )
  for ( nitems in testsizes )  {
    combmn <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, minitems=nitems, 
                            scpMatrix=scpMatrix, depvar=ncol(scpMatrix),
                            adjusted=adjusted, both=both)
    combmx <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, maxitems=nitems,
                            scpMatrix=scpMatrix, depvar=ncol(scpMatrix),
                            adjusted=adjusted, both=both)
    combined <- ( combmx * nitems + combmn * (size1-nitems) ) / size1
    test_that ( paste0 ( "Combo min/max match for nitems=", nitems, 
                         "ipass=", ipass ), {
      expect_equal ( combined, comb1 )
    })
  
    ordmn <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix, minitems=nitems, 
                                  adjusted=adjusted, both=both)
    ordmx <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix, maxitems=nitems, 
                                  adjusted=adjusted, both=both)
    test_that ( paste0 ( "Combo/orders mins match for nitems=", nitems ), {
      expect_equal ( ordmn, combmn )
    })
    test_that ( paste0 ( "Combo/orders maxs match for nitems=", nitems ), {
      expect_equal ( ordmx, combmx )
    })
  }
}


# ==== Test SVsByOrders1
#      Note that value function must accept 2-dim array orderings
macall <- SVsByOrders1 ( ordall, KDRsolveCombos1_cpp, scpMatrix=scpMatrix,
                         depvar=ncol(scpMatrix) )
test_that("SVsByOrders1 on all orders matches RelaImpo", {
  expect_equal ( rel1, macall )
})

# ==== Check both SVsByOrders1 and SVsByCOmbos for partials growing 
#      as maxitems goes down, and for matching each other
macallb <- macall                       # prior size results
for ( keep in (size1-1):1 )  {
  combp <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE,
                         maxitems=keep, minitems=0, 
                         scpMatrix=scpMatrix, depvar=ncol(scpMatrix) )
  ordallp <- ordall[1:keep,,drop=FALSE]
  macallp <- SVsByOrders1 ( ordallp, KDRsolveCombos1_cpp, scpMatrix=scpMatrix,
                          depvar=ncol(scpMatrix) )
  test_that("SVsByOrders1 matches SVsByCombos on partials", {
    expect_equal ( macallp, combp )
  })
  test_that("SVsByOrders1 partials exceed fulls", {
    expect_gt ( sum(macallp), sum(macallb) )
  })
  macallb <- macallp
}
# Final partials of size 1 should be just univariate r2's
unicor2 <- cor ( testdd ) [-(size1+1),size1+1] ^ 2
names(unicor2) <- NULL
test_that("SVsByOrders1 partial 1s match r2s", {
  expect_equal ( macallb, unicor2 )
})

# ==== Check both SVsByOrders1 and SVsByCOmbos for partials shrinking 
#      as minitems goes up, and for matching each other
macallb <- macall                       # prior size results
for ( minitems in 0:(size1-1) )  {
  combp <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, minitems=minitems, 
                         scpMatrix=scpMatrix, depvar=ncol(scpMatrix) )
  macallp <- SVsByOrders1 ( ordall, KDRsolveCombos1_cpp, minitems=minitems,
                            scpMatrix=scpMatrix, depvar=ncol(scpMatrix) )
  test_that("SVsByOrders1 matches SVsByCombos on minitems", {
    expect_equal ( combp, macallp ) 
  })
  if ( minitems > 0 )  {
    test_that("SVsByOrders1 partials higher min exceeds lower", {
      expect_lt ( sum(macallp), sum(macallb) )
    })
    macallb <- macallp
  }
}

# ==== Check orderings approaches for maxitems matching combos
maxitems <- floor ( size1 / 2 )
comb1mx <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, maxitems=maxitems, 
                         scpMatrix=scpMatrix, depvar=ncol(scpMatrix),
                         adjusted=TRUE, both=TRUE )
ord2mx <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix, 
                               adjusted=TRUE, both=TRUE,
                               maxitems=maxitems )
test_that("Combos and orders on maxitems", {
  expect_equal ( comb1mx, ord2mx )
})

minitems <- 3 
comb1mn <- SVsByCombos ( , size1, multi=TRUE, silent=TRUE, minitems=minitems, 
                         scpMatrix=scpMatrix, depvar=ncol(scpMatrix),
                         adjusted=TRUE, both=TRUE )
ord2mn <- SVsForKDRsByOrders ( ordall, scpMatrix=scpMatrix, 
                               adjusted=TRUE, both=TRUE,
                               minitems=minitems )
test_that("Combos and orders on minitems", {
  expect_equal ( comb1mn, ord2mn )
})


