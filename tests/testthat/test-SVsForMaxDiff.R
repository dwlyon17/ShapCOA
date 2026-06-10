# Test SVsOnMaxDiffs2.R

# ==== Get 0/1 TURF data to test on, select random vars to test
data_path <- system.file ( "TURFtest.RDS", package="ShapCOA" )
cat ( "data path is ", data_path, "\n" )
testd <- readRDS ( data_path )

size1 <- floor ( runif(1) * 5 ) + 4    # how many variables to test on, 4-8
      # == If range 5:8 is changed, do NOT include any non-primes requiring
      #      drops of more than 1 item from COAs, else max/min won't match
ivars <- sample ( ncol(testd), size1 )            # pick variables to test
cat ( "Testing 0/1 TURF on variables:", ivars, "\n" )
testd1 <- testd[,ivars]

ordall <- genAllOrders ( size1 )

# ==== Test all-orders small problem matches 
# Try to match standard combos and orders on simple-minded value function
comb1 <- SVsByCombos ( turf1, size1, multi=FALSE, silent=TRUE,
                       tdata=testd1, depth=1 )
ord1 <- SVsByOrders1 ( ordall, turf1, tdata=testd1, depth=1 )
# 1 combo at a time
test_that("Pre-check: dumb orders driver matches Combos", {
  expect_equal ( comb1, ord1 )
})
mxd1 <- SVsForMaxDiff ( ordall, testd1 )
test_that("MD driver for TURF: matches Combos", {
  expect_equal ( comb1, mxd1 )
})
# See if depths behave correctly
for ( depth in seq(1,size1-1,2) )  {
  combd <- SVsByCombos ( turf1, size1, multi=FALSE, silent=TRUE,
                         tdata=testd1, depth=depth )
  mxdd <- SVsForMaxDiff ( ordall, testd1, threshold=depth )
  test_that(paste ( "Depth", depth, "MD driver for TURF: matches Combos" ), {
    expect_equal ( combd, mxdd )
  })
}
test_that ( "Depth==# items gives same reach for all items", {
  mxdd <- SVsForMaxDiff ( ordall, testd1, threshold=size1 )
  expect_equal ( length(unique(mxdd)), 1 )
})

# ==== Draw a bigger test set and check max/min item congruence 
pp <- primeSize ( 30 )                   # all primes up to 30
equipos <- unique ( c ( pp, pp-1 ) )     # all n's w/ equipositioning
equipos <- equipos[equipos>=10]
size2 <- sample ( equipos, 1 )           # pick an equi-posnd size to test
ivars <- sample ( ncol(testd), size2 )            # pick variables to test
cat ( "Testing 0/1 TURF on variables:", ivars, "\n" )
testd2 <- testd[,ivars]
coa <- getCOA ( size2 )

# See if max and min reconcile exactly.  This only works for equi-positioned
#   "true" COAs!
for ( depth in c(2,4) )  {
  testsizes <- sort ( unique ( c(depth, size2-depth, sample ( size2-1, 3 )) ))
  mxdda <- SVsForMaxDiff ( coa, testd2, threshold=depth )
  for ( sizecut in testsizes )  {
    mxddn <- SVsForMaxDiff ( coa, testd2, threshold=depth, minitems=sizecut )
    mxddx <- SVsForMaxDiff ( coa, testd2, threshold=depth, maxitems=sizecut )
    test_that ( paste("TURF","depth=",depth,"sizecut=",sizecut), {
      expect_equal ( ( sizecut*mxddx + (size2-sizecut) * mxddn ) / size2,
                     mxdda )
    })
  }
}

# Be sure different COAs relate well.
ncoa <- 4 
mxddc <- matrix ( NA, size2, ncoa )
for ( i in 1:ncoa )  {
  mxddc[,i] <- SVsForMaxDiff ( , testd2, threshold=depth, 
                               random=TRUE, silent=TRUE, nocheck=TRUE )
}
test_that ( "High correlations between COAs", {
  expect_gt ( min(cor(mxddc)), 0.997 )
})


# ========================== Now get some real MaxDiff data.  It's anchored.
#     and is not pre-zeroed on the anchor!
data_path <- system.file ( "MDtester.RDS", package="ShapCOA" )
testd <- readRDS ( data_path )

cat ( "Testing anchored MD\n" )
mdasa <- list()
for ( ipass in 1:6 )  {
  threshold <- list ( TRUE, FALSE, 0.5, 0.6, 0.0001, 0.9999 )[[ipass]]
  # Have callee handle the anchor, and implicit xform
  mdasa[[ipass]] <- 
           SVsForMaxDiff ( , testd, anchor=ncol(testd), threshold=threshold,
                           silent=TRUE,  nocheck=TRUE )
  # Anchor it ourselves and get same answer ...
  if ( ipass %in% 2:3 )  {
    testa <- ( testd - testd[,ncol(testd)] ) [,-ncol(testd)] 
    mdasb <- SVsForMaxDiff ( , testa, xform=2, anchor=TRUE, threshold=threshold,
                             silent=TRUE,  nocheck=TRUE )
    test_that ( "We anchor/funct anchors come out same", {
      expect_equal ( mdasa[[ipass]], mdasb )
    })
  }
  # tasksize should not matter
  if ( ipass %in% c(2,4) )  {
    mdasc <- SVsForMaxDiff ( , testa, xform=2, anchor=TRUE, threshold=threshold,
                             silent=TRUE,  nocheck=TRUE, tasksize=18 )
    test_that ( "Tasksize does not affect anchored analysis", {
      expect_equal ( mdasa[[ipass]], mdasc )
    })
  }
}
test_that ( "Anchor threshold higher-> results lower", {
  expect_gt ( sum(mdasa[[3]]), sum(mdasa[[4]]) ) 
})
test_that ( "Anchor threshold defaults to 0.5", {
  expect_equal ( mdasa[[1]], mdasa[[3]] )
})

# ==== Now try MNL analyses, using the anchored data (with anchor removed)
cat ( "Testing MNL Maxdiff\n" )
testdd <- testd[,-ncol(testd)]     # get rid of the anchor
mdmsa <- list()
for ( ipass in 1:7 )  {
  threshold <- list ( TRUE, FALSE, TRUE, 0.5, 0.6, 0.0001, 0.9999 )[[ipass]]
  if ( ipass %in% 3:5 )  { maxitems <- ceiling ( size2/2 )
  } else maxitems <- NULL         # For MNL, threshold doesn't matter w/o max
  # Basic MNL analysis
  mdmsa[[ipass]] <- 
           SVsForMaxDiff ( , testdd, xform=0, threshold=threshold, 
                           maxitems=maxitems, silent=TRUE,  nocheck=TRUE ) 
  
  # Perturb utilities and get same answer ...
  if ( ipass %in% 1:2 )  {
    testp <- testdd + runif(nrow(testdd)) 
    mdmsb <- SVsForMaxDiff ( , testp, xform=0, threshold=threshold,
                             maxitems=maxitems, silent=TRUE,  nocheck=TRUE )
    test_that ( "MNL not affected by centering ", {
      expect_equal ( mdmsa[[ipass]], mdmsb )
    })
    # tasksize should not matter
  }
  if ( ipass %in% 1:2 )  {
    mdmsc <- SVsForMaxDiff ( , testp, xform=0, threshold=threshold,
                             maxitems=maxitems, silent=TRUE,  nocheck=TRUE,
                             tasksize=18 )
    test_that ( "Tasksize does not affect MNL analysis", {
      expect_equal ( mdmsb, mdmsc )
    })
  }
}
test_that ( "MNL threshold defaults to 0.5", {
  expect_equal ( mdmsa[[3]], mdmsa[[4]] )
})
# For MNL, we always get to 1.00 eventually ...
ttm <- sapply ( mdmsa, cbind, simplify=TRUE )
test_that ( "MNL threshold defaults to 0.5", {
  expect_equal ( colSums(ttm)[-(3:5)], rep(1,ncol(ttm))[-(3:5)] )
})
test_that ( "MNL threshold higher-> results lower", {
  expect_gt ( sum(mdasa[[5]]), sum(mdasa[[4]]) ) 
})

# Now treat it as regular non-anchored, default scoring, Maxdiff
cat ( "Testing default MaxDiff\n" )
mddsa <- list()
for ( ipass in 1:7 )  {
  threshold <- list ( TRUE, FALSE, 0.5, 0.9, 0.0001, 0.9999, 0.5 )[[ipass]]
  tasksize <-  c ( rep(4,6),   5 )[ipass]     
  mddsa[[ipass]] <- 
    SVsForMaxDiff ( , testdd, xform=1, threshold=threshold, tasksize=tasksize,
                    silent=TRUE,  nocheck=TRUE )
  # Perturb utilities and get same answer ...
  if ( ipass %in% 1:2 )  {
    testp <- testdd + runif(nrow(testdd)) 
    mddsb <- SVsForMaxDiff ( , testp, xform=1, threshold=threshold,
                             tasksize=tasksize,
                             silent=TRUE,  nocheck=TRUE )
    test_that ( "Default not affected by pre-centering ", {
      expect_equal ( mddsa[[ipass]], mddsb )
    })
  }
}
test_that ( "Default threshold higher-> results lower", {
  expect_gt ( sum(mddsa[[3]]), sum(mddsa[[4]]) ) 
})
test_that ( "Default threshold defaults to 0.9", {
  expect_equal ( mddsa[[1]], mddsa[[4]] )
})
test_that ( "Default No thresh never exceeds 1",  {
  expect_lt ( sum(mddsa[[2]]), 1.0 )
})
test_that ( "Default tasksize matters",  {
  expect_false ( isTRUE ( all.equal ( mddsa[[3]], mddsa[[7]] ) ) )
})






