# Use usethis::use_test() to create these files

# ==== Get data to test on, select random vars to test
data_path <- system.file ( "GSSdata.RDS", package="ShapCOA" )
cat ( "data path is ", data_path, "\n" )
testd <- readRDS ( data_path )
# cat ( "WD is ", getwd(), "\n" )
# testd <- readRDS ( "inst/GSSdata.RDS" )     # 16-variable KDR data, plus depvar
# testd <- readRDS ( "../../inst/GSSdata.RDS" )     # 16-variable KDR data, plus depvar
            # Note that tests are run from tests/testthat
testd <- data.frame ( as.matrix ( testd ) )

size1 <- floor ( runif(1) * 7 ) + 7         # how many variables to test on
ivars <- sample ( 16, size1 )               # pick variables to test
cat ( "Testing KDRs on variables:", ivars, "\n" )

# ==== Setup weighted stuff for testing
testdd <- testd[,c(ivars,17)]
wtest <- sample ( nrow(testdd), floor(nrow(testdd)/2 ) )  # half the rows

weights <- rep ( 0, nrow(testdd) )
weights[wtest] <- 4        # will test scale invariance et al
testunw <- testdd[wtest,]  # will test weighting for subsetting
testww <- testdd
testww$weight <- weights

# ==== Now basic cpBuild tests
CPwithC <- cpBuild ( testdd, precenter=FALSE )
CPwithout <- cpBuild ( testdd )
test_that("Pre-center and constant-out give same CP", {
  expect_equal ( CPwithout, cpConstantOut ( CPwithC ) )
})
# CPsweptout <- cpConstantOut ( CPwithC )

# Test weight scale invariance
CPwtd1 <- cpBuild ( testdd, weights=rep(1,nrow(testdd)) )
test_that("Weights of 1 make no difference",  {
  expect_equal ( CPwithout, CPwtd1 )
  
})
CPwtd2 <- cpBuild ( testdd, weights=rep(2,nrow(testdd)) )
test_that("Weights of 2 make no difference",  {
  expect_equal ( CPwithout, CPwtd1 )
})
  
# Weights as subsetters
CPunwtd <- cpBuild ( testunw )
CPwtd <- cpBuild ( testww, weights="weight" )
test_that("Weights as subsetters works", {
  expect_equal ( CPunwtd, CPwtd/4 )
})

# Finding weight and dependent columns 
nc <- ncol(testww)
testwwback <- testww[,c(1,nc,2,nc-1,3:(nc-2))]      # reverse order of columns
CPwtdb <- cpBuild ( testww, weights="weight", depvar="DepVar" )
test_that("Specifying depvar and reordering still works", {
  expect_equal ( CPwtd, CPwtdb )
} )

# Multiple dependent variables
deps <- c(2,3,which(colnames(testdd)=="DepVar"))
ndep <- length(deps)
nondeps <- (1:ncol(testdd)) [-deps]
depnames <- colnames(testdd)[deps]
CPmultDV <- cpBuild ( testdd, depvar=depnames )
for ( i in seq_along(deps) )  {
  cols <- c ( 1:length(nondeps), length(nondeps)+i )
  CPmpart <- CPmultDV[cols,cols]         # our multi-DV part for just one DV
  ocols <- c ( nondeps, deps[i] )        # pull one DV from original
  ocols <- colnames(CPmpart)
  CP1part <- CPwithout[ocols,ocols]      # 
  test_that("Multiple DVs work right", {
    expect_equal ( CPmpart, CP1part )
  } )
}

# ==== Depvars as either names or numbers
test_that("deps= numbers same as deps=names", {
  expect_equal( cpBuild ( testdd, depvar=c(2,4) ),
                cpBuild ( testdd, depvar=colnames(testdd)[c(2,4)] ) )
})
test_that("deps= numbers same as default", {
  expect_equal( cpBuild ( testdd ),
                cpBuild ( testdd, depvar=ncol(testdd) ) )
})
 
# ==== Test checkweights
test_that("checkweights picks up on wrong lengths", {
  expect_false ( checkWeights ( runif(23), 24, "Length Test" ) )
} )
test_that("checkweights picks up on negativess", {
  expect_false ( checkWeights ( c(runif(22),-1,1,-1,1), 26, "Non-neg Test" ) )
} )


# ========== Now test KDR solutions
nvtot <- ncol(testd) - 1
cvars <- sample ( nvtot, floor(2/3*nvtot) )  # A combination
testcc <- testd[,c(cvars,ncol(testd))]
lmcoef <- lm ( DepVar ~ ., testcc )

test_that("KDRsolve1combo matches lm r-squared", {
  expect_equal (  summary ( lm ( DepVar ~ ., testcc ) )$r.squared,
                  KDRsolve1Combo ( cvars, cpBuild(testd), depvar="DepVar" ) )
} )
   
# Multiple combo solves, for 1/more dependents
tsize <- floor ( nvtot * ( runif(1) * 0.5 + 0.3 ) )
cat ( "Testing multiples of size", tsize, "in KDRsolveComboR\n" )
combos <- matrix ( NA, tsize, 8 )
for ( i in 1:ncol(combos) )  combos[,i] <- sample ( 1:nvtot, tsize )
kdrsq <- KDRsolveComboR ( combos, cpBuild(testd), depvar="DepVar" )
# kdrsq <- KDRsolveComboR ( combos, cpBuild(testd), depvar="DepVar" )[,1]
lmrsq <- rep(NA,ncol(combos) )
for ( i in 1:ncol(combos) )  {
  testcc <- testd[,c(combos[,i],ncol(testd))]
  lmrsq[i] <- summary ( lm ( DepVar ~ ., testcc ) )$r.squared
}
test_that("KDRsolveComboR matches lm, 1 dependent variable", {
  expect_equal ( kdrsq, lmrsq )
})

# Now multiple DVs
testd2 <- testd
testd2$DepVar2 <- testd$DepVar * runif(nrow(testd))
kdrsq <- KDRsolveComboR ( combos, cpBuild(testd2), depvar=c("DepVar","DepVar2") )
lmrsq <- matrix(NA,ncol(combos),2 )
for ( i in 1:ncol(combos) )  {
  testcc <- testd2[,c(combos[,i],ncol(testd))]
  lmrsq[i,1] <- summary ( lm ( DepVar ~ ., testcc ) )$r.squared
  testcc <- testd2[,c(combos[,i],ncol(testd2))]
  lmrsq[i,2] <- summary ( lm ( DepVar2 ~ ., testcc ) )$r.squared
}
test_that("KDRsolveComboR matches lm, 2 dependent variables", {
  expect_equal ( kdrsq, lmrsq )
})


# ==== Now the "real" KDR/combo solver
rsqstr <- KDRsolveCombos ( combos, cpBuild(testd),  depvar="DepVar" )
rsqadj <- KDRsolveCombos ( combos, cpBuild(testd), depvar="DepVar", 
                           adjusted=TRUE )
rsqbot <- KDRsolveCombos ( combos, cpBuild(testd), depvar="DepVar", 
                           both=TRUE )
test_that("KDRsolveCombos straight/adjusted match both", {
  expect_equal ( cbind(rsqstr,rsqadj), rsqbot )
})
lmrsq <- matrix(NA,ncol(combos),2 )
for ( i in 1:ncol(combos) )  {
  testcc <- testd[,c(combos[,i],ncol(testd))]
  sumtemp <-  summary ( lm ( DepVar ~ ., testcc ) )
  lmrsq[i,1] <- sumtemp$r.squared
  lmrsq[i,2] <- sumtemp$adj.r.squared
}
test_that("KDRsolveCombos straight/adjusted match lm", {
  expect_equal ( rsqbot, lmrsq )
})

# Now repeat with 2 dependents
deps <- c("DepVar","DepVar2")
rsqstr <- KDRsolveCombos ( combos, cpBuild(testd2,depvar=deps), depvar=deps )
rsqadj <- KDRsolveCombos ( combos, cpBuild(testd2,depvar=deps), depvar=deps, 
                           adjusted=TRUE )
rsqbot <- KDRsolveCombos ( combos, cpBuild(testd2,depvar=deps), depvar=deps, 
                           both=TRUE )
test_that("KDRsolveCombos straight/adj match both for 2 DVs", {
  expect_equal ( cbind(rsqstr,rsqadj), rsqbot )
})

lmrsq <- matrix(NA,ncol(combos),4 )
for ( i in 1:ncol(combos) )  {
  testcc <- testd[,c(combos[,i],ncol(testd))]
  sumtemp <-  summary ( lm ( DepVar ~ ., testcc ) )
  lmrsq[i,1] <- sumtemp$r.squared
  lmrsq[i,3] <- sumtemp$adj.r.squared
  testcc <- testd2[,c(combos[,i],ncol(testd2))]
  sumtemp <-  summary ( lm ( DepVar2 ~ ., testcc ) )
  lmrsq[i,2] <- sumtemp$r.squared
  lmrsq[i,4] <- sumtemp$adj.r.squared
}
test_that("KDRsolveCombos straight/adj match lm for 2 DVs", {
  expect_equal ( rsqbot, lmrsq )
})
