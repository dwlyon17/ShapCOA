# Test file for UtilsN.R

# ======== Basic min/max relationships first
test_that("Legit min/max work OK", {
  expect_equal( ckMinMax ( 1, 7,, 8 ), c(1,7) )
})

test_that("Defaulted min/max work OK", {
  expect_equal( ckMinMax ( ,,, 8 ), c(0,8) )
})

test_that("Defaulted min/max 0 works OK", {
  expect_equal( ckMinMax ( , 0, , 8 ), c(0,8) )
})

test_that("max exceeds total", {
  expect_error ( ckMinMax ( 0, 9, , 8 ), "exceeds total", 
                 inherit=FALSE )
})

test_that("min exceeds total", {
  expect_error ( ckMinMax ( 9, 8, , 8 ), "must be strictly less", 
                 inherit=FALSE )
})

test_that("min equal to total", {
  expect_error ( ckMinMax ( 8, 8, , 8 ), "must be strictly less", 
                 inherit=FALSE )
})

test_that("min exceeds max", {
  expect_error ( ckMinMax ( 7, 2, , 8 ), "exceeds maxitems", 
                 inherit=FALSE )
})

test_that("min is negative", {
  expect_error ( ckMinMax ( -1, 8, , 8 ), "negative or non-integer", 
                 inherit=FALSE )
})

test_that("min/max fractional", {
  expect_error ( ckMinMax ( 1, 7.1, , 8 ), "negative or non-integer", 
                 inherit=FALSE )
})

test_that("max exceeds total", {
  expect_error ( ckMinMax ( 0, 9, , 8 ), "exceeds total", 
                 inherit=FALSE )
})

# ======== Now deducing things from a data/orderings/SCP matrix

tscp <- matrix ( 1:100, 10, 10 )        # it's square; all that matters here
tcoa <- getCOA ( 12 )
tdata <- matrix ( runif ( 900 ), 100, 9 )  # pretend data matrix
tbad <- rep(1,32)                       # not a matrix

test_that("max exceeds total", {
  expect_error ( ckMinMax ( , , tbad ), "cannot recognize", 
                 inherit=FALSE )
})

test_that("Defaults from COA work", {
  expect_equal( ckMinMax ( , , tcoa ), c(0,nrow(tcoa) ) )
})

test_that("Defaults from size-limited COA work", {
  expect_message( minmax <- ckMinMax ( , , tcoa[1:5,] ), "implied by" )
  expect_equal ( minmax, c(0,5) )
})

test_that("orders match neither total nor length", {
  expect_error ( ckMinMax ( 0, 4, tcoa[1:5,], 8 ), "is neither", 
                 inherit=FALSE )
})

test_that("Defaults from data work", {
  expect_equal( ckMinMax ( , , tdata ), c(0,ncol(tdata) ) )
})

test_that("Defaults from SCP work", {
  expect_equal( ckMinMax ( , , tscp ), c(0,nrow(tscp)-1)  )
})

test_that("Defaults from size-limited COA work", {
  expect_equal( ckMinMax ( , , tcoa[1:5,] ), c(0,5) )
})
 
# ==== Test checkweights
test_that("checkweights picks up on wrong lengths", {
  expect_error ( checkWeights ( runif(23), 24 ), "rows but only",
                 inherit=FALSE )
} )
test_that("checkweights picks up on negativess", {
  expect_error ( checkWeights ( c(runif(22),-1,1,-1,1), 26 ),
                 "Negative weights are not allowed", 
                 inherit=FALSE )
} )




