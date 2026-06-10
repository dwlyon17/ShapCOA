
# Testing for Orders 1.R

ckOrders <- function ( label, orders, ckrows=TRUE )  {
  test_that(paste0(label," in range"), {
    expect_true ( min(orders) == 1 && max(orders) == size1 )
  })
  test_that(paste0(label," all items"), {  # all cols contain all items
    expect_equal ( 1, length(unique(colSums(orders))) )
  })
  if ( ckrows )  {
    test_that(paste0(label," all balanced"), {
      expect_equal ( 1, length(unique(rowSums(orders))) )
    })
  }
}

size1 <- floor ( runif(1) * 4 ) + 5

# ==== genAllOrders
genall <- genAllOrders ( size1 )
ckOrders ( "genAllOrders", genall, TRUE )

# ==== genyclic
ncycles <- floor ( runif(1) * 10 + 5 )
norders <- floor ( runif(1) * 10 + 5 ) * size1
norders3 <- floor ( runif(1) * 100 ) + 100
gencyc1 <- genCyclic ( size1, norders=norders )
gencyc1x <- genCyclic ( size1, norders=norders )
gencyc1y <- genCyclic ( size1, norders=norders )
gencyc2 <- genCyclic ( size1, ncycles=ncycles )
gencyc3 <- genCyclic ( size1, norders=norders3 )
test_that("genCyclic right size by norders", {
  expect_equal ( norders, ncol(gencyc1) )
})
test_that("genCyclic right size by ncycles", {
  expect_equal ( ncycles * size1, ncol(gencyc2) )
})
test_that("genCyclic norders rounded up right", {
  expect_true ( ncol(gencyc3) == norders3 ||
                ncol(gencyc3) %/% size1 == norders3 %/% size1 + 1 )
})
test_that("genCyclic reproduces", {
  expect_identical ( genCyclic ( size1 ), genCyclic ( size1 ) )
})
test_that("genCyclic randomizes when asked", {
  expect_false ( identical ( genCyclic ( 29, random=TRUE ),
                             genCyclic ( 29, random=TRUE ) ) )
})                   # chance of accidental duplication for 29 is nil!
test_that("genCyclic randomizes automatically", {
  expect_false ( identical ( gencyc1x, gencyc1y ) )
})
ckOrders ( "genCyclic 1", gencyc1, TRUE )
ckOrders ( "genCyclic 2", gencyc2, TRUE )
ckOrders ( "genCyclic 3", gencyc3, TRUE )


# ==== genWilliams
ncycles <- floor ( runif(1) * 10 + 5 )
norders <- floor ( runif(1) * 10 + 5 ) * 2 * size1
norders3 <- floor ( runif(1) * 100 ) + 100
genwill1 <- genWilliams ( size1, norders=norders )
genwill1x <- genWilliams ( size1, norders=norders )
genwill1y <- genWilliams ( size1, norders=norders )
genwill2 <- genWilliams ( size1, ncycles=ncycles )
genwill3 <- genWilliams ( size1, norders=norders3 )
test_that("genWilliams right size by norders", {
  expect_equal ( norders, ncol(genwill1) )
})
test_that("genWilliams right size by ncycles", {
  expect_equal ( ncycles * size1 *( 2), ncol(genwill2) )
})
test_that("genWilliams norders rounded up right", {
  expect_true ( ncol(genwill3) == norders3 ||
                ncol(genwill3) %/% (2*size1) == norders3 %/% (2*size1) + 1 )
})
test_that("genWilliams randomizes automatically", {
  expect_false ( identical ( genwill1x, genwill1y ) )
})
test_that("genWilliams reproduces", {
  expect_identical ( genWilliams ( size1, ncycles=1 ), 
                     genWilliams ( size1, ncycles=1 ) )
})
test_that("genWilliams randomizes when asked", {
  expect_false ( identical ( genWilliams ( 29, ncycles=1, random=TRUE ),
                             genWilliams ( 29, ncycles=1, random=TRUE ) ) )
})                   # chance of accidental duplication for 29 is nil!
ckOrders ( "genWilliams 1", genwill1, TRUE )
ckOrders ( "genWilliams 2", genwill2, TRUE )
ckOrders ( "genWilliams 3", genwill3, TRUE )
test_that("genWilliams passes assessCOA", {
  expect_true ( assessCOA ( designCode ( genwill1 ) ) )
})

# ==== genRandOrders
norders <- floor ( runif(1) * 1000 + 400 ) 
genrand <- genRandOrders ( size1, norders )
test_that("genRandOrders right size by norders", {
  expect_equal ( norders, ncol(genrand) )
})
ckOrders ( "genRandOrders", genrand, FALSE )

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
