# ==== primeSize
cat ( "\n" )
primes <- primeSize ( 175, all=TRUE )
test_that("primeSize looks OK up to 175", {
  expect_equal ( c (length(primes), sum(primes), sum(log(primes)) ),
                 c ( 54, 3945, 207.553335 ) )
})
 
psmall <- primes[primes>10&primes<50]         # Smallish ones for faster tests
testsize <- psmall[order(runif(length(psmall)))[1]]
cat ( "Testing COA of size", testsize, "\n" )

testCOA <- getCOA ( testsize )
testOA <- getOA ( testsize )

nonprimes <- min(psmall):max(psmall) 
nonprimes <- nonprimes[is.na(match(nonprimes,psmall))]
psize <- length(nonprimes)
for ( i in seq_along(nonprimes) )  {
  psize[i] <- primeSize ( nonprimes[i] )    
}
test_that("primeSize grows non-primes sizes to primes", {
  expect_in ( psize, primes )
})


# ==== assessCOA, getCOA
test_that("Basic COA checks out", {
  expect_true ( assessCOA ( testCOA ) )
})

test_that("Non-prime COA is right size", {
  expect_equal ( dim (getCOA(14)), c(14,240) )
})

test_that("COAs are reproducible", {
  expect_identical ( getCOA(testsize), getCOA(testsize) )
})

test_that("COAs are randomoized when asked", {
  expect_false ( identical ( getCOA(testsize, random=TRUE ),
                             getCOA(testsize, random=TRUE ) ) )
})

# ==== permuteCOA
test_that("Permuted COA still checks out", {
  expect_true ( ckCOA ( permuteCOA ( designCode ( testCOA ) )[[2]] ) )
})

nperm <- floor ( runif(1)*5 ) + 3
cat ( "Testing", nperm, "permutations.\n" )

perm3 <- permuteCOA ( designCode ( testCOA ), nperm )
test_that("permuteCOA delivers correct length", {
  expect_equal ( nperm, length(perm3) )
})

test_that("permuteCOA preserves original COA", {
  expect_identical ( perm3[[1]], designCode ( testCOA ) )
})

coastack <- getCOA ( testsize, nperm * 2 )
test_that("Non-duplication in permuted COAs", {
  expect_false ( any ( duplicated ( coastack ) ) )
})

# ==== dropCOA
test_that("Drop-1 COA still checks out", {
  expect_true ( ckCOA ( dropCOA ( designCode ( getCOA (16) ), 15 ) ) )
})
test_that("Drop-2 COA does not check out", {
  expect_false ( ckCOA ( dropCOA ( designCode ( getCOA (16) ), 14 ) ) )
})

# ==== makeCOA is tested indirectly, through getCOA()

# ==== userCode and designCode
test_that("userCode inverts designCode for a COA", {
  expect_equal ( testCOA, userCode ( designCode ( testCOA ) ) )
})
dtest <- designCode ( testCOA )
test_that("designCode inverts userCode for a COA", {
  expect_equal ( dtest, designCode ( userCode ( dtest ) ) )
})

# ==== assessOA, getOA
test_that("Basic OA checks out", {
  expect_true ( assessOA ( designCode ( testOA ) ) )
})

# ==== userCode, designCode on 0-1 transposed designs
otest <- userCode ( testOA )
test_that("userCode inverts designCode for an OA", {
  expect_equal ( otest, userCode ( designCode ( userCode ( otest ) ) ) )
})
test_that("designCode inverts userCode for an OA", {
  expect_equal ( testOA, designCode ( userCode ( testOA ) ) )
})

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

