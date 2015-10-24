library(Rinbix)
context("Utilities")

test_that("fisherRtoZ check range from -1 to 1", {
  testValues <- seq(from=-1, to=1, by=0.25)
  lapply(testValues, function(x) { expect_equal(fisherRtoZ(x), 0.5 * log((1 + x) / (1 - x))) })
})

test_that("scaleAB adjusts a vector to the correct scale", {
  expect_equal(scaleAB(seq(from=0, to=10), 0, 1), 
               c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
})

test_that("sumOfPowers computes the power series of a matrix", {
  require(expm)
  A <- matrix(abs(rnorm(10*10)), nrow=10, ncol=10)
  g <- A + A%^%2 + A%^%3 + A%^%4
  expect_equal(sumOfPowers(A, 4), g)
})