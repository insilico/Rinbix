context("Fisher R to Z Transformation tests")
  test_that("Check range from -1 to 1", {
    require(Rinbix)
    testValues <- seq(from=-1, to=1, by=0.25)
    lapply(testValues, function(x) { expect_that(fisherRtoZ(x), equals(0.5 * log((1+x)/(1-x)))) })
  })
