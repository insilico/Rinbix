context("Predict Expression")

test_that("", {
  require(Rinbix)
  data(testdata10)
  expect_that(rinbixRegain, equals(rinbixCppRegain$reGAIN, tolerance=0.002))
})
