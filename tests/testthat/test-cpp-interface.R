context("GAIN tests vs C++ interface")

test_that("R reGAIN same as C++", {
  require(Rinbix)
  data(testdata10)
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  #print(rinbixRegain)
  rinbixCppRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)
  #print(rinbixCppRegain)
  expect_that(rinbixRegain, equals(rinbixCppRegain$reGAIN, tolerance=0.002))
})
