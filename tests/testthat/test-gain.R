library(Rinbix)
data("testdata10")
data("testdata100ME4")

context("GAIN")

test_that("dcGAIN", {
  inbixDcgain <- dcgainInbix(testdata100ME4)
  rinbixDcgain <- dcgain(testdata100ME4)
  expect_equal(inbixDcgain$scores, rinbixDcgain$scores, tolerance = 0.002)
})

test_that("dmGAIN", {
  inbixDmgain <- as.matrix(read.table("testdata10.dmgain", header = T))
  rinbixDmgain <- dmgain(testdata10)
  expect_equal(inbixDmgain, rinbixDmgain$scores, tolerance = 0.002)
})

test_that("R reGAIN front end same as inbix C++ reGAIN", {
  rinbixRegain <- regain(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)
  rinbixCppRegain <- regainInbix(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)$reGAIN
  expect_equal(object = rinbixRegain, expected = rinbixCppRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = TRUE, absBetas = TRUE", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)
  rinbixCppRegain <- regainInbix(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)$reGAIN
  expect_equal(object = rinbixRegain, expected = rinbixCppRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = FALSE, absBetas = TRUE", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = FALSE, absBetas = TRUE)
  rinbixCppRegain <- regainInbix(testdata100ME4, stdBetas = FALSE, absBetas = TRUE)$reGAIN
  expect_equal(object = rinbixRegain, expected = rinbixCppRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = TRUE, absBetas = FALSE", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = TRUE, absBetas = FALSE)
  rinbixCppRegain <- regainInbix(testdata100ME4, stdBetas = TRUE, absBetas = FALSE)$reGAIN
  expect_equal(object = rinbixRegain, expected = rinbixCppRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = FALSE, absBetas = FALSE", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = FALSE, absBetas = FALSE)
  rinbixCppRegain <- regainInbix(testdata100ME4, stdBetas = FALSE, absBetas = FALSE)$reGAIN
  expect_equal(object = rinbixRegain, expected = rinbixCppRegain, tolerance = 0.02)
})

test_that("GAIN to Simple SIF", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)
  gainSIF <- gainToSimpleSIF(rinbixRegain)
  expect_equal(nrow(gainSIF), choose(nrow(rinbixRegain), 2))
})
