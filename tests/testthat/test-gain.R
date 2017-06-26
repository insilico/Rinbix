library(Rinbix)
data(testdata10)
context("GAIN")

test_that("dcGAIN", {
  inbixDcgain <- as.matrix(read.table("testdata10.dcgain", header=T))
  rinbixDcgain <- dcgain(testdata10)
  expect_equal(inbixDcgain, rinbixDcgain$scores, tolerance=0.002)
})

test_that("dmGAIN", {
  inbixDmgain <- as.matrix(read.table("testdata10.dmgain", header=T))
  rinbixDmgain <- dmgain(testdata10)
  expect_equal(inbixDmgain, rinbixDmgain$scores, tolerance=0.002)
})

test_that("R reGAIN front end same as inbix C++ reGAIN", {
  rinbixRegain <- regain(testdata10, stdBetas=TRUE, absBetas=TRUE)
  rinbixCppRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)$reGAIN
  expect_equal(object=rinbixRegain, expected=rinbixCppRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=TRUE, absBetas=TRUE", {
  inbixRegain <- as.matrix(read.table("testdata10-abs-zval.block.regain", header=T))
  rinbixRegain <- regain(testdata10, stdBetas=TRUE, absBetas=TRUE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=FALSE, absBetas=TRUE", {
  inbixRegain <- as.matrix(read.table("testdata10-abs-beta.block.regain", header=T))
  rinbixRegain <- regain(testdata10, stdBetas=FALSE, absBetas=TRUE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=TRUE, absBetas=FALSE", {
  inbixRegain <- as.matrix(read.table("testdata10-noabs-zval.block.regain", header=T))
  rinbixRegain <- regain(testdata10, stdBetas=TRUE, absBetas=FALSE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=FALSE, absBetas=FALSE", {
  inbixRegain <- as.matrix(read.table("testdata10-noabs-beta.block.regain", header=T))
  rinbixRegain <- regain(testdata10, stdBetas=FALSE, absBetas=FALSE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("GAIN to Simple SIF", {
  rinbixRegain <- regain(testdata10, stdBetas=TRUE, absBetas=TRUE)
  gainSIF <- gainToSimpleSIF(rinbixRegain)
  expect_equal(nrow(gainSIF), choose(nrow(rinbixRegain), 2))
})
