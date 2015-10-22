
library(Rinbix)
data(testdata10)
context("Feature Selection")

test_that("Rinbix SNPrank from reGAIN stdBetas=TRUE, absBetas=TRUE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  inbixSnpranksDF <- read.table("testdata10-abs-zval.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object=rinbixSnpranks, expected=inbixSnpranks, tolerance=0.05)
})

test_that("Rinbix SNPrank from reGAIN stdBetas=FALSE, absBetas=TRUE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=TRUE)
  inbixSnpranksDF <- read.table("testdata10-abs-beta.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object=rinbixSnpranks, expected=inbixSnpranks, tolerance=0.05)
})

test_that("Rinbix SNPrank from reGAIN stdBetas=TRUE, absBetas=FALSE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=FALSE)
  inbixSnpranksDF <- read.table("testdata10-noabs-zval.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object=rinbixSnpranks, expected=inbixSnpranks, tolerance=0.05)
})

test_that("Rinbix SNPrank from reGAIN stdBetas=FALSE, absBetas=FALSE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=FALSE)
  inbixSnpranksDF <- read.table("testdata10-noabs-beta.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object=rinbixSnpranks, expected=inbixSnpranks, tolerance=0.05)
})
