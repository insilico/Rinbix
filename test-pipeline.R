context("Genetic association interaction network (GAIN) Pipeline tests")

test_that("SNPrank", {
  require(Rinbix)
  data(testdata10)
  inbixSnpranksDF <- read.table("testdata10.ranks", header=T)
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  rinbixSnpranksDF <- snprank(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_that(inbixSnpranks, equals(rinbixSnpranks, tolerance=0.002))
})

test_that("Modularity", {
  require(Rinbix)
  data(testdata10)
  inbixModulesDF <- read.table("testdata10.modules", header=F)
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  rinbixModulesDF <- Rinbix::modularity(rinbixRegain)
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF[,2])
  expect_that(inbixModules, equals(rinbixModules))
})
