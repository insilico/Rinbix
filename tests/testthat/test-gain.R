context("Genetic association interaction network (GAIN) tests")

test_that("dcGAIN", {
  require(Rinbix)
  data(testdata10)
  inbixDcgain <- as.matrix(read.table("testdata10.dcgain", header=T))
  rinbixDcgain <- dcgain(testdata10)
  expect_that(inbixDcgain, equals(rinbixDcgain, tolerance=0.002))
})

test_that("dmGAIN", {
  require(Rinbix)
  data(testdata10)
  inbixDmgain <- as.matrix(read.table("testdata10.dmgain", header=T))
  rinbixDmgain <- dmgain(testdata10)
  expect_that(inbixDmgain, equals(rinbixDmgain, tolerance=0.002))
})

test_that("reGAIN Workflow", {
  require(Rinbix)
  data(testdata10)
  inbixRegain <- as.matrix(read.table("testdata10.block.regain", header=T))
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  expect_that(inbixRegain, equals(rinbixRegain, tolerance=0.002))

  inbixSnpranksDF <- read.table("testdata10.ranks", header=T)
  rinbixSnpranksDF <- snprank(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_that(inbixSnpranks, equals(rinbixSnpranks, tolerance=0.002))

  inbixModulesDF <- read.table("testdata10.modules", header=F)
  rinbixModulesDF <- Rinbix::modularity(rinbixRegain)
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF[,2])
  expect_that(inbixModules, equals(rinbixModules))
})
