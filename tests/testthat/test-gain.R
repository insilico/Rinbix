context("Genetic association interaction network (GAIN) tests")

test_that("dcGAIN", {
  require(Rinbix)
  data(testdata10)
  inbixDcgain <- as.matrix(read.table("testdata10.dcgain", header=T))
  rinbixDcgain <- dcgain(testdata10)
  expect_that(inbixDcgain, equals(rinbixDcgain$scores, tolerance=0.002))
})

test_that("dmGAIN", {
  require(Rinbix)
  data(testdata10)
  inbixDmgain <- as.matrix(read.table("testdata10.dmgain", header=T))
  rinbixDmgain <- dmgain(testdata10)
  expect_that(inbixDmgain, equals(rinbixDmgain$scores, tolerance=0.002))
})

test_that("reGAIN Workflow", {
  require(Rinbix)
  data(testdata10)
  inbixRegain <- as.matrix(read.table("testdata10.block.regain", header=T))
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  expect_that(all(inbixRegain-rinbixRegain < 0.005), is_true())
  
  inbixSnpranksDF <- read.table("testdata10.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(inbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_that(all(inbixSnpranks-rinbixSnpranks < 0.005), is_true())

  inbixModulesDF <- read.table("testdata10.modules", header=F)
  rinbixModulesDF <- Rinbix::modularity(inbixRegain)
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF[,2])
  expect_that(all(inbixModules == rinbixModules), is_true())
})

test_that("R reGAIN same as C++", {
  require(Rinbix)
  data(testdata10)
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  print(rinbixRegain)
  rinbixCppRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)
  print(rinbixCppRegain)
  expect_that(rinbixRegain, equals(rinbixCppRegain$reGAIN, tolerance=0.002))
})
