context("GAIN tests vs Standalone")

test_that("dcGAIN", {
  require(Rinbix)
  data(testdata10)
  inbixDcgain <- as.matrix(read.table("testdata10.dcgain", header=T))
  inbixDcgainP <- as.matrix(read.table("testdata10.pvals.dcgain", header=T))
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

test_that("reGAIN Workflow, stdBetas=TRUE, absBetas=TRUE", {
  require(Rinbix)
  data(testdata10)
  inbixRegain <- as.matrix(read.table("testdata10-abs-zval.block.regain", header=T))
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  expect_that(inbixRegain, equals(rinbixRegain, tolerance=0.002))
  
  inbixSnpranksDF <- read.table("testdata10-abs-zval.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(inbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_that(all(inbixSnpranks-rinbixSnpranks < 0.005), is_true())

  inbixModulesDF <- read.table("testdata10-abs-zval.modules", header=F)
  rinbixModulesDF <- Rinbix::modularity(inbixRegain)
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF$groups[,2])
  expect_that(all(inbixModules == rinbixModules), is_true())
})

test_that("reGAIN Workflow, stdBetas=FALSE, absBetas=TRUE", {
  require(Rinbix)
  data(testdata10)
  inbixRegain <- as.matrix(read.table("testdata10-abs-beta.block.regain", header=T))
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=TRUE)
  expect_that(inbixRegain, equals(rinbixRegain, tolerance=0.002))
  
  inbixSnpranksDF <- read.table("testdata10-abs-beta.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(inbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_that(all(inbixSnpranks-rinbixSnpranks < 0.005), is_true())
  
  inbixModulesDF <- read.table("testdata10-abs-beta.modules", header=F)
  rinbixModulesDF <- Rinbix::modularity(inbixRegain)
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF[,2])
  expect_that(all(inbixModules == rinbixModules), is_true())
})

test_that("reGAIN Workflow, stdBetas=TRUE, absBetas=FALSE", {
  require(Rinbix)
  data(testdata10)
  inbixRegain <- as.matrix(read.table("testdata10-noabs-zval.block.regain", header=T))
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=FALSE)
  expect_that(inbixRegain, equals(rinbixRegain, tolerance=0.002))
  
  inbixSnpranksDF <- read.table("testdata10-noabs-zval.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(inbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_that(all(inbixSnpranks-rinbixSnpranks < 0.005), is_true())
  
  inbixModulesDF <- read.table("testdata10-noabs-zval.modules", header=F)
  rinbixModulesDF <- Rinbix::modularity(inbixRegain)
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF[,2])
  expect_that(all(inbixModules == rinbixModules), is_true())
})

test_that("reGAIN Workflow, stdBetas=FALSE, absBetas=FALSE", {
  require(Rinbix)
  data(testdata10)
  inbixRegain <- as.matrix(read.table("testdata10-noabs-beta.block.regain", header=T))
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=FALSE)
  expect_that(inbixRegain, equals(rinbixRegain, tolerance=0.002))
  
  inbixSnpranksDF <- read.table("testdata10-noabs-beta.ranks", header=T)
  inbixSnpranks <- inbixSnpranksDF[,2]
  rinbixSnpranksDF <- snprank(inbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  print(inbixSnpranks)
  print(rinbixSnpranks)
  expect_that(all(inbixSnpranks-rinbixSnpranks < 0.005), is_true())
  
  inbixModulesDF <- read.table("testdata10-noabs-beta.modules", header=F)
  colnames(inbixModulesDF) <- c("Gene", "Group")
  rinbixModulesDF <- as.data.frame(Rinbix::modularity(inbixRegain))
  inbixModules <- inbixModulesDF[,2]
  rinbixModules <- as.integer(rinbixModulesDF[,2])
  print(inbixModulesDF)
  print(rinbixModulesDF)
  expect_that(all(inbixModules == rinbixModules), is_true())
})
