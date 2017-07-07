library(Rinbix)
data(testdata10)
context("GAIN tests vs C++ interface")

test_that("R reGAIN same as inbix C++ reGAIN", {
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  rinbixCppRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)$reGAIN
  expect_equal(object=rinbixRegain, expected=rinbixCppRegain, tolerance=0.02)
})

test_that("R dcGAIN same as inbix C++ dcGAIN", {
  rinbixDcgain <- dcgain(testdata10)
  rinbixCppDcgain <- dcgainInbix(testdata10)
  expect_equal(object=rinbixDcgain, expected=rinbixCppDcgain, tolerance=0.02)
})

test_that("R modularity same as inbix C++ modularity", {
	corMatrix <- cor(testdata10[, -ncol(testdata10)])
	rinbixModulesDF <- Rinbix::modularity(corMatrix)
	moduleListGrps <- as.data.frame(rinbixModulesDF$groups)
	moduleListGrps$Gene <- as.character(moduleListGrps$Gene)
	moduleListGrps$Group <- as.integer(moduleListGrps$Group)
	moduleListGrpsKey <- as.integer(substr(moduleListGrps$Gene, 4, length(moduleListGrps$Gene)))
	moduleListGrps <- moduleListGrps[order(moduleListGrpsKey), ]
	rinbixCppModulesDF <- modularityInbix(corMatrix)
	moduleListGrpsCpp <- data.frame(Gene=as.character(rinbixCppModulesDF$Var),
	                                Group=as.integer(rinbixCppModulesDF$Module))
	moduleListGrpsCppKey <- as.integer(substr(moduleListGrpsCpp$Gene, 4, length(moduleListGrpsCpp$Gene)))
	moduleListGrpsCpp <- moduleListGrpsCpp[order(moduleListGrpsCppKey), ]
	expect_equal(fossil::rand.index(moduleListGrpsCpp$Group, moduleListGrps$Group), 1)
})

test_that("reGAIN stdBetas=TRUE, absBetas=TRUE", {
  inbixRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)$reGAIN
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=FALSE, absBetas=TRUE", {
  inbixRegain <- regainInbix(testdata10, stdBetas=FALSE, absBetas=TRUE)$reGAIN
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=TRUE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=TRUE, absBetas=FALSE", {
  inbixRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=FALSE)$reGAIN
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=FALSE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("reGAIN stdBetas=FALSE, absBetas=FALSE", {
  inbixRegain <- regainInbix(testdata10, stdBetas=FALSE, absBetas=FALSE)$reGAIN
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=FALSE)
  expect_equal(object=rinbixRegain, expected=inbixRegain, tolerance=0.02)
})

test_that("SNPrank Rinbix vs C++ from reGAIN stdBetas=TRUE, absBetas=TRUE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  inbixSnpranksDF <- snprankInbix(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[, 2]
  expect_equal(object=rinbixSnpranks, expected=inbixSnpranks, tolerance=0.05)
})

test_that("SNPrank Rinbix vs C++ from reGAIN stdBetas=TRUE, absBetas=TRUE", {
  inbixRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=FALSE)
  inbixSnpranksDF <- snprankInbix(inbixRegain$reGAIN)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(inbixRegain$reGAIN)
  rinbixSnpranks <- rinbixSnpranksDF[, 2]
  expect_equal(object=rinbixSnpranks, expected=inbixSnpranks, tolerance=0.05)
})

test_that("Read/write inbix data files", {
	writeRegressionDataAsInbixNumeric(testdata10, "foo")
	foodata10 <- readInbixNumericAsRegressionData("foo")
	file.remove(c("foo.num", "foo.pheno"))
  expect_equal(object=foodata10, expected=testdata10, tolerance=0.05)
})

test_that("Permute GAIN inbix", {
  inbixRegainThresholds <- permuteGainInbix(testdata10)
  expect_equal(object=nrow(inbixRegainThresholds), expected=ncol(testdata10)-1)
})
