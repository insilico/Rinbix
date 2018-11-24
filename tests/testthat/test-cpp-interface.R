library(Rinbix)
context("Test Rinbix versus C++ interface to inbix")

data(testdata100ME4)
data(testdata10)

test_that("R reGAIN same as inbix C++ reGAIN", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = FALSE, absBetas = TRUE)
  rinbixCppRegainResults <- regainInbix(testdata100ME4, stdBetas = FALSE, 
                                        absBetas = TRUE, verbose = FALSE)
  rinbixCppRegain <- rinbixCppRegainResults$reGAIN
  expect_equal(object = rinbixRegain, expected = rinbixCppRegain, 
               tolerance = 0.02)
})

test_that("R dcGAIN same as inbix C++ dcGAIN", {
  rinbixDcgain <- dcgain(testdata100ME4)
  rinbixCppDcgain <- dcgainInbix(testdata100ME4)
  expect_equal(object = rinbixDcgain$scores, expected = rinbixCppDcgain$scores, 
               tolerance = 0.02)
})

test_that("R modularity same as inbix C++ modularity", {
	corMatrix <- cor(testdata100ME4[, -ncol(testdata100ME4)])
	rinbixModulesDF <- Rinbix::modularity(corMatrix)
	moduleListGrps <- as.data.frame(rinbixModulesDF$groups)
	moduleListGrps$Gene <- as.character(moduleListGrps$Gene)
	moduleListGrps$Group <- as.integer(moduleListGrps$Group)
	moduleListGrpsKey <- as.integer(substr(moduleListGrps$Gene, 5, 
	                                       length(moduleListGrps$Gene)))
	moduleListGrps <- moduleListGrps[order(moduleListGrpsKey), ]
	rinbixCppModulesDF <- modularityInbix(corMatrix)
	moduleListGrpsCpp <- data.frame(Gene = as.character(rinbixCppModulesDF$Var),
	                                Group = as.integer(rinbixCppModulesDF$Module))
	moduleListGrpsCppKey <- as.integer(substr(moduleListGrpsCpp$Gene, 5, 
	                                          length(moduleListGrpsCpp$Gene)))
	moduleListGrpsCpp <- moduleListGrpsCpp[order(moduleListGrpsCppKey), ]
	expect_equal(fossil::rand.index(moduleListGrpsCpp$Group, 
	                                moduleListGrps$Group), 1)
})

test_that("reGAIN stdBetas = TRUE, absBetas = TRUE", {
  inbixRegain <- regainInbix(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)$reGAIN
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)
  expect_equal(object = rinbixRegain, expected = inbixRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = FALSE, absBetas = TRUE", {
  inbixRegain <- regainInbix(testdata100ME4, stdBetas = FALSE, absBetas = TRUE)$reGAIN
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = FALSE, absBetas = TRUE)
  expect_equal(object = rinbixRegain, expected = inbixRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = TRUE, absBetas = FALSE", {
  inbixRegain <- regainInbix(testdata100ME4, stdBetas = TRUE, absBetas = FALSE)$reGAIN
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = TRUE, absBetas = FALSE)
  expect_equal(object = rinbixRegain, expected = inbixRegain, tolerance = 0.02)
})

test_that("reGAIN stdBetas = FALSE, absBetas = FALSE", {
  inbixRegain <- regainInbix(testdata100ME4, stdBetas = FALSE, absBetas = FALSE)$reGAIN
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = FALSE, absBetas = FALSE)
  expect_equal(object = rinbixRegain, expected = inbixRegain, tolerance = 0.02)
})

test_that("SNPrank Rinbix vs C++ from reGAIN stdBetas = TRUE, absBetas = TRUE", {
  rinbixRegain <- regainParallel(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)
  inbixSnpranksDF <- snprankInbix(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[, 2]
  expect_equal(object = rinbixSnpranks, expected = inbixSnpranks, 
               tolerance = 0.05)
})

test_that("SNPrank Rinbix vs C++ from reGAIN stdBetas = TRUE, absBetas = TRUE", {
  # NOTE: does not match when using absBetas = FALSE
  inbixRegain <- regainInbix(testdata100ME4, stdBetas = TRUE, absBetas = TRUE)
  inbixSnpranksDF <- snprankInbix(inbixRegain$reGAIN)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(inbixRegain$reGAIN)
  rinbixSnpranks <- rinbixSnpranksDF[, 2]
  expect_equal(object = rinbixSnpranks, expected = inbixSnpranks, 
               tolerance = 0.05)
})

test_that("Read/write inbix data files", {
	writeRegressionDataAsInbixNumeric(testdata100ME4, "foo")
	foodata100ME4 <- readInbixNumericAsRegressionData("foo")
	rownames(foodata100ME4) <- rownames(testdata100ME4)
	file.remove(c("foo.num", "foo.pheno"))
  expect_equal(object = foodata100ME4, expected = testdata100ME4, 
               tolerance = 0.05)
})

test_that("Permute GAIN inbix", {
  inbixRegainThresholds <- permuteGainInbix(testdata100ME4, numPerms = 5, 
                                            method = 'regain')
  expect_equal(object = nrow(inbixRegainThresholds),
               expected = ncol(testdata100ME4) - 1)
})

test_that("Permute dcgain inbix", {
  inbixRegainThresholds <- permuteGainInbix(testdata100ME4, numPerms = 5, 
                                            method = 'dcgain')
  expect_equal(object = nrow(inbixRegainThresholds),
               expected = ncol(testdata100ME4) - 1)
})
