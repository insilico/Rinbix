library(Rinbix)
data(testdata10)
data(testdata100ME4)
context("Feature Selection")

test_that("GeneRank", {
	geneRankResults <- rankGeneRank(testdata100ME4)
	expect_equal(ncol(geneRankResults), 2)
	expect_equal(nrow(geneRankResults), ncol(testdata100ME4) - 1)
})

test_that("dcGAIN + SNPrank", {
	dcgainResults <- rankDcgainSnprank(testdata100ME4)
	expect_equal(ncol(dcgainResults), 2)
	expect_equal(nrow(dcgainResults), ncol(testdata100ME4) - 1)
})

test_that("glmnet", {
	glmnetResults <- rankGlmnet(testdata100ME4)
	expect_equal(ncol(glmnetResults), 2)
	expect_equal(nrow(glmnetResults), 2)
})

test_that("Iterative Relief", {
  bestReliefResults <- rankIterativeRelieff(testdata100ME4)
  expect_equal(ncol(bestReliefResults), 2)
  expect_equal(nrow(bestReliefResults), ncol(testdata100ME4) - 1)
})

test_that("Lasso", {
	lassoResults <- rankLasso(testdata100ME4)
	expect_equal(ncol(lassoResults), 2)
	expect_equal(nrow(lassoResults), 2)
})

test_that("Limma", {
	limmaResults <- rankLimma(testdata100ME4)
	expect_equal(ncol(limmaResults), 2)
	expect_equal(nrow(limmaResults), ncol(testdata100ME4) - 1)
})

test_that("RandomForest", {
	randomForestResults <- rankRandomForest(testdata100ME4)
	expect_equal(ncol(randomForestResults), 2)
	expect_equal(nrow(randomForestResults), ncol(testdata100ME4) - 1)
})

test_that("Relief", {
  reliefResults <- rankRelieff(testdata100ME4)
  expect_equal(ncol(reliefResults), 2)
  expect_equal(nrow(reliefResults), ncol(testdata100ME4) - 1)
})

test_that("t-test", {
	tTestResults <- rankTTest(testdata100ME4)
	expect_equal(ncol(tTestResults), 2)
	expect_equal(nrow(tTestResults), ncol(testdata100ME4) - 1)
})

test_that("univariate regression", {
	univariateRegressionResults <- rankUnivariateRegression(testdata100ME4)
	expect_equal(ncol(univariateRegressionResults), 2)
	expect_equal(nrow(univariateRegressionResults), ncol(testdata100ME4) - 1)
})

# reGAIN parameter combinations + snprank
test_that("Rinbix SNPrank from reGAIN stdBetas = TRUE, absBetas = TRUE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas = TRUE, absBetas = TRUE)
  inbixSnpranksDF <- snprankInbix(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[, 2]
  expect_equal(object = rinbixSnpranks, expected = inbixSnpranks, tolerance = 0.05)
})

test_that("Rinbix SNPrank from reGAIN stdBetas = FALSE, absBetas = TRUE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas = FALSE, absBetas = TRUE)
  inbixSnpranksDF <- snprankInbix(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object = rinbixSnpranks, expected = inbixSnpranks, tolerance = 0.05)
})

test_that("Rinbix SNPrank from reGAIN stdBetas = TRUE, absBetas = FALSE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas = TRUE, absBetas = FALSE)
  inbixSnpranksDF <- snprankInbix(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object = rinbixSnpranks, expected = inbixSnpranks, tolerance = 0.05)
})

test_that("Rinbix SNPrank from reGAIN stdBetas = FALSE, absBetas = FALSE", {
  rinbixRegain <- regainParallel(testdata10, stdBetas = FALSE, absBetas = FALSE)
  inbixSnpranksDF <- snprankInbix(rinbixRegain)
  inbixSnpranks <- inbixSnpranksDF[, 2]
  rinbixSnpranksDF <- snprank(rinbixRegain)
  rinbixSnpranks <- rinbixSnpranksDF[,2]
  expect_equal(object = rinbixSnpranks, expected = inbixSnpranks, tolerance = 0.05)
})

# r2VIM
test_that("r2VIM", {
	predictors <- as.matrix(testdata100ME4[, -ncol(testdata100ME4)])
	response <- factor(testdata100ME4[, ncol(testdata100ME4)])
	r2vimResults <- r2VIMorig(predictors = predictors, response = response, 
	                          verbose = FALSE)
	expect_equal(((r2vimResults$votes["gene0001"] == 10) && 
	              (r2vimResults$votes["gene0005"] == 10)), TRUE)
})
