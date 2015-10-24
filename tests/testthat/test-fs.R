
library(Rinbix)
data(testdata10)
data(testdata100ME4)
context("Feature Selection")

test_that("GeneRank", {
	success <- FALSE
	geneRankResults <- rankGeneRank(testdata100ME4)
	if(nrow(geneRankResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("dcGAIN + SNPrank", {
	success <- FALSE
	dcgainResults <- rankDcgainSnprank(testdata100ME4)
	if(nrow(dcgainResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("GeneRank", {
	success <- FALSE
	rankGeneRankResults <- rankGeneRank(testdata100ME4)
	if(nrow(rankGeneRankResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("glmnet", {
	success <- TRUE
	rankGlmnetResults <- rankGlmnet(testdata100ME4)
	if(nrow(rankGlmnetResults) == 2) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("Lasso", {
	success <- TRUE
	rankLassoResults <- rankLasso(testdata100ME4)
		if(nrow(rankLassoResults) == 2) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("Limma", {
	success <- FALSE
	rankLimmaResults <- rankLimma(testdata100ME4)
	if(nrow(rankLimmaResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("RandomForest", {
	success <- FALSE
	rankRandomForestResults <- rankRandomForest(testdata100ME4)
	if(nrow(rankRandomForestResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("ReliefSeq", {
	success <- FALSE
	rankReliefSeqResults <- rankReliefSeq(testdata100ME4)
	if(nrow(rankReliefSeqResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("SAM", {
	success <- FALSE
	rankSamResults <- rankSam(testdata100ME4)
	if(nrow(rankSamResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("t-test", {
	success <- FALSE
	rankTTestResults <- rankTTest(testdata100ME4)
	if(length(rankTTestResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

test_that("univariate regression", {
	success <- FALSE
	rankUnivariateRegressionResults <- rankUnivariateRegression(testdata100ME4)
	if(nrow(rankUnivariateRegressionResults) == (ncol(testdata100ME4) - 1)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})

# snprank
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

# r2VIM
test_that("r2VIM", {
	success <- FALSE
	predictors <- as.matrix(testdata100ME4[, -ncol(testdata100ME4)])
	response <- factor(testdata100ME4[, ncol(testdata100ME4)])
	r2vimResults <- r2VIMorig(predictors=predictors, response=response, verbose=FALSE)
	if((r2vimResults$votes["gene0001"] == 10) && (r2vimResults$votes["gene0005"] == 10)) {
		success <- TRUE
	}
	expect_equal(success, TRUE)
})
