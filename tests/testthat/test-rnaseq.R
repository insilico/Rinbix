library(Rinbix)
data(simrnaseq)
context("RNA-Seq")

# test_that("DESeq", {
# 	require(DESeq)
# 	data(simrnaseq)
# 	X <- t(predictorsTrain)
# 	y <- as.factor(responseTrain - 1)
# 	rDeseq <- rankDeseq(X, y)
#   expect_equal(nrow(rDeseq), ncol(predictorsTrain))
# })

test_that("DESeq2", {
	require(DESeq2)
	data(simrnaseq)
	X <- t(predictorsTrain)
	y <- as.factor(responseTrain - 1)
	rDeseq2 <- rankDeseq2(X, y)
  expect_equal(nrow(rDeseq2), ncol(predictorsTrain))
})

test_that("edgeR", {
	require(edgeR)
	data(simrnaseq)
	X <- t(predictorsTrain)
	y <- as.factor(responseTrain - 1)
	rEdger <- rankEdgeR(X, y)
  expect_equal(nrow(rEdger), ncol(predictorsTrain))
})

test_that("edgeR GLM", {
	require(edgeR)
	data(simrnaseq)
	X <- t(predictorsTrain)
	y <- as.factor(responseTrain - 1)
	covars <- rnorm(ncol(X))
  rEdgerGlm <- rankEdgeRGlm(X, y, covars)
  expect_equal(nrow(rEdgerGlm), ncol(predictorsTrain))
})
