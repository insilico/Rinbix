library(Rinbix)
data(simrnaseq)
context("Predict RNA-Seq")

test_that("predictRnaseq", {
	data(simrnaseq)
	set.seed(1965)
	predictResult <- predictRnaseq(rnaseqCountsTrain=predictorsTrain, 
                                 groupLabelsTrain=responseTrain, 
                                 rnaseqCountsTest=predictorsTest, 
                                 groupLabelsTest=responseTest, 
                                 preprocessMethod="none", 
                                 filterMethod="randomforests", 
                                 topN=2, 
                                 classifierMethod="svm",
                                 verbose=FALSE)
  expect_equal(predictResult$stats$statsTrain$TP, 19)
  expect_equal(predictResult$stats$statsTest$TP, 0)
})

test_that("preprocess", {
	preprocessResult <- preprocess(method="log2", 
	                               predictorsTrain, 
	                               predictorsTest, 
	                               verbose=FALSE)
	expect_equal(length(preprocessResult), 2)
})

test_that("filterGenes", {
	filteredGenes <- filterGenes(method="randomforests", 
	                             predictorsTrain, 
	                             responseTrain, 
	                             predictorsTest, 
	                             responseTest,
	                             nTopGenes=10, 
	                             verbose=FALSE)
	expect_equal(length(filteredGenes), 2)
})

test_that("classify", {
	classifyStats <- classify(method="svm", 
	                          predictorsTrain, 
	                          responseTrain, 
	                          predictorsTest, 
	                          responseTest,
	                          verbose=FALSE)
	expect_equal(length(classifyStats), 2)
})
