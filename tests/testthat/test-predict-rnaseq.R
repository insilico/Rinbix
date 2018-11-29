library(Rinbix)
data(simrnaseq)
context("Predict RNA-Seq")

test_that("predictRnaseq", {
	set.seed(1965)
	predictResult <- predictRnaseq(rnaseqCountsTrain = predictorsTrain, 
                                 groupLabelsTrain = responseTrain, 
                                 rnaseqCountsTest = predictorsTest, 
                                 groupLabelsTest = responseTest, 
                                 preprocessMethod = "none", 
                                 filterMethod = "randomforests", 
                                 topN = 2, 
                                 classifierMethod = "svm",
                                 verbose = FALSE)
  expect_equal(predictResult$stats$statsTrain$TP, 19)
  expect_equal(predictResult$stats$statsTest$TP, 0)
})

test_that("preprocessRnaseq", {
	preprocessResult <- preprocessRnaseq(method = "log2", 
	                                     predictorsTrain, 
	                                     predictorsTest, 
	                                     verbose = FALSE)
	expect_equal(length(preprocessResult), 2)
})

test_that("filterRnaseq", {
  filteredGenes <- filterRnaseq(method = "randomforests", 
                                predictorsTrain, 
                                responseTrain, 
                                predictorsTest, 
                                responseTest,
                                nTopGenes = 10, 
                                verbose = FALSE)
  expect_equal(length(filteredGenes), 2)
})

test_that("classifyRnaseq", {
  classifyStats <- classifyRnaseq(method = "svm", 
                                  predictorsTrain, 
                                  responseTrain, 
                                  predictorsTest, 
                                  responseTest,
                                  verbose = FALSE)
  expect_equal(length(classifyStats), 2)
})
