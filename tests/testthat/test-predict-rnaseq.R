library(Rinbix)
data(simrnaseq)  
context("Predict RNA-Seq")

test_that("Simulated RNASeq", {
  cat("FOO!\n")
  predictResult <- predictRnaseq(rnaseqCountsTrain=predictorsTrain, 
                                 groupLabelsTrain=responseTrain, 
                                 rnaseqCountsTest=predictorsTest, 
                                 groupLabelsTest=responseTest, 
                                 preprocessMethod="none", 
                                 filterMethod="randomforests", 
                                 topN=2, 
                                 classifierMethod="svm",
                                 verbose=TRUE)
  expect_equal(TRUE, TRUE)
})
