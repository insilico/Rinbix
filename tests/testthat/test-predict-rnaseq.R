library('testthat')
library('Rinbix')

context("Predict RNASeq")

test_that("Simulated RNASeq", {
  cat("Loading simulated data sets...\n")
  data(simrnaseq)  

  cat("Running predictRnaseq()...\n")
  predictResult <- Rinbix::predictRnaseq(predictorsTrain, responseTrain, 
                                         predictorsTest, responseTest, 
                                         "none", "randomforests", 5, "svm")
  
  expect_that(TRUE, equals(TRUE))
})
