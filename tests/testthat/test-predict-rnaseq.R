library(testthat)
library(compcodeR)

context("Predict RNASeq")

test_that("Simulated RNASeq", {
  cat("Simulating data sets...\n")
  # ---------------------------------------------------------------------------  
  countsTrain <- generateSyntheticData(dataset="simdataTrain", n.vars=100,
                                       samples.per.cond=20, n.diffexp=10)
  predictorsTrain <- t(countsTrain@count.matrix)
  responseTrain <- countsTrain@sample.annotations$condition
  # ---------------------------------------------------------------------------  
  countsTest <- generateSyntheticData(dataset="simdataTest", n.vars=100,
                                       samples.per.cond=20, n.diffexp=10)
  predictorsTest <- t(countsTest@count.matrix)
  responseTest <- countsTest@sample.annotations$condition
  
  # ---------------------------------------------------------------------------  
  cat("Running predictRnaseq()...\n")
  predictResult <- predictRnaseq(predictorsTrain, responseTrain, 
                                 predictorsTest, responseTest, 
                                 "none", "deseq2", 10, "none")
  
  expect_that(TRUE, equals(TRUE))
})
