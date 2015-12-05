library(Rinbix)
context("Visualization")

test_that("adjacencyToNetList", {
	data("testdata10")
	predictors <- testdata10[, -ncol(testdata10)]
	Acorr <- cor(predictors)
	netlist <- adjacencyToNetList(Acorr, 
	                              thresholdType="hard", 
	                              thresholdValue=0.2, 
	                              useAbs=TRUE, 
	                              useWeighted=TRUE, verbose=FALSE)
	filteredNetwork <- netlist$net
  expect_equal(nrow(netlist$nodes), 10)
  expect_equal(nrow(netlist$links), 23)
})

test_that("netListToSimpleD3", {
  data("testdata10")
  predictors <- testdata10[, -ncol(testdata10)]
  Acorr <- cor(predictors)
  netlist <- adjacencyToNetList(Acorr, 
                                thresholdType="hard", 
                                thresholdValue=0.2, 
                                useAbs=TRUE, 
                                useWeighted=TRUE)
  net <- netListToSimpleD3(netlist)
  expect_equal(length(net), 7)
})

test_that("netListToForceD3", {
  data("testdata10")
  predictors <- testdata10[, -ncol(testdata10)]
  Acorr <- cor(predictors)
  netlist <- adjacencyToNetList(Acorr, 
                                thresholdType="hard", 
                                thresholdValue=0.2, 
                                useAbs=TRUE, 
                                useWeighted=TRUE)
  net <- netListToForceD3(netlist)
  expect_equal(length(net), 7)
})
