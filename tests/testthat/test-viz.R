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
	                              useWeighted=TRUE, 
	                              verbose=FALSE, 
	                              groups=NULL)
	filteredNetwork <- netlist$net
  expect_equal(nrow(netlist$nodes), 9)
  expect_equal(nrow(netlist$links), 12)
})
