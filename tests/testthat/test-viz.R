library(Rinbix)
context("Visualization")

test_that("adjacencyMatrixToChordDiagram", {
  data("testdata10")
  predictors <- testdata10[, -ncol(testdata10)]
  Acorr <- cor(predictors)
  chord <- adjacencyMatrixToChordDiagram(Acorr)
  expect_equal(TRUE, TRUE)
})

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

test_that("igraphToChordDiagram", {
  data("testdata10")
  predictors <- testdata10[, -ncol(testdata10)]
  Acorr <- cor(predictors)
  netlist <- adjacencyToNetList(Acorr, 
                                thresholdType="hard", 
                                thresholdValue=0.2, 
                                useAbs=TRUE, 
                                useWeighted=TRUE)
  chord <- igraphToChordDiagram(netlist$net, netlist$nodes$Name)
  expect_equal(TRUE, TRUE)
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
  # convert nodes to zero based indexes
  netlinks <- data.frame(Source=netlist$links$Source - 1,
                         Target=netlist$links$Target - 1,
                         Value=netlist$links$Value)
  netnodes <- data.frame(NodeID=netlist$nodes$NodeID -1,
                         Name=netlist$nodes$Name,
                         Group=netlist$nodes$Group,
                         Size=netlist$nodes$Size)
  net <- netListToSimpleD3(list(links=netlinks, nodes=netnodes, groups=NULL))
  expect_equal(length(net), 8)
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
  data("testdata10")
  predictors <- testdata10[, -ncol(testdata10)]
  Acorr <- cor(predictors)
  netlist <- adjacencyToNetList(Acorr,
                                thresholdType="hard",
                                thresholdValue=0.2,
                                useAbs=TRUE,
                                useWeighted=TRUE)
  # convert nodes to zero based indexes
  netlinks <- data.frame(Source=netlist$links$Source - 1,
                         Target=netlist$links$Target - 1,
                         Value=netlist$links$Value)
  netnodes <- data.frame(NodeID=netlist$nodes$NodeID -1,
                         Name=netlist$nodes$Name,
                         Group=netlist$nodes$Group,
                         Size=netlist$nodes$Size)
  net <- netListToForceD3(list(links=netlinks, nodes=netnodes, groups=NULL))
  expect_equal(length(net), 8)
})
