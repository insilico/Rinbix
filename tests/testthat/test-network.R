library(Rinbix)
context("Network")

test_that("getIGraphStats", {
	require(igraph)
	g <- erdos.renyi.game(1000, 1/1000)
	igstats <- getIGraphStats(g)
  expect_equal(igstats$nnodes, 1000)
})

test_that("randomNetworkSim", {
	net <- randomNetworkSim(n=1000)
  expect_equal(ncol(net), 1000)
  expect_equal(nrow(net), 1000)
})

test_that("scaleFreeSim", {
	net <- scaleFreeSim(n=1000)
  expect_equal(ncol(net), 1000)
  expect_equal(nrow(net), 1000)
})
