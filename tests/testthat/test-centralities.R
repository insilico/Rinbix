library(Rinbix)
context("Centralities")

test_that("EpistasisRank", {
  data(testdata10)
  rinbixRegain <- regainParallel(testdata10, stdBetas = TRUE, absBetas = TRUE)
  rinbixER_DF <- EpistasisRank(rinbixRegain)
  expect_equal(ncol(rinbixER_DF), 2)
  expect_equal(nrow(rinbixER_DF), ncol(testdata10) - 1)
})

test_that("EpistasisKatz", {
  data(testdata10)
  rinbixRegain <- regainParallel(testdata10, stdBetas = TRUE, absBetas = TRUE)
  alpha <- 1 / mean(colSums(rinbixRegain))
  beta <- diag(rinbixRegain)
  rinbixEK_DF <- EpistasisKatz(rinbixRegain, alpha, beta)
  expect_equal(length(rinbixEK_DF), ncol(testdata10) - 1)
})

test_that("PageRank", {
  data(testdata10)
  adj <- ifelse((cor(testdata10)) > 0, 1, 0)
  damping <- 0.85
  myPageRank <- PageRank(adj, damping)
  expect_equal(nrow(myPageRank), ncol(testdata10))
})
