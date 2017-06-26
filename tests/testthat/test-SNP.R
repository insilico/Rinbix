library(Rinbix)
context("SNP")

test_that("readPlinkBinary", {
  plinkObj <- readPlinkBinary("wgas1")
  expect_equal(nrow(plinkObj$genotypes), 90)
  expect_equal(ncol(plinkObj$genotypes), 228694)
  expect_equal(nrow(plinkObj$fam), 90)
  expect_equal(ncol(plinkObj$fam), 6)
  expect_equal(nrow(plinkObj$map), 228694)
  expect_equal(ncol(plinkObj$map), 6)
})
