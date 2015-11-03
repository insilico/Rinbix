library(Rinbix)
data(geneListSymbols)
context("GSEA")

test_that("getReactomePathways", {
	reactomePathDesc <- getReactomePathways(geneListSymbols)
  expect_equal(class(reactomePathDesc)[1], "enrichResult")
})

test_that("getKEGGAnalysis", {
	keggEnrichment <- getKEGGAnalysis(geneListSymbols)
  expect_equal(class(keggEnrichment)[1], "enrichResult")
})

test_that("getGOAnalysis", {
	goEnrichment <- getGOAnalysis(geneListSymbols)
  expect_equal(class(goEnrichment)[1], "enrichResult")
})
