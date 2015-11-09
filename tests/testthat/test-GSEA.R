library(Rinbix)
context("GSEA")

test_that("getReactomePathways", {
	data(geneListSymbols)
	reactomePathDesc <- getReactomePathways(geneListSymbols)
  expect_equal(class(reactomePathDesc)[1], "enrichResult")
})

test_that("getKEGGAnalysis", {
	data(geneListSymbols)
	keggEnrichment <- getKEGGAnalysis(geneListSymbols)
  expect_equal(class(keggEnrichment)[1], "enrichResult")
})

test_that("getGOAnalysis", {
	data(geneListSymbols)
	goEnrichment <- getGOAnalysis(geneListSymbols)
  expect_equal(class(goEnrichment)[1], "enrichResult")
})
