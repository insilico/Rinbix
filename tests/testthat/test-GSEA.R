library(Rinbix)
context("GSEA")

test_that("Map SNPs to Genes Biomart", {
  snp.info <- mapSNPsToGenesBiomart(c("rs1048194", "rs13126272"))
  expect_equal(nrow(snp.info), 2)
})

test_that("Lookup Gene Descriptions Biomart", {
  gene.info <- lookupGeneDescBiomart(c("RXRA", "ATMIN"))
  expect_equal(nrow(gene.info), 2)
})

test_that("getReactomePathways", {
	data(geneListSymbols)
	reactomePathDesc <- getReactomePathways(geneListSymbols)
  expect_equal(class(reactomePathDesc)[1], "enrichResult")
})

test_that("getGOAnalysis", {
	data(geneListSymbols)
	goEnrichment <- getGOAnalysis(geneListSymbols)
  expect_equal(class(goEnrichment)[1], "enrichResult")
})
