## ---- eval=FALSE---------------------------------------------------------
#  install.packages("path_to_file/Rinbix-VERSIONSTRING.tar.gz",
#                   repos = NULL, type = "source")

## ---- eval=FALSE---------------------------------------------------------
#  R CMD INSTALL Rinbix-VERSIONSTRING.tar.gz

## ------------------------------------------------------------------------
library(Rinbix)
data("scaleFreeNetwork")
datasetObj <- createDiffCoexpMatrixNoME(M = 10, N = 40, meanExpression = 7, 
                                        A = scaleFreeNetwork,
                                        randSdNoise = 0.5, sdNoise = 0.05,
                                        sampleIndicesInteraction = seq(1:3))
dataset <- as.data.frame(datasetObj$regressionData)
colnames(dataset) <- colnames(datasetObj$regressionData)
rownames(dataset) <- rownames(datasetObj$regressionData)
cat("Simulated genes:", colnames(dataset)[1:3], "\n")

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(dataset)

## ------------------------------------------------------------------------
dcgain_snprank_results <- rankDcgainSnprank(dataset)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(dcgain_snprank_results)

## ------------------------------------------------------------------------
regain_snprank_results <- rankRegainSnprank(dataset)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(regain_snprank_results)

