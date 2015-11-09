## ---- eval=FALSE---------------------------------------------------------
#  install.packages("path_to_file/Rinbix-1.0.tar.gz", repos=NULL, type="source")

## ---- eval=FALSE---------------------------------------------------------
#  R CMD INSTALL Rinbix-1.0.tar.gz

## ---- echo=FALSE---------------------------------------------------------
library(Rinbix)

## ------------------------------------------------------------------------
data("scaleFreeNetwork")
dataset <- createDiffCoexpMatrixNoME(M=10, N=40, meanExpression=7, 
                                     A=scaleFreeNetwork,
                                     randSdNoise=0.5, sdNoise=0.05,
                                     sampleIndicesInteraction=seq(1:3))
cat("Simulated genes:", colnames(dataset$regressionData)[1:3], "\n")

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(dataset$regressionData))

## ------------------------------------------------------------------------
dcgain_snprank_results <- rankDcgainSnprank(dataset$regressionData)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(dcgain_snprank_results)

## ------------------------------------------------------------------------
regain_snprank_results <- rankRegainSnprank(dataset$regressionData)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(regain_snprank_results)

