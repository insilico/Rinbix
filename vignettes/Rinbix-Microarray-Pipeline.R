## ----libraries, echo=FALSE, message=FALSE, warning=FALSE-----------------
# library(BiocParallel)
# register(DoparParam())
# register(MulticoreParam())
# library(affy)
# library(affyPLM)
library(preprocessCore)
library(leukemiasEset)
library(clusterProfiler)
library(broom)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Rinbix)
library(biomaRt)
knitr::opts_chunk$set(fig.width = 12, fig.height = 8) 

## ----dataset-------------------------------------------------------------
data(leukemiasEset)  # load bone marrow samples
# objects in the data
#str(leukemiasEset)
#slotNames(leukemiasEset) # phenoData, etc
dim(leukemiasEset) # 20172 genes x 60 samples
Biobase::sampleNames(leukemiasEset) # cell file names
allPheno <- Biobase::pData(leukemiasEset)
#head(allPheno)  # look at phenotype info
leukPheno <- allPheno$LeukemiaType  # abbreviated leukemia types
summary(leukPheno)
# ALL AML CLL CML NoL 
# 12  12  12  12  12 
#featureNames(leukemiasEset)[1:5] # first 5 gene ensemble ids
leukExprData <- Biobase::exprs(leukemiasEset) # exprs is an affy function to extract expression data from eset
colnames(leukExprData) <- leukPheno  # add phenotype names to matrix

## ----preprocess----------------------------------------------------------
# quantiles function needs eset to operate on
leukExprData_quantile <- preprocessCore::normalize.quantiles(leukExprData)
boxplot(leukExprData_quantile,range = 0,ylab = "raw intensity", main = "Quantile Normalized")
#leukExprData_quantileLog2 <- log2(exprs(leukExprData_quantile))
leukExprData_quantileLog2 <- log2(leukExprData_quantile)
colnames(leukExprData_quantileLog2) <- leukPheno  # add phenotype names to matrix
boxplot(leukExprData_quantileLog2,range = 0,ylab = "log2 intensity", 
        main = "Quantile Normalized Log2")
# compare the ALL and NoL groups; find genes that are differentially expressed between groups
# use t.test
# 1. create subset of data for the two groups
ALL.NoL.mask <- colnames(leukExprData) == "ALL" | colnames(leukExprData) == "NoL"
ALL.NoL.Data <- leukExprData[,ALL.NoL.mask]

## ----cov-----------------------------------------------------------------
# there are a lot of genes that have very low signal to noise that we can get rid of.
coef.of.vars <- apply(ALL.NoL.Data, 1, function(x) { sd(x) / abs(mean(x)) })
# the smaller the threshold, the higher the experimental effect relative to the measurement precision
sum(coef.of.vars < 0.05)  # 5,378 genes
# filter the data matrix
ALL.NoL.Data.filter <- ALL.NoL.Data[coef.of.vars < 0.05,]
dim(ALL.NoL.Data.filter)
phenos <- as.integer(factor(colnames(ALL.NoL.Data.filter)))

## ----ttest---------------------------------------------------------------
# put it all together. apply to all genes
# i is the data row or gene index
filteredGenes <- rownames(ALL.NoL.Data.filter)
filteredDataset <- cbind(t(ALL.NoL.Data.filter), phenos - 1)
colnames(filteredDataset) <- c(filteredGenes, "Class")
ttest.results <- rankTTest(filteredDataset)
knitr::kable(head(ttest.results$variable), row.names = FALSE, digits = 6)

