# inbixRnaseq.R - Bill White - 10/12/15
#
# Rinbix package RNA-Seq functions.

#' Rank using DESeq.
#' 
#' \code{rankDeseq} 
#' 
#' @family RNA-Seq functions
#' @family feature selection functions
#' @param rnaExpr \code{matrix} of rna transcripts: gene rows by subject columns
#' @param groups \code{vector} of case-control groups 0/1
#' @return \code{data.frame} differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' rDeseq <- rankDeseq(X, y)
#' @export
rankDeseq <- function(rnaExpr, groups) {
	deData <- DESeq::newCountDataSet(rnaExpr, as.factor(t(groups)))
	deData <- DESeq::estimateSizeFactors(deData)
	deData <- DESeq::estimateDispersions(deData)
	results <- DESeq::nbinomTest(deData, "0", "1")
	resultsSorted <- results[order(results$pval), ]
	resultsSave <- cbind(resultsSorted$pval, resultsSorted[, 1])
	colnames(resultsSave) <- c("pval", "gene")
	resultsSave
}

#' Rank using DESeq2.
#' 
#' \code{rankDeseq2} 
#' 
#' @family RNA-Seq functions
#' @family feature selection functions
#' @param rnaExpr \code{matrix} of rna transcripts: gene rows by subject columns
#' @param groups \code{vector} of case-control groups 0/1
#' @return \code{data.frame} differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' rDeseq2 <- rankDeseq2(X, y)
#' @export
rankDeseq2 <- function(rnaExpr, groups) {
	colData <- data.frame(groups = groups)
	deData <- DESeq2::DESeqDataSetFromMatrix(countData = rnaExpr, colData = colData, design = ~groups)
	deData <- DESeq2::DESeq(deData, betaPrior = FALSE)
	results <- DESeq2::results(deData)
	resultsSorted <- results[order(results$pvalue), ]
}

#' Rank using edgeR.
#' 
#' \code{rankEdgeR} 
#' 
#' @family RNA-Seq functions
#' @family feature selection functions
#' @param rnaExpr \code{matrix} of rna transcripts: gene rows by subject columns
#' @param groups \code{vector} of case-control groups 0/1
#' @return \code{data.frame} differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' rEdger <- rankEdgeR(X, y)
#' @export
rankEdgeR <- function(rnaExpr, groups) {
	y <- edgeR::DGEList(counts = rnaExpr, group = groups)
	y <- edgeR::estimateCommonDisp(y)
	d <- edgeR::exactTest(y)
	results <- d$table
	results[order(results$PValue), ]
	results
}

#' Rank using edgeR GLM model.
#' 
#' \code{rankEdgeR}Glm 
#' 
#' @family RNA-Seq functions
#' @family feature selection functions
#' @param rnaExpr \code{matrix} of rna transcripts: gene rows by subject columns
#' @param groups \code{vector} of case-control groups 0/1
#' @param covariates \code{vector} of covariates to add to the GLM model
#' @return \code{data.frame} differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' covars <- rnorm(ncol(X))
#' rEdgerGlm <- rankEdgeRGlm(X, y, covars)
#' @export
rankEdgeRGlm <- function(rnaExpr, groups, covariates) {
	design <- model.matrix(~groups+covariates)
	# make the DGElist object as above
	d <- edgeR::DGEList(counts = rnaExpr, group = groups)
	# work with the DEGlist as above
	d <- edgeR::estimateGLMCommonDisp(d, design)
	d <- edgeR::estimateGLMTrendedDisp(d, design)
	d <- edgeR::estimateGLMTagwiseDisp(d, design)   
	# fit the general linear model
	# EdgeR glmfit
	fit <- edgeR::glmFit(d, design)
	# conduct the likelihood ratio test
	lrt <- edgeR::glmLRT(fit)
	lrt
}
