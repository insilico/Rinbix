# ------------------------------------------------------------------------------------
# inbixRnaseq.R - Bill White - 10/12/15
#
# Rinbix package RNA-Seq functions.

# ------------------------------------------------------------------------------------
#' Rank using DESeq.
#' 
#' \code{rankDeseq} 
#' 
#' @param rnaExpr matrix of rna transcripts: gene rows by subject columns
#' @param groups vector of case-control groups 0/1
#' @return differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' rDeseq <- rankDeseq(X, y)
#' @export
rankDeseq <- function(rnaExpr, groups) {
	deData <- newCountDataSet(rnaExpr, as.factor(t(groups)))
	deData <- estimateSizeFactors(deData)
	deData <- estimateDispersions(deData)
	results <- nbinomTest(deData, "0", "1")
	resultsSorted <- results[order(results$pval), ]
	resultsSave <- cbind(resultsSorted$pval, resultsSorted[, 1])
	colnames(resultsSave) <- c("pval", "gene")
	resultsSave
}

# ------------------------------------------------------------------------------------
#' Rank using DESeq2.
#' 
#' \code{rankDeseq2} 
#' 
#' @param rnaExpr matrix of rna transcripts: gene rows by subject columns
#' @param groups vector of case-control groups 0/1
#' @return differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' rDeseq2 <- rankDeseq2(X, y)
#' @export
rankDeseq2 <- function(rnaExpr, groups) {
	colData <- data.frame(groups=groups)
	deData <- DESeqDataSetFromMatrix(countData=rnaExpr, colData=colData, design=~groups)
	deData <- DESeq(deData, betaPrior=F)
	results <- results(deData)
	resultsSorted <- results[order(results$pvalue), ]
}

# ------------------------------------------------------------------------------------
#' Rank using edgeR.
#' 
#' \code{rankEdgeR} 
#' 
#' @param rnaExpr matrix of rna transcripts: gene rows by subject columns
#' @param groups vector of case-control groups 0/1
#' @return differentially expressed genes ordered by p-value
#' @examples
#' data(simrnaseq)
#' X <- t(predictorsTrain)
#' y <- as.factor(responseTrain - 1)
#' rEdger <- rankEdgeR(X, y)
#' @export
rankEdgeR <- function(rnaExpr, groups) {
	y <- DGEList(counts=rnaExpr, group=groups)
	y <- estimateCommonDisp(y)
	d <- exactTest(y)
	results <- d$table
	results[order(results$PValue), ]
	results
}

# ------------------------------------------------------------------------------------
#' Rank using edgeR GLM model.
#' 
#' \code{rankEdgeR}Glm 
#' 
#' @param rnaExpr matrix of rna transcripts: gene rows by subject columns
#' @param groups vector of case-control groups 0/1
#' @param covariates covariates to add to the GLM model
#' @return differentially expressed genes ordered by p-value
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
	d <- DGEList(counts=rnaExpr, group=groups)
	# work with the DEGlist as above
	d <- estimateGLMCommonDisp(d, design)
	d <- estimateGLMTrendedDisp(d, design)
	d <- estimateGLMTagwiseDisp(d, design)   
	# fit the general linear model
	# EdgeR glmfit
	fit <- glmFit(d, design)
	# conduct the likelihood ratio test
	lrt <- glmLRT(fit)
	lrt
}
