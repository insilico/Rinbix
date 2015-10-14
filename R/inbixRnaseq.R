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
#' @export
rankDeseq <- function(rnaExpr, groups) {
	deData <- newCountDataSet(rnaExpr, as.factor(t(groups)))
	deData <- estimateSizeFactors(deData)
	deData <- estimateDispersions(deData)
	results <- nbinomTest(deData, "0", "1")
	resultsSorted <- results[order(results$pval), ]
	resultsSave <- cbind(resultsSorted$pval, resultsSorted[, 1])
	colnames(resultsSave) <- c("pval", "gene")
}

# ------------------------------------------------------------------------------------
#' Rank using DESeq2.
#' 
#' \code{rankDeseq2} 
#' 
#' @param rnaExpr matrix of rna transcripts: gene rows by subject columns
#' @param groups vector of case-control groups 0/1
#' @return differentially expressed genes ordered by p-value
#' @export
rankDeseq2 <- function(rnaExpr, groups) {
	countData <- rnaExpr
	colGroups <- ifelse(groups > 0, 1, 0)
	colData <- data.frame(groups=factor(colGroups, levels=c(0, 1)))
	dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~groups)
	dds <- DESeq(dds)
	res <- results(dds)
	resultsSorted <- res[order(res$pvalue), ]
}

# ------------------------------------------------------------------------------------
#' Rank using edgeR.
#' 
#' \code{rankEdgeR} 
#' 
#' @param rnaExpr matrix of rna transcripts: gene rows by subject columns
#' @param groups vector of case-control groups 0/1
#' @return differentially expressed genes ordered by p-value
#' @export
rankEdgeR <- function(rnaExpr, groups) {
	# Call edgeR with a design matrix
	# design <- model.matrix(~Class, data=edgerData)
	# y <- estimateGLMCommonDisp(t(edgerData), design)
	# y <- estimateGLMTrendedDisp(edgerData, design)
	# y <- estimateGLMTagwiseDisp(edgerData, design)
	# fit <- glmFit(y, design)
	# lrt <- glmLRT(fit, coef=2)
	# resultsTable <- topTags(lrt)
	y <- DGEList(counts=rnaExpr, group=groups)
	y <- estimateCommonDisp(y)
	dgeResult <- exactTest(y)
	resultsTable <- dgeResult$table
	resultsTable[order(resultsTable$PValue), ]
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
#' @export
rankEdgeRGlm <- function(rnaExpr, groups, covariates) {
	design <- model.matrix(~groups + covariates)
	# make the DGElist object as above
	d <- DGEList(counts=rnaExpr, group=phenos)
	# work with the DEGlist as above
	d <- estimateGLMCommonDisp(d, design)
	d <- estimateGLMTrendedDisp(d, design)
	d <- estimateGLMTagwiseDisp(d, design)   
	# fit the general linear model
	# EdgeR glmfit
	fit <- glmFit(d, design)
	# conduct the likelihood ratio test
	lrt <- glmLRT(fit)
}
