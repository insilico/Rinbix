# -----------------------------------------------------------------------------
# inbixFS.R - Bill White - 10/10/15
#
# Rinbix package feature selection algorithms.

# -----------------------------------------------------------------------------
#' Implements the GeneRank algorithm.
#'
#' http://www.biodatamining.org/content/8/1/2
#' There are three implementations based on the number of genes.
#' 
#' \code{geneRank} 
#' 
#' @param X Matrix n rows of genes by m columns of subjects.
#' @param a Numeric probability?.
#' @param eps Numeric tolerance.
#' @param maxit Numeric maximum nuber of iterations.
#' @return ranked list of genes.
#' @examples
#' data(testdata100ME4)
#' geneRankResults <- rankGeneRank(testdata100ME4)
#' @export
geneRank <- function(X, a=0.9, eps=0.0001, maxit=10, verbose=FALSE) {
  n <- nrow(X)
  if(n <= 5000) {
    if(verbose) cat("Calling geneRank with built-in R function eigen()\n")
    leftMaxEigen1(X, a)
  } else {
    if(n <= 10000) {
      if(verbose) cat("Calling geneRank with power method\n")
      leftMaxEigen2(X, a, eps, maxit)
    } else {
      if(verbose) cat("Calling geneRank with power method, no R2 matrix storage\n")
      leftMaxEigen3(X, a, eps, maxit)
    }
  }
}

# -----------------------------------------------------------------------------
# From the paper Appendix - Computation of the left maximum eigenvector
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' 1. Using R's built in eigen() function
#'
#' NOTE: I had to adapt this code for matrix operations from the copy-and-pasted
#' version from the PDF. Looks like it was orginally Matlab "converted" to R syntax?
#'
#' \code{leftMaxEigen2} 
#'
#' @param X Matrix n rows of genes by m columns of subjects.
#' @param a Numeric probability?.
#' @return left maximum eigenvector.
#' @keywords internal
leftMaxEigen1 <- function(X, a=0.9) {
  n <- nrow(X)
  R2 <- cor(t(X))^2
  sum.row <- rowSums(R2)
  R2.star <- 1/sum.row * R2
  H <- (1-a)/n + a*R2.star
  p <- abs(eigen(t(H))$vectors[,1])
  return(p)
}

# -----------------------------------------------------------------------------
#' 2. The power method for moderate 5000 < n < 10000.
#'
#' NOTE: I had to adapt this code for matrix operations from the copy-and-pasted
#' version from the PDF. Looks like it was orginally Matlab "converted" to R syntax?
#'
#' \code{leftMaxEigen2} 
#' 
#' @param X Matrix n rows of genes by m columns of subjects.
#' @param a Numeric probability?.
#' @param eps Numeric epsilon tolerance stopping criteria.
#' @param maxit Numeric maximum number of iterations.
#' @return left maximum eigenvector.
#' @keywords internal
leftMaxEigen2 <- function(X, a=0.9, eps=0.0001, maxit=10) {
  m <- ncol(X)
  n <- nrow(X)
  R2 <- cor(t(X))^2
  sum.row <- rowSums(R2)
  R2.star <- 1/sum.row * R2
  H <- (1-a)/n + a*R2.star
  tH=t(H)
  p=sum.row/sqrt(sum(sum.row^2))
  # make iteratively better
  for(it in 1:maxit) {
    p.new=tH%*%p
    p.new=p.new/sqrt(sum(p.new^2))
    adiff=max(abs(p-p.new))
    if(adiff<eps) break
    p=p.new[,1]
  }
  return(p.new[,1])
}

# -----------------------------------------------------------------------------
#' 3. for large n > 10000.
#'
#' NOTE: This code worked directly from cut-and-paste from the PDF file.
#'
#' \code{leftMaxEigen3} 
#' 
#' @param X Matrix n rows of genes by m columns of subjects.
#' @param a Numeric probability?.
#' @param eps Numeric epsilon tolerance stopping criteria.
#' @param maxit Numeric maximum number of iterations.
#' @return left maximum eigenvector.
#' @keywords internal
leftMaxEigen3 <- function(X, a=0.9, eps=0.0001, maxit=10) {
  m <- ncol(X)
  n <- nrow(X)
  # compute normalized gene expression matrix
  x.bar=rowMeans(X)
  Xsub.mean=X-x.bar
  sdX=sqrt(rowSums(Xsub.mean^2))
  Z=(1/sdX)*Xsub.mean
  # compute Rstar^2'p without computing Rstar^2
  sumR2=rep(0,n)
  for(i in 1:m) {
    for(j in 1:m)
    {
      qij=sum(Z[,i]*Z[,j])
      sumR2=sumR2+qij*Z[,i]*Z[,j]
    }
  }
  p=sumR2/sqrt(sum(sumR2^2))
  # make iteratively better
  for(it in 1:maxit) {
    tR2p.fast=rep(0,n)
    for(i in 1:m)
      for(j in 1:m)
      {
        hij=sum(Z[,i]*Z[,j]*p/sumR2)
        tR2p.fast=tR2p.fast+hij*Z[,i]*Z[,j]
      }
    p.new=(1-a)/n*sum(p)+a*tR2p.fast
    p.new=p.new/sqrt(sum(p.new^2))
    adiff=max(abs(p-p.new))
    if(adiff<eps) break
    p=p.new
  }
  return(p.new)
}

# -----------------------------------------------------------------------------
#' Relative recurrency variable importance metric (r2VIM).
#' There is a new version. This is based on the earlier PSB 2015 version.
#' 
#' \code{r2VIMorig} 
#' 
#' @param predictors Matrix independent variables.
#' @param response Vector response vector, case-control.
#' @param numRfRuns Numeric of randomForest runs.
#' @param thresholdVIM Numeric threshold importance.
#' @param thresholdMedianRIS Numeric threshold importance RIS score.
#' @param thresholdProb Numeric threshold probability seen in top 10.
#' @param verbose Flag write verbose messages to console/stdout.
#' @return list with: run statistics, votes matrix, network matrix, 
#' importance distribution, importance distribution RIS.
#' @examples
#' data(testdata100ME4)
#' predictors <- as.matrix(testdata100ME4[, -ncol(testdata100ME4)])
#' response <- factor(testdata100ME4[, ncol(testdata100ME4)])
#' r2vimResults <- r2VIMorig(predictors=predictors, response=response, verbose=TRUE)
#' @export
r2VIMorig <- function(predictors=NULL, 
                      response=NULL, 
                      numRfRuns=10,
                      thresholdVIM=0,
                      thresholdMedianRIS=2,
                      thresholdProb=0.2,
                      verbose=FALSE) {
  require(randomForest)
  if(is.null(predictors ) || is.null(response)) {
    stop("r2VIMorig: predictors and response are required parameters")
  }
  numGenes <- ncol(predictors)
  geneIDs <- colnames(predictors)
  # -----------------------------------------------------------------------------
  # compute the variable importances using random forests numRfRuns times,
  # accumulating the importance vectors into a "Dist"ribution data frames
  importanceDistribution <- NULL
  importanceDistributionRIS <- NULL
  # variable votes
  votes <- as.data.frame(matrix(nrow=1, ncol=numGenes, data=c(0)))
  colnames(votes) <- geneIDs
  # GAIN matrix for interaction votes
  networkMatrix <- matrix(ncol=numGenes, nrow=numGenes, data=c(0))
  colnames(networkMatrix) <- geneIDs
  rownames(networkMatrix) <- geneIDs
  if(verbose) cat("Running", numRfRuns, "randomForest models\n")
  runStats <- NULL
  for(testIdx in 1:numRfRuns) {
    # determine variable importance with random forest
    rfResult <- randomForest(predictors, response, importance=TRUE)
    # accumulate importances
    thisRunImportances <- rfResult$importance[, 3]
    estimatedVariance <- abs(min(thisRunImportances))
    thisRunImportancesRIS <- thisRunImportances / estimatedVariance
    importanceDistribution <- rbind(importanceDistribution, thisRunImportances)
    importanceDistributionRIS <- rbind(importanceDistributionRIS, thisRunImportancesRIS)
    # count top ten
    top10Genes <- names(sort(thisRunImportances, decreasing=T)[1:10])
    votes[1, top10Genes] <- votes[1, top10Genes] + 1
    # count all pairs in top 10
    for(pair1idx in 1:9) {
      for(pair2idx in (pair1idx + 1):10) {
        gene1 <- top10Genes[pair1idx]
        gene2 <- top10Genes[pair2idx]
        # add 1/45 to each pair? 45=choose(10, 2)
        networkMatrix[gene1, gene2] <- networkMatrix[gene1, gene2] + 1
        networkMatrix[gene2, gene1] <- networkMatrix[gene1, gene2] + 1
      }
    }
    # classification
    confusionMatrix <- rfResult$confusion
    TN <- confusionMatrix[1, 1]
    FN <- confusionMatrix[2, 1]
    FP <- confusionMatrix[1, 2]
    TP <- confusionMatrix[2, 2]
    SENS <- TP / (TP + FN) # TPR/recall/hit rate/sensitivity
    SPEC <- TN / (TN + FP) # TNR/specificity/SPC
    runStats <- rbind(runStats, data.frame(run=testIdx, 
                                           confusion=confusionMatrix,
                                           sensitivity=SENS,
                                           specificity=SPEC))
    if(verbose) cat(testIdx, " / ", numRfRuns, ": SENS/SPEC: ", round(SENS, 6), 
                    " / ", round(SPEC, 6), "\n", sep="")
  }
  diag(networkMatrix) <- as.numeric(votes[1, ])
 
  list(run.stats=runStats,
       votes=votes,
       net.matrix=networkMatrix, 
       importance.dist=importanceDistribution, 
       importance.dist.ris=importanceDistributionRIS)
}

# ----------------------------------------------------------------------------
#' Rank by dcGAIN + SNPrank.
#' 
#' \code{rankDcgainSnprank} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @param gamma Numeric gamma value for SNPrank
#' @param saveMatrix Flag to save the dcGAIN matrix to text file.
#' @return data frame with SNPrank results: gene, SNPrank, degree.
#' @examples
#' data(testdata100ME4)
#' dcgainResults <- rankDcgainSnprank(testdata100ME4)
#' @export
rankDcgainSnprank <- function(regressionData, gamma=0.85, saveMatrix=FALSE) {
  numPredictors <- ncol(regressionData) - 1
  
  # make predictors (genes from SAM's POV) into rows
  predictors <- t(as.matrix(regressionData[, 1:numPredictors]))
  # phenotype is response
  response <- regressionData[,ncol(regressionData)]
  
  # run Rinbix version of dcGAIN
  #dcResult <- dcgain(ds$regressionData) <- BUG 5/13/14
  dcResult <- dcgainInbix(regressionData)
  #dcResult <- dcgain(regressionData)
  if(saveMatrix) {
    write.table(dcResult$scores, file="test.dcgain", sep="\t", 
      quote=FALSE, row.names=FALSE, col.names=TRUE)
  }

  # SNPrank
  scores <- dcResult$scores
  pvals <- dcResult$pvals
  # threshold based on p-value
  #scores[pvals > 0.05] <- 0
  #cat("p-value threshold sets", length(scores[pvals > 0.05]), "values to zero\n")
  snprankResults <- snprankInbix(scores, gamma=gamma)
}

# ----------------------------------------------------------------------------
#' Rank by GeneRank algorithm.
#' 
#' \code{rankGeneRank} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @return data frame with gene and gene rank score, ordered by score.
#' data(testdata100ME4)
#' rankGeneRankResults <- rankGeneRank(testdata100ME4)
#' @export
rankGeneRank <- function(regressionData) {
  numPredictors <- nrow(regressionData)
  predictors <- as.matrix(regressionData[, 1:numPredictors])
  varNames <- colnames(predictors)
  ranks <- geneRank(t(predictors))
  ranks_df <- data.frame(gene=varNames, GRscore=ranks)
  sorted_ranks <- ranks_df[order(ranks_df$GRscore, decreasing=TRUE),]
}

# -----------------------------------------------------------------------------
#' Rank by glmnet.
#' 
#' \code{rankGlmnet} 
#' 
#' @param regressionData Data frame of predictors and final class column,
#' @param verbose Flag send verbose messages to stdout.
#' @return data frame with gene indices, gene names and coefficients.
#' @examples
#' data(testdata100ME4)
#' rankGlmnetResults <- rankGlmnet(testdata100ME4)
#' @export
rankGlmnet <- function(regressionData, verbose=FALSE) {
  require(glmnet)
  predictors <- regressionData[, -ncol(regressionData)]
  predictor_names <- colnames(predictors)
  response <- factor(regressionData[, ncol(regressionData)], levels=c(0, 1))
  alpha <- 0.5
  glmnet_fit_cv <- cv.glmnet(as.matrix(predictors), response, alpha=alpha, family="binomial")
  lambda_min <- glmnet_fit_cv$lambda.min
  if(verbose) cat("Minimum lambda [", lambda_min, "]\n")
  res <- predict(glmnet_fit_cv, s=lambda_min, type="coef")
  # nonzero indices minus intercept
  nzi <- (res@i)[-1]
  if(verbose) cat("glmnet non-zero indices:", nzi, "\n")
  nzv <- abs(res[nzi+1])
  if(verbose) cat("glmnet non-zero values:", nzv, "\n")
  nznames <- predictor_names[nzi]
  if(verbose) cat("glmnet non-zero names:", nznames, "\n")
  
  resultDF <- data.frame(indices=as.integer(nzi), 
             names=as.character(nznames), 
             values=as.numeric(nzv),
             stringsAsFactors=FALSE)
  resultDF[order(resultDF$values, decreasing=TRUE), ]
}

# ----------------------------------------------------------------------------
#' Rank by Interact algorithm.
#' 
#' \code{rankInteract} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @return data frame with gene and gene and q-value, ordered by q-value.
#' @examples
#' data(testdata100ME4)
#' rankInteractResults <- rankInteract(testdata100ME4)
#' @export
rankInteract <- function(regressionData) {
  require(Interact)
  numPredictors <- ncol(regressionData) - 1
  numPerms <- 100
  topFDR <- choose(numPredictors, 2)
  # make predictors (genes from SAM's POV) into rows
  predictors <- t(as.matrix(regressionData[, 1:numPredictors]))
  # phenotype is response
  response <- regressionData[,ncol(regressionData)]
  # rank by Interact FDR
  fit <- interact(predictors, response, numPerm=numPerms, numFDR=topFDR)
  #print(fit)
  #png(paste(outPrefix, ".png", sep=""), width=1024, height=768)
  #plot(fit)
  #dev.off()
  #write.table(fit$interaction.ordered, file=paste(outPrefix, ".interact.ranks", sep=""), 
  #            col.names=T, row.names=F, quote=F, sep="\t")
  # feature1 feature2 qval
  fit$interaction.ordered
}

# ----------------------------------------------------------------------------
#' Rank by Lasso.
#' 
#' \code{rankLasso} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @return data frame with gene and non-zero coefficients, ordered by coefficient.
#' @examples
#' data(testdata100ME4)
#' rankLassoResults <- rankLasso(testdata100ME4)
#' @export
rankLasso <- function(regressionData) {
  require(glmnet)
  numPredictors <- ncol(regressionData) - 1
  predictors <- as.matrix(regressionData[, 1:numPredictors])
  varNames <- colnames(predictors)
  # scale the predictors! lesson learned 5/2/14
  scaledPredictors <- base::scale(predictors)
  response <- factor(regressionData$Class, levels=c(0, 1))
  cv <- cv.glmnet(scaledPredictors, response, alpha=1, nfolds=5, family="binomial")
  l <- cv$lambda.min
  alpha=0.5
  # fit the model
  fits <- glmnet(scaledPredictors, response, family="binomial", alpha=alpha, nlambda=100)
  res <- predict(fits, s=l, type="coefficients")
  # nonzero indices minus intercept
  nzi <- (res@i)[-1]
  #cat("nzi:", nzi, "\n")
  nzv <- res[nzi+1]
  #cat("nzv:", nzv, "\n")
  returnDF <- data.frame(gene=varNames[nzi], coef=as.vector(nzv))
  returnDF[order(returnDF$coef, decreasing=TRUE), ]
}

# ----------------------------------------------------------------------------
#' Rank by Lasso interact (glinternet).
#' 
#' \code{rankLassoInteract} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @return data frame with gene and non-zero coefficients, ordered by coefficient.
#' @examples
#' data(testdata100ME4)
#' rankLassoInteractResults <- rankLassoInteract(testdata100ME4)
#' @export
rankLassoInteract <- function(regressionData) {
  require(glinternet)
  numPredictors <- ncol(regressionData) - 1
  predictors <- as.matrix(regressionData[, 1:numPredictors])
  varNames <- colnames(predictors)
  response <- regressionData[,ncol(regressionData)]
  cv <- glinternet.cv(predictors, response, family="binomial", numLevels=1)
  l <- cv$lambdaHat
  # fit the model
  fits <- glinternet(predictors, response, family="binomial", numLevels=1, lambda=l)
  res <- predict(fits, s=l, type="coefficients")
  # nonzero indices minus intercept
  nzi <- (res@i)[-1]
  #cat("nzi:", nzi, "\n")
  nzv <- res[nzi+1]
  #cat("nzv:", nzv, "\n")
  returnDF <- data.frame(gene=varNames[nzi], coef=as.vector(nzv))
  returnDF[order(returnDF$score, decreasing=TRUE), ]
}

# ----------------------------------------------------------------------------
#' Rank by limma - linear methods for microarrays.
#' 
#' \code{rankLimma} 
#' 
#' @param regressionData Data frame of predictors and final class column,
#' @return data frame with gene and p-value, ordered by p-value,
#' @examples
#' data(testdata100ME4)
#' rankLimmaResults <- rankLimma(testdata100ME4)
#' @export
rankLimma <- function(regressionData) {
  require(limma)
  # set up the data
  numPredictors <- ncol(regressionData) - 1
  # make predictors (genes from Limma's POV) into rows
  predictors <- t(as.matrix(regressionData[, 1:numPredictors]))
  # phenotype is response
  response <- regressionData[,ncol(regressionData)]
  design <- model.matrix(~response)
  # proceed with analysis
  fit <- lmFit(predictors, design)
  fit2 <- eBayes(fit)
  tT <- topTable(fit2, coef=2, adjust="fdr", sort.by="p", n=Inf)
  returnDF <- data.frame(gene=rownames(tT), score=tT$P.Value)
  returnDF[order(returnDF$score), ]
}

# ----------------------------------------------------------------------------
#' Rank by random forests.
#' 
#' \code{rankRandomForest} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @return data frame with gene and importance score, ordered by importance.
#' @examples
#' data(testdata100ME4)
#' rankRandomForestResults <- rankRandomForest(testdata100ME4)
#' @export
rankRandomForest <- function(regressionData) {
  require(randomForest)
  numPredictors <- ncol(regressionData) - 1
  # make predictors
  predictors <- as.matrix(regressionData[, 1:numPredictors])
  varNames <- colnames(predictors)
  #predictors <- t(predictors)
  # phenotype is response
  response <- as.factor(regressionData[,ncol(regressionData)])
  rfResult <- randomForest(predictors, response, importance=TRUE)
  returnDF <- data.frame(gene=rownames(rfResult$importance), score=rfResult$importance[, 3])
  returnDF[order(returnDF$score, decreasing=TRUE), ]
}

# ----------------------------------------------------------------------------
#' Rank by reGAIN + SNPrank.
#' 
#' \code{rankRegainSnprank} 
#' 
#' @param regressionData Data frame of predictors and final class column,
#' @return data frame with SNPrank results: gene, SNPrank, degree,
#' @examples
#' data(testdata100ME4)
#' rankRegainSnprankResults <- rankRegainSnprank(testdata100ME4)
#' @export
rankRegainSnprank <- function(regressionData, gamma=0.85, saveMatrix=FALSE,  pThresh=1) {
  numPredictors <- ncol(regressionData) - 1
  # run inbix C++ version of reGAIN
  rgResult <- regainInbix(regressionData, stdBetas=FALSE, absBetas=TRUE, pThreshold=pThresh)
  #rgResult <- regainParallel(as.data.frame(regressionData), stdBetas=FALSE, absBetas=TRUE)
  if(saveMatrix) {
    write.table(rgResult$reGAIN, file="test.regain", sep="\t",
      quote=FALSE, row.names=FALSE, col.names=TRUE)  
  }
  # if(length(rgResult$warningsText) > 0) {
  #   print(rgResult$warningsText)
  # }
  # if(length(rgResult$failuresText) > 0) {
  #   print(rgResult$failuresText)
  # }
  # SNPrank
  snprankResults <- snprankInbix(rgResult$reGAIN, gamma=gamma)
}

# ----------------------------------------------------------------------------
#' Rank by ReliefSeq.
#' 
#' \code{rankReliefSeq} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @param outPrefix String output file prefix for temporary files.
#' @param k Numeric k-nearest neighbors in the Relief-F algorithm.
#' @return data frame with ReliefSeq results: gene, score.
#' @examples
#' data(testdata100ME4)
#' rankReliefSeqResults <- rankReliefSeq(testdata100ME4)
#' @export
rankReliefSeq <- function(regressionData, outPrefix="Rinbix", k=10) {
  if(Sys.which("reliefseq") == "") {
    stop("reliefseq is not in the PATH")
  }
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, "Rinbix")
  # run reliefseq command
  relieffCmd <- paste("reliefseq -n Rinbix.num -a Rinbix.pheno -k", k , "-o", outPrefix)
  #cat("Running relieff command:", relieffCmd, "\n")
  relieffStdout <- system(relieffCmd, intern=TRUE)
  relieffRankings <- read.table("Rinbix.reliefseq", head=FALSE, sep="\t")
  file.remove(c("Rinbix.num", "Rinbix.pheno", "Rinbix.reliefseq"))
  data.frame(gene=relieffRankings[, 2], score=relieffRankings[, 1] )
}

# ----------------------------------------------------------------------------
#' Rank by SAM - sequential analysis of microarrays.
#' 
#' \code{rankSam} 
#' 
#' @param regressionData Data frame of predictors and final class column.
#' @return data frame with gene and p-value, ordered by p-value.
#' @examples
#' data(testdata100ME4)
#' rankSamResults <- rankSam(testdata100ME4)
#' @export
rankSam <- function(regressionData) {
  require(samr)
  numPredictors <- ncol(regressionData) - 1
  # make predictors (genes from SAM's POV) into rows
  predictors <- t(as.matrix(regressionData[, 1:numPredictors]))
  # phenotype is response
  response <- regressionData[,ncol(regressionData)] + 1
  capture.output(samFit <- SAM(x=predictors, y=response, resp.type="Two class unpaired", 
                 genenames=rownames(predictors)))
  pVals <- samr.pvalues.from.perms(samFit$samr.obj$tt, samFit$samr.obj$ttstar)
  returnDF <- data.frame(gene=colnames(predictors), score=pVals)
  returnDF[order(returnDF$score), ]
}

# -----------------------------------------------------------------------------
#' Rank genes by SNPrank algorithm.
#' 
#' \code{snprank}
#' 
#' @param G Matrix genetic association interaction network.
#' @param gamma Numeric weighting, interactions (closer to 1) versus main effects (closer to 0).
#' @return sortedTable Data frame with gene, SNPrank, diagonal and degree columns.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' rinbixSnpranksDF <- snprank(rinbixRegain)
#' @export
snprank <- function(G, gamma=0.85) {
  n <- nrow(G);
  geneNames <- colnames(G)
  Gdiag <- diag(G)
  Gtrace <- sum(Gdiag)
  colsum <- colSums(G);
  diag(G) <- 0
  Gtrace <- Gtrace * n;
  colsumG <- colSums(G)
  rowSumG <- rowSums(G)
  rowsum_denom <- matrix(0, n, 1);
  for(i in 1:n) {
    localSum <- 0;
    for(j in 1:n) {
      factor <- ifelse(G[i, j] != 0, 1, 0);
      localSum <- localSum + (factor * colsumG[j]);
    }
    rowsum_denom[i] <- localSum;
  }
  gamma_vec <- rep(gamma, n);
  gamma_matrix <- matrix(nrow=n, ncol=n, data=rep(gamma_vec, n))
  if(Gtrace) {
    b <- ((1.0 - gamma_vec) / n) + (Gdiag / Gtrace)  
  } else {
    b <- ((1.0 - gamma_vec) / n)
  }
  D <- matrix(nrow=n, ncol=n, data=c(0))
  diag(D) <- 1 / colsumG
  I <- diag(n)
  temp <- I - gamma_matrix * G %*% D
  r <- solve(temp, b)
  snpranks <- r / sum(r)
  saveTable <- data.frame(gene=geneNames, snprank=snpranks)
  sortedTable <- saveTable[order(saveTable$snprank, decreasing=TRUE),]  
}

# ----------------------------------------------------------------------------
#' Rank by t-test.
#' 
#' \code{rankReliefSeq} 
#' 
#' @param regressionData Data frame of predictors and final class column,
#' @return data frame with t-test results: gene, p-value.
#' @examples
#' data(testdata100ME4)
#' rankTTestResults <- rankTTest(testdata100ME4)
#' @export
# ----------------------------------------------------------------------------
rankTTest <- function(regressionData) {
  numPredictors <- ncol(regressionData) - 1
  # make predictors
  predictors <- as.matrix(regressionData[, 1:numPredictors])
  varNames <- colnames(predictors)
  predictors <- t(predictors)
  # phenotype is response
  response <- regressionData[,ncol(regressionData)]
  stats_matrix <- NULL
  for(i in 1:nrow(predictors)) {
    g1_data <- predictors[i, response == 0]
    g2_data <- predictors[i, response == 1]
    t_result <- t.test(g1_data, g2_data)
    stats_matrix <- rbind(stats_matrix, t_result$p.value)
  }
  rownames(stats_matrix) <- varNames
  stats_matrix[order(stats_matrix[, 1], decreasing=FALSE), ]
}

# -----------------------------------------------------------------------------
#' Rank genes by univariate regression.
#' 
#' \code{rankUnivariateRegression}
#' 
#' @param regressionData Data frame with gene in columns and samples in rows;
#' the last column should be labeled 'Class' and be 0 or 1 values.
#' @return Data frame with gene, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns.
#' @examples
#' data(testdata100ME4)
#' rankUnivariateRegressionResults <- rankUnivariateRegression(testdata100ME4)
#' @export
rankUnivariateRegression <- function(regressionData) {
  # calculate the logistic regression coefficient for each variable in 
  # regressionData versus the "Class" column (last column)
  colNames <- colnames(regressionData)[1:(ncol(regressionData)-1)]
  numVars <- length(colNames)
  results <- NULL
  regressionData$Class <- factor(regressionData$Class)
  for(i in 1:numVars) {
    glmFormula <- paste("Class ~ ", colNames[i], sep="")
    #print(glmFormula)
    interactionModel <- glm(as.formula(glmFormula), family="binomial", data=regressionData)
    #print(summary(interactionModel))
    glmConverged <- interactionModel$converged
    intCoef <- summary(interactionModel)$coefficients[2, "Estimate"]
    intStdErr <- summary(interactionModel)$coefficients[2, "Std. Error"]
    intPval <- summary(interactionModel)$coefficients[2, "Pr(>|z|)"]
    results <- rbind(results, 
                     cbind(colNames[i],
                           glmConverged, 
                           round(intCoef, 6), 
                           round(intPval, 6), 
                           round(intStdErr, 6), 
                           round(intCoef / intStdErr, 6)))
  }
  colnames(results) <- c("Variable", "Converged", "Beta", "p", "se", "std beta")
  sortedResults <- results[order(results[, 4]), ]
}
