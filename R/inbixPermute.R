# ----------------------------------------------------------------------------
# inbixPermute.R - Bill White - 10/10/15
#
# Rinbix package permutation functions.

# ----------------------------------------------------------------------------
#' Permute SNPrank gene scores and determine classification accuracy.
#'
#' \code{permuteSnpranks}
#' 
#' @param ds Data frame data set.
#' @param genesHit Vector gene indices of simulated genes.
#' @param snprankResults Data frame genes and their SNPrank scores.
#' @param M Numeric number of genes.
#' @param N Numeric number of subjects.
#' @param rankNull Flag rank as null.
#' @param qqPlot Flag plot QQ.
#' @param qqPlotFilename String QQ plot filename.
#' @param numPerms Numeric number of permutation.
#' @param threshold Numeric p-value significance threshold.
#' @param method String method to permute.
#' @param gamma Numeric SNPrank gamma.
#' @param pTh Numeric p-value threshold for method.
#' @return data frame of classification statistics.
#' @export
permuteSnpranks <- function(ds, genesHit, snprankResults, M, N, rankNull=FALSE,
  qqPlot=FALSE, qqPlotFilename="foo.png", numPerms=10, threshold=0.05, 
  method="regain", gamma=0.85, pTh=1) {
  ranksByGene <- snprankResults[order(snprankResults$gene), ]
  
  if(qqPlot) {
    png(qqPlotFilename, width=800, height=800)
    rQQ <- qqnorm(ranksByGene$snprank, plot.it=T, main="SNPRank Q-Q Plot", 
      xlab="Theoretical", ylab="SNPrank")
    qqline(ranksByGene$snprank, distribution=qnorm)
    degreeColors <- rep("white", numGenes)
    degreeColors[genesHit] <- "red"
    points(rQQ$x, rQQ$y, col=degreeColors)
    dev.off()
  }

  # -----------------------------------------------------------------------------
  # get rankings on passed data set  
  data1 <- ds$regressionData
  snprankResultsSim <- ranksByGene

  # -----------------------------------------------------------------------------
  # permute the data set numPerms times, recording snpranks
  perms <- NULL
  n1 <- N / 2
  n2 <- N / 2
  for(permIdx in 1:numPerms) {
    #cat(permIdx, " ", sep="") 
    # permute Class column labels
    classLabels <- sample(c(rep(0, n1), rep(1, n2)), N, replace=FALSE)
    newData <- cbind(data1[, 1:M], classLabels)
    colnames(newData) <- colnames(data1)

    # run ranker
    if(method == "regain") {
      sr <- rankRegainSnprank(newData, gamma, pThresh=pTh)
    }
    if(method == "dcgain") {
      sr <- rankDcgainSnprank(newData, gamma)
    }
    if(method == "relieff") {
      sr <- rankReliefF(newData)
    }
    if(method == "randomforest") {
      sr <- rankRandomForest(newData)
    }
    sr <- sr[order(sr$gene), ]  

    # save scores for this permutation
    perms <- rbind(perms, sr$snprank)
  }

  # -----------------------------------------------------------------------------
  #cat("Finding genes that pass the", threshold, "threshold from permutations\n")
  thresholdIdx <- as.integer(numPerms * (1-threshold))
  thresholdVals <- c()
  rankedClassification <- factor(rep(1, M), levels=c(1,2))
  hitClass <- 2
  for(geneIdx in 1:M) {
    thisGeneScores <- as.numeric(perms[, geneIdx])
    #print(thisGeneScores)
    thresholdScore <- sort(thisGeneScores)[thresholdIdx]
    #print(thresholdScore)
    thresholdVals <- c(thresholdVals, thresholdScore)
    if(snprankResultsSim[geneIdx, 2] >= thresholdScore) {
      #cat("Gene index:", geneIdx, snprankResultsSim[geneIdx, 2] , " >= threshold score:", thresholdScore, "\n") 
      rankedClassification[geneIdx] <- hitClass
    }
  }

  # cat("Summary of genes threshold values\n")
  # print(summary(thresholdVals))
  # cat("Variance:", var(thresholdVals), "\n")

  # -----------------------------------------------------------------------------
  # cat("Ranker found [", length(rankedClassification[rankedClassification == hitClass]), 
  #   "] significant genes out of [", M, "]\n")
  simClassification <- factor(rep(1, M), levels=c(1,2))
  if(!rankNull) {
    simClassification[genesHit] <- hitClass
  }

  # cat(length(trueClassification[trueClassification == 1]), "\n")
  # cat(length(trueClassification[trueClassification == 2]), "\n")
  # cat(length(rankedClassification[rankedClassification == 1]), "\n")
  # cat(length(rankedClassification[rankedClassification == 2]), "\n")

  classificationMetrics <- computeFeatureMetrics(rankedClassification, simClassification)
}

# ----------------------------------------------------------------------------
#' Permute SNPrank gene scores and determine classification accuracy using inbix.
#'
#' \code{permuteSnpranksInbix}
#' 
#' @param ds Data frame data set.
#' @param genesHit Vector gene indices of simulated genes.
#' @param snprankResults Data frame genes and their SNPrank scores.
#' @param numPerms Numeric number of permutation.
#' @param threshold Numeric p-value significance threshold.
#' @param method String method to permute.
#' @param pTh Numeric p-value threshold for method.
#' @param rankNull Flag rank as null.
#' @return data frame of classification statistics.
#' @export
permuteSnpranksInbix <- function(ds, genesHit, snprankResults, numPerms=10, 
  threshold=0.05, method="regain", pTh=1, rankNull=FALSE) {

  # call inbix C++ permutation
  geneThresholds <- permuteGainInbix(ds$regressionData, method=method, 
    numPerms=numPerms, threshold=threshold, pThresh=pTh)
  M <- nrow(geneThresholds)

  hitClass <- 2
  rankedClassification <- factor(rep(1, M), levels=c(1,2))

  ranksByGene <- snprankResults[order(snprankResults$gene), ]
  for(geneIdx in 1:M) {
    thisGeneScore <- ranksByGene$snprank[geneIdx]
    thresholdScore <- geneThresholds[geneIdx, 2]
    #cat(thisGeneScore, thresholdScore, "\n")
    if(thisGeneScore >= thresholdScore) {
      rankedClassification[geneIdx] <- hitClass
    }
  }

  # -----------------------------------------------------------------------------
  # cat("Ranker found [", length(rankedClassification[rankedClassification == hitClass]), 
  #   "] significant genes out of [", M, "]\n")
  simClassification <- factor(rep(1, M), levels=c(1,2))
  if(!rankNull) {
    simClassification[genesHit] <- hitClass
  }

  # print(simClassification)
  # cat(length(trueClassification[trueClassification == 1]), "\n")
  # cat(length(trueClassification[trueClassification == 2]), "\n")
  # cat(length(rankedClassification[rankedClassification == 1]), "\n")
  # cat(length(rankedClassification[rankedClassification == 2]), "\n")

  classificationMetrics <- computeFeatureMetrics(rankedClassification, simClassification)
}

# ----------------------------------------------------------------------------
#' Permute SNPrank gene scores and determine classification accuracy - Simple.
#'
#' \code{permuteSnpranksSimple}
#' 
#' @param ds Data frame data set.
#' @param method String method to permute.
#' @param numPerms Numeric number of permutation.
#' @param threshold Numeric p-value significance threshold.
#' @return data frame gene and empirical thresholds.
#' @export
permuteSnpranksSimple <- function(ds, method, numPerms, threshold) {
  ranksByGene <- snprankResults[order(snprankResults$gene), ]
  
  # -----------------------------------------------------------------------------
  # get rankings on passed data set  
  data1 <- ds

  # -----------------------------------------------------------------------------
  # permute the data set numPerms times, recording snpranks
  M <- ncol(ds) - 1
  N <- nrow(ds)
  n1 <- N / 2
  n2 <- N / 2
  perms <- NULL
  for(permIdx in 1:numPerms) {
    cat(permIdx, " ", sep="") 
    # permute Class column labels
    classLabels <- sample(c(rep(0, n1), rep(1, n2)), N, replace=FALSE)
    newData <- cbind(data1[, 1:M], classLabels)
    colnames(newData) <- colnames(data1)

    # run ranker
    if(method == "regain") {
      sr <- rankRegainSnprank(newData, gamma, pThresh=pTh)
    }
    if(method == "dcgain") {
      sr <- rankDcgainSnprank(newData, gamma)
    }
    if(method == "relieff") {
      sr <- rankReliefF(newData, k=5)
    }
    if(method == "randomforest") {
      sr <- rankRandomForest(newData)
    }
    sr <- sr[order(sr$gene), ]  

    # save scores for this permutation
    perms <- rbind(perms, sr$snprank)
  }
  cat("\n")

  # -----------------------------------------------------------------------------
  thresholdIdx <- as.integer(numPerms * (1-threshold))
  thresholdVals <- c()
  for(geneIdx in 1:M) {
    thisGeneScores <- as.numeric(perms[, geneIdx])
    thresholdScore <- sort(thisGeneScores)[thresholdIdx]
    cat(geneIdx, thresholdScore, "\n")
    print(summary(thisGeneScores))
    thresholdVals <- c(thresholdVals, thresholdScore)
  }

  data.frame(gene=ranksByGene$gene, threshold=thresholdVals)
}
