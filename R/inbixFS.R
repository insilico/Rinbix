# inbixFS.R - Bill White - 10/10/15
#
# Rinbix package feature selection/ranking algorithms.
# Refactored ranker algorithms to return common data frame - bcw 7/7/17

#' Rank by dcGAIN + SNPrank.
#'
#' Run dcGAIN on the labelledDataFrame creating a GAIN matrix.
#' Then run SNPrank on the GAIN matrix.
#' 
#' \code{rankDcgainSnprank} 
#' 
#' @family inbix interface functions
#' @keywords models GAIN differential correlation
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @param gammaParam \code{numeric} gamma value for SNPrank
#' @param saveMatrix \code{logical} to save the dcGAIN matrix to text file.
#' @return \code{data.frame} with: variable, score.
#' @examples
#' data(testdata10)
#' dcgainResults <- rankDcgainSnprank(testdata10)
#' @export
rankDcgainSnprank <- function(labelledDataFrame, gammaParam=0.85, saveMatrix=FALSE) {
  # run Rinbix version of dcGAIN
  #dcResult <- dcgain(ds$labelledDataFrame) <- BUG 5/13/14
  dcResult <- dcgainInbix(labelledDataFrame)
  if (saveMatrix) {
    write.table(dcResult$scores, file = "test.dcgain", sep = "\t", 
                quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  # SNPrank
  scores <- dcResult$scores
  # threshold based on p-value
  #scores[pvals > 0.05] <- 0
  #cat("p-value threshold sets", length(scores[pvals > 0.05]), "values to zero\n")
  snprankResults <- snprank(scores, alg.gamma = gammaParam)
  snprankResults
}

#' Rank by GeneRank algorithm.
#' 
#' \code{rankGeneRank} 
#' 
#' @keywords models eigenvectors
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @return \code{data.frame} with gene and gene rank score, ordered by score.
#' data(testdata10)
#' rankGeneRankResults <- rankGeneRank(testdata10)
#' @export
rankGeneRank <- function(labelledDataFrame) {
  numPredictors <- nrow(labelledDataFrame)
  predictors <- as.matrix(labelledDataFrame[, 1:numPredictors])
  geneRank(t(predictors))
}

#' Rank by glmnet.
#' 
#' \code{rankGlmnet} 
#' 
#' @keywords models regression glm reduction
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column,
#' @param verbose \code{logical} send verbose messages to stdout.
#' @return \code{data.frame} with non-zero variable indices, variable names and coefficients.
#' @examples
#' data(testdata10)
#' rankGlmnetResults <- rankGlmnet(testdata10)
#' @export
rankGlmnet <- function(labelledDataFrame, verbose = FALSE) {
  predictors <- labelledDataFrame[, -ncol(labelledDataFrame)]
  predictor_names <- colnames(predictors)
  response <- factor(labelledDataFrame[, ncol(labelledDataFrame)], levels = c(0, 1))
  alpha <- 0.5
  glmnet_fit_cv <- glmnet::cv.glmnet(as.matrix(predictors), response, alpha = alpha, family = "binomial")
  lambda_min <- glmnet_fit_cv$lambda.min
  if (verbose) cat("Minimum lambda [", lambda_min, "]\n")
  res <- predict(glmnet_fit_cv, s = lambda_min, type = "coef")
  # nonzero indices minus intercept
  nzi <- (res@i)[-1]
  if (verbose) cat("glmnet non-zero indices:", nzi, "\n")
  nzv <- abs(res[nzi + 1])
  if (verbose) cat("glmnet non-zero values:", nzv, "\n")
  nznames <- predictor_names[nzi]
  if (verbose) cat("glmnet non-zero names:", nznames, "\n")
  
  resultDF <- data.frame(indices = as.integer(nzi), 
                         names = as.character(nznames), 
                         values = as.numeric(nzv),
                         stringsAsFactors = FALSE)
  resultDF[order(resultDF$values, decreasing = TRUE), ]
  data.frame(variable = resultDF$names, 
             score = resultDF$values, 
             stringsAsFactors = FALSE)
}

#' Rank by Iterative Relief-F.
#' 
#' \code{rankIterativeRelieff} 
#' 
#' @family feature selection functions
#' @keywords relief nearest neighbors interactions reduction
#' @seealso \link{rankRelieff}
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @param percentRemovePerIteration \code{numeric} percent of attributes to remover per iteration.
#' @param targetNumAttributes \code{numeric} target number of attributes.
#' @param verbose \code{logical} to send messages to stdout.
#' @return \code{list} with \code{data.frame} of reduced labelledDataFrame 
#'   with the target number of attributes and a \code{data.frame} of scores.
#' @examples
#' data(testdata10)
#' irelieffResults <- rankIterativeRelieff(testdata10)
#' @export
rankIterativeRelieff <- function(labelledDataFrame, 
                                 percentRemovePerIteration = 10, 
                                 targetNumAttributes = 10,
                                 verbose = FALSE) {
  classCol <- as.integer(labelledDataFrame[, ncol(labelledDataFrame)])
  curData <- labelledDataFrame[, -ncol(labelledDataFrame)]
  curNumAttributes <- ncol(curData)
  iterDone <- FALSE
  iteration <- 0
  scoresDF <- NULL
  while (!iterDone) {
    iteration <- iteration + 1
    # run Relief-F on this subset
    curData$Class <- classCol
    thisRanking <- rankRelieff(curData)
    thisRanking <- thisRanking[order(thisRanking$score), ]
    thisNumToRemove <- as.integer(curNumAttributes * (percentRemovePerIteration / 100))
    if ((thisNumToRemove > 0) && (curNumAttributes - thisNumToRemove) > 0) {
      attributesToRemove <- head(thisRanking, n = thisNumToRemove)
      scoresDF <- rbind(scoresDF, data.frame(variable = attributesToRemove$variable, 
                                             score = attributesToRemove$score))
      if (verbose) cat("Number of attributes:", curNumAttributes, 
                      ", removing:", thisNumToRemove, "\n")
      curData <- curData[, -which(colnames(curData) %in% attributesToRemove$variable)]
      curNumAttributes <- ncol(curData) - 1
    } else {
      scoresDF <- rbind(scoresDF, data.frame(variable = thisRanking$variable, 
                                             score = thisRanking$score,
                                             stringsAsFactors = FALSE))
      iterDone <- TRUE
    }
  }
  scoresDF <- scoresDF[order(scoresDF$score, decreasing = TRUE), ]
  if (verbose) cat("Iterative Relief-F complete in [", iteration, "] iterations\n")
  scoresDF
}

#' Rank by Lasso.
#' 
#' \code{rankLasso} 
#' 
#' @family feature selection functions reduction
#' @keywords glm regression
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @return \code{data.frame} with non-zero variable coefficients, ordered by coefficient.
#' @examples
#' data(testdata10)
#' rankLassoResults <- rankLasso(testdata10)
#' @export
rankLasso <- function(labelledDataFrame) {
  numPredictors <- ncol(labelledDataFrame) - 1
  predictors <- as.matrix(labelledDataFrame[, 1:numPredictors])
  varNames <- colnames(predictors)
  # scale the predictors! lesson learned 5/2/14
  scaledPredictors <- base::scale(predictors)
  response <- factor(labelledDataFrame$Class, levels = c(0, 1))
  cv <- glmnet::cv.glmnet(scaledPredictors, response, alpha = 1, nfolds = 5, 
                          family = "binomial")
  l <- cv$lambda.min
  alpha <- 0.5
  # fit the model
  fits <- glmnet::glmnet(scaledPredictors, response, family = "binomial", 
                         alpha = alpha, nlambda = 100)
  res <- predict(fits, s = l, type = "coefficients")
  # nonzero indices minus intercept
  nzi <- (res@i)[-1]
  #cat("nzi:", nzi, "\n")
  nzv <- res[nzi + 1]
  #cat("nzv:", nzv, "\n")
  returnDF <- data.frame(variable = varNames[nzi], 
                         score = as.vector(nzv),
                         stringsAsFactors = FALSE)
  returnDF[order(returnDF$score, decreasing = TRUE), ]
  returnDF
}

#' Rank by limma - linear methods for microarrays.
#' 
#' \code{rankLimma} 
#' 
#' @keywords models linear microarray
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column,
#' @return \code{data.frame} with variable and p-value, ordered by p-value,
#' @examples
#' data(testdata10)
#' rankLimmaResults <- rankLimma(testdata10)
#' @export
rankLimma <- function(labelledDataFrame) {
  # set up the data
  numPredictors <- ncol(labelledDataFrame) - 1
  # make predictors (variables from Limma's POV) into rows
  predictors <- t(as.matrix(labelledDataFrame[, 1:numPredictors]))
  # phenotype is response
  response <- labelledDataFrame[, ncol(labelledDataFrame)]
  design <- model.matrix(~response)
  # proceed with analysis
  fit <- limma::lmFit(predictors, design)
  fit2 <- limma::eBayes(fit)
  tT <- limma::topTable(fit2, coef = 2, adjust.method = "fdr", sort.by = "p", 
                        number = Inf)
  returnDF <- data.frame(variable = rownames(tT), 
                         score = tT$P.Value,
                         stringsAsFactors = FALSE)
  returnDF[order(returnDF$score), ]
  returnDF
}

#' Rank by random forests.
#' 
#' \code{rankRandomForest} 
#' 
#' @keywords models ensembles decision trees
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @return \code{data.frame} with variable and importance score, ordered by importance.
#' @examples
#' data(testdata10)
#' rankRandomForestResults <- rankRandomForest(testdata10)
#' @export
rankRandomForest <- function(labelledDataFrame) {
  numPredictors <- ncol(labelledDataFrame) - 1
  # make predictors
  predictors <- as.matrix(labelledDataFrame[, 1:numPredictors])
  #varNames <- colnames(predictors)
  #predictors <- t(predictors)
  # phenotype is response
  response <- as.factor(labelledDataFrame[,ncol(labelledDataFrame)])
  rfResult <- randomForest::randomForest(predictors, response, importance = TRUE)
  returnDF <- data.frame(variable = rownames(rfResult$importance), 
                         score = rfResult$importance[, 3],
                         stringsAsFactors = FALSE)
  returnDF[order(returnDF$score, decreasing = TRUE), ]
  returnDF
}

#' Rank by reGAIN + SNPrank.
#' 
#' \code{rankRegainSnprank} 
#' 
#' @keywords models regression GAIN interactions
#' @family feature selection functions
#' @family inbix interface functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column,
#' @param paramGamma \code{numeric} SNPrank gamma, NOTE! built in function gamma!
#' @param saveMatrix \code{logical} to save reGAIN matrix to file.
#' @param pThresh \code{numeric} probability threshold for reGAIN models.
#' @param useAbs \code{logical} take absolute value of reGAIN matrix
#' @param useStdCoef \code{logical} use coefficient statstic (T) or raw coefficient (F)
#' @return \code{data.frame} with SNPrank results: variable, SNPrank, degree,
#' @examples
#' data(testdata10)
#' rankRegainSnprankResults <- rankRegainSnprank(testdata10)
#' @export
rankRegainSnprank <- function(labelledDataFrame, paramGamma = 0.85, 
                              saveMatrix = FALSE, pThresh = 1, 
                              useAbs = TRUE, useStdCoef = TRUE) {
  # run inbix C++ version of reGAIN
  #rgResult <- regainInbix(labelledDataFrame, stdBetas = TRUE, absBetas = TRUE, 
  #                        pThreshold = pThresh)
  rgResult <- regainParallel(labelledDataFrame, 
                             stdBetas = useStdCoef, 
                             absBetas = useAbs,
                             numCores = 2, 
                             verbose = FALSE)
  if (saveMatrix) {
    write.table(rgResult, file = "test.regain", sep = "\t",
      quote = FALSE, row.names = FALSE, col.names = TRUE)  
  }
  # if(length(rgResult$warningsText) > 0) {
  #   print(rgResult$warningsText)
  # }
  # if(length(rgResult$failuresText) > 0) {
  #   print(rgResult$failuresText)
  # }
  # SNPrank
  snprankResults <- snprank(rgResult, alg.gamma = paramGamma)
  snprankResults
}

#' Rank by Relief-F.
#' 
#' \code{rankRelieff} 
#' 
#' @family feature selection functions
#' @keywords relief nearest neighbors interactions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @param corelearn.est \code{character} CORElearn estimator infoCore(what="attrEval")
#' @return \code{data.frame} with Relief-F results: variable, score.
#' @examples
#' data(testdata100ME4)
#' rankRelieffResults <- rankRelieff(testdata100ME4)
#' @export
rankRelieff <- function(labelledDataFrame, corelearn.est = "ReliefFbestK") {
  labelledDataFrame$Class <- factor(labelledDataFrame$Class, levels = c(0, 1))
  relieffRankings <- CORElearn::attrEval(Class ~ ., labelledDataFrame, 
                                         estimator = corelearn.est)
  retDF <- data.frame(variable = names(relieffRankings), 
                      score = relieffRankings,
                      stringsAsFactors = FALSE)
  retDF <- retDF[order(retDF$score, decreasing = TRUE), ]
  retDF
}

#' Rank by t-test.
#' 
#' \code{rankTTest} 
#' 
#' @keywords models univariate
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column,
#' @return \code{data.frame} with t-test results: variable, p-value.
#' @examples
#' data(testdata10)
#' rankTTestResults <- rankTTest(testdata10)
#' @export
rankTTest <- function(labelledDataFrame) {
  numPredictors <- ncol(labelledDataFrame) - 1
  # make predictors
  predictors <- as.matrix(labelledDataFrame[, 1:numPredictors])
  varNames <- colnames(predictors)
  predictors <- t(predictors)
  # phenotype is response
  response <- labelledDataFrame[, ncol(labelledDataFrame)]
  stats_matrix <- NULL
  for (i in 1:nrow(predictors)) {
    g1_data <- predictors[i, response == 0]
    g2_data <- predictors[i, response == 1]
    t_result <- t.test(g1_data, g2_data)
    stats_matrix <- rbind(stats_matrix, t_result$p.value)
  }
  rownames(stats_matrix) <- varNames
  returnDF <- data.frame(variable = varNames, 
                         score = stats_matrix[, 1], 
                         stringsAsFactors = FALSE)
  returnDF <- returnDF[order(returnDF$score, decreasing = FALSE), ]
  returnDF
}

#' Rank variables by univariate regression.
#' 
#' \code{rankUnivariateRegression}
#' 
#' @keywords models univariate regression
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} with variable in columns and samples in rows;
#' the last column should be binary and labelled 'Class'.
#' @param numCores \code{numeric} number of processor cores to use in mclapply
#' @param sortOrder \code{character} sort retruned data frame by named column
#' @return \code{data.frame} with variable, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns.
#' @examples
#' data(testdata10)
#' rankUnivariateRegressionResults <- rankUnivariateRegression(testdata10)
#' @export
rankUnivariateRegression <- function(labelledDataFrame, numCores=2, sortOrder="pval") {
  # calculate the logistic regression coefficient for each variable in 
  # labelledDataFrame versus the "Class" column (last column)
  colNames <- colnames(labelledDataFrame)[1:(ncol(labelledDataFrame) - 1)]
  numVars <- length(colNames)
  labelledDataFrame$Class <- factor(labelledDataFrame$Class)
  results <- parallel::mclapply(1:numVars, 
                                mc.cores = numCores, 
                                FUN = function(i) {
                                  glmFormula <- paste("Class ~ ", colNames[i], sep = "")
                                  interactionModel <- glm(as.formula(glmFormula), 
                                                          family = "binomial", 
                                                          data = labelledDataFrame)
                                  glmConverged <- interactionModel$converged
                                  fitVarStats <- broom::tidy(interactionModel)[2, ]
                                  data.frame(Variable = colNames[i], 
                                             Converged = glmConverged, 
                                             Beta = round(fitVarStats$estimate, 6), 
                                             p = round(fitVarStats$p.value, 6), 
                                             SE = round(fitVarStats$std.error, 6), 
                                             Stat = round(fitVarStats$statistic, 6), 
                                             stringsAsFactors = FALSE)
                                })
  # put all the parallel results into a data frame
  resultsDF <- do.call(rbind, results)
  
  # return the sorted results complying with the feature selection interface
  returnDF <- NULL
  if (sortOrder == "pval") {
    sortedResults <- resultsDF[order(resultsDF$p), ]
    returnDF <- data.frame(variable = sortedResults$Variable, 
                           score = sortedResults$p)
  }
  if (sortOrder == "coef") {
    sortedResults <- resultsDF[order(resultsDF$Beta, decreasing = TRUE), ]
    returnDF <- data.frame(variable = sortedResults$Variable, 
                           score = sortedResults$Beta)
  }
  if (sortOrder == "stat") {
    sortedResults <- resultsDF[order(resultsDF$Stat, decreasing = TRUE), ]
    returnDF <- data.frame(variable = sortedResults$Variable, 
                           score = sortedResults$Stat)
  }
  returnDF  
}

#' Implements the GeneRank algorithm.
#'
#' http://www.biodatamining.org/content/8/1/2
#' There are three implementations based on the number of genes.
#' 
#' \code{geneRank} 
#' 
#' @keywords models eigenvectors
#' @family feature selection functions
#' @param X \code{matrix} n rows of genes by m columns of subjects.
#' @param a \code{numeric} probability?.
#' @param eps \code{numeric} tolerance.
#' @param maxit \code{numeric} maximum nuber of iterations.
#' @param verbose \code{logical} to send messages to stdout.
#' @return ordered \code{data.frame} with variable and score, ordered by score.
#' @examples
#' data(testdata10)
#' geneRankResults <- rankGeneRank(testdata10)
#' @export
geneRank <- function(X, a=0.9, eps=0.0001, maxit=10, verbose=FALSE) {
  n <- nrow(X)
  if (n <= 5000) {
    if (verbose) cat("Calling geneRank with built-in R function eigen()\n")
    eigenvec <- leftMaxEigen1(X, a)
  } else {
    if (n <= 10000) {
      if (verbose) cat("Calling geneRank with power method\n")
      eigenvec <- leftMaxEigen2(X, a, eps, maxit)
    } else {
      if (verbose) cat("Calling geneRank with power method, no R2 matrix storage\n")
      eigenvec <- leftMaxEigen3(X, a, eps, maxit)
    }
  }
  # return standard data frame for ranker algorithms
  returnDF <- data.frame(variable = rownames(X), 
                         score = eigenvec, 
                         stringsAsFactors = FALSE)
  returnDF <- returnDF[order(returnDF$score, decreasing = TRUE), ]
  returnDF
}

# From the paper Appendix - Computation of the left maximum eigenvector

#' 1. Using R's built in eigen() function
#'
#' NOTE: I had to adapt this code for matrix operations from the copy-and-pasted
#' version from the PDF. Looks like it was orginally Matlab "converted" to R syntax?
#'
#' \code{leftMaxEigen2} 
#'
#' @param X \code{matrix} n rows of genes by m columns of subjects.
#' @param a \code{numeric} probability?.
#' @return left maximum eigenvector.
#' @keywords internal
leftMaxEigen1 <- function(X, a=0.9) {
  n <- nrow(X)
  R2 <- cor(t(X))^2
  sum.row <- rowSums(R2)
  R2.star <- 1/sum.row * R2
  H <- (1 - a) / n + a * R2.star
  p <- abs(eigen(t(H))$vectors[, 1])
  return(p)
}

#' 2. The power method for moderate 5000 < n < 10000.
#'
#' NOTE: I had to adapt this code for matrix operations from the copy-and-pasted
#' version from the PDF. Looks like it was orginally Matlab "converted" to R syntax?
#'
#' \code{leftMaxEigen2} 
#' 
#' @param X \code{matrix} n rows of genes by m columns of subjects.
#' @param a \code{numeric} probability?.
#' @param eps \code{numeric} epsilon tolerance stopping criteria.
#' @param maxit \code{numeric} maximum number of iterations.
#' @return left maximum eigenvector.
#' @keywords internal
leftMaxEigen2 <- function(X, a=0.9, eps=0.0001, maxit=10) {
  #m <- ncol(X)
  n <- nrow(X)
  R2 <- cor(t(X))^2
  sum.row <- rowSums(R2)
  R2.star <- 1/sum.row * R2
  H <- (1 - a)/n + a * R2.star
  tH = t(H)
  p = sum.row/sqrt(sum(sum.row^2))
  # make iteratively better
  for (it in 1:maxit) {
    p.new = tH %*% p
    p.new = p.new / sqrt(sum(p.new^2))
    adiff = max(abs(p - p.new))
    if (adiff < eps) break
    p = p.new[, 1]
  }
  return(p.new[, 1])
}

#' 3. for large n > 10000.
#'
#' NOTE: This code worked directly from cut-and-paste from the PDF file.
#'
#' \code{leftMaxEigen3} 
#' 
#' @param X \code{matrix} n rows of genes by m columns of subjects.
#' @param a \code{numeric} probability?.
#' @param eps \code{numeric} epsilon tolerance stopping criteria.
#' @param maxit \code{numeric} maximum number of iterations.
#' @return left maximum eigenvector.
#' @keywords internal
leftMaxEigen3 <- function(X, a = 0.9, eps = 0.0001, maxit = 10) {
  m <- ncol(X)
  n <- nrow(X)
  # compute normalized gene expression matrix
  x.bar = rowMeans(X)
  Xsub.mean = X - x.bar
  sdX = sqrt(rowSums(Xsub.mean^2))
  Z = (1/sdX)*Xsub.mean
  # compute Rstar^2'p without computing Rstar^2
  sumR2 = rep(0, n)
  for (i in 1:m) {
    for (j in 1:m)
    {
      qij = sum(Z[, i]*Z[, j])
      sumR2 = sumR2 + qij * Z[, i] * Z[, j]
    }
  }
  p = sumR2/sqrt(sum(sumR2^2))
  # make iteratively better
  for (it in 1:maxit) {
    tR2p.fast = rep(0,n)
    for (i in 1:m)
      for (j in 1:m)
      {
        hij = sum(Z[,i]*Z[,j]*p/sumR2)
        tR2p.fast = tR2p.fast + hij * Z[, i] * Z[, j]
      }
    p.new = (1 - a) / n * sum(p) + a * tR2p.fast
    p.new = p.new / sqrt(sum(p.new^2))
    adiff = max(abs(p - p.new))
    if (adiff < eps) break
    p = p.new
  }
  return(p.new)
}

#' Relative recurrency variable importance metric (r2VIM).
#' 
#' There is a new version. This is based on the earlier PSB 2015 version.
#' 
#' \code{r2VIMorig} 
#' 
#' @keywords models
#' @family feature selection functions
#' @param predictors \code{matrix} independent variables.
#' @param response \code{vector} response vector, case-control.
#' @param numRfRuns \code{numeric} of randomForest runs.
#' @param thresholdVIM \code{numeric} threshold importance.
#' @param thresholdMedianRIS \code{numeric} threshold importance RIS score.
#' @param thresholdProb \code{numeric} threshold probability seen in top 10.
#' @param verbose \code{logical} write verbose messages to console/stdout.
#' @return \code{list} with: run statistics, votes matrix, network matrix, 
#' importance distribution, importance distribution RIS.
#' @examples
#' data(testdata10)
#' predictors <- as.matrix(testdata10[, -ncol(testdata10)])
#' response <- factor(testdata10[, ncol(testdata10)])
#' r2vimResults <- r2VIMorig(predictors = predictors, response = response, verbose = TRUE)
#' @export
r2VIMorig <- function(predictors = NULL, 
                      response = NULL, 
                      numRfRuns = 10,
                      thresholdVIM = 0,
                      thresholdMedianRIS = 2,
                      thresholdProb = 0.2,
                      verbose = FALSE) {
  if (is.null(predictors ) || is.null(response)) {
    stop("r2VIMorig: predictors and response are required parameters")
  }
  numvariables <- ncol(predictors)
  variableIDs <- colnames(predictors)
  # compute the variable importances using random forests numRfRuns times,
  # accumulating the importance vectors into a "Dist"ribution data frames
  importanceDistribution <- NULL
  importanceDistributionRIS <- NULL
  # variable votes
  votes <- as.data.frame(matrix(nrow = 1, ncol = numvariables, data = c(0)))
  colnames(votes) <- variableIDs
  # GAIN matrix for interaction votes
  networkMatrix <- matrix(ncol = numvariables, nrow = numvariables, data = c(0))
  colnames(networkMatrix) <- variableIDs
  rownames(networkMatrix) <- variableIDs
  if (verbose) cat("Running", numRfRuns, "randomForest models\n")
  runStats <- NULL
  for (testIdx in 1:numRfRuns) {
    # determine variable importance with random forest
    rfResult <- randomForest::randomForest(predictors, response, importance = TRUE)
    # accumulate importances
    thisRunImportances <- rfResult$importance[, 3]
    estimatedVariance <- abs(min(thisRunImportances))
    thisRunImportancesRIS <- thisRunImportances / estimatedVariance
    importanceDistribution <- rbind(importanceDistribution, thisRunImportances)
    importanceDistributionRIS <- rbind(importanceDistributionRIS, thisRunImportancesRIS)
    # count top ten
    top10variables <- names(sort(thisRunImportances, decreasing = TRUE)[1:10])
    votes[1, top10variables] <- votes[1, top10variables] + 1
    # count all pairs in top 10
    for (pair1idx in 1:9) {
      for (pair2idx in (pair1idx + 1):10) {
        variable1 <- top10variables[pair1idx]
        variable2 <- top10variables[pair2idx]
        # add 1/45 to each pair? 45 = choose(10, 2)
        networkMatrix[variable1, variable2] <- networkMatrix[variable1, variable2] + 1
        networkMatrix[variable2, variable1] <- networkMatrix[variable1, variable2] + 1
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
    runStats <- rbind(runStats, data.frame(run = testIdx, 
                                           confusion = confusionMatrix,
                                           sensitivity = SENS,
                                           specificity = SPEC))
    if (verbose) cat(testIdx, " / ", numRfRuns, ": SENS/SPEC: ", round(SENS, 6), 
                     " / ", round(SPEC, 6), "\n", sep = "")
  }
  diag(networkMatrix) <- as.numeric(votes[1, ])
  
  list(run.stats = runStats,
       votes = votes,
       net.matrix = networkMatrix, 
       importance.dist = importanceDistribution, 
       importance.dist.ris = importanceDistributionRIS)
}

#' Rank variables by SNPrank algorithm.
#' 
#' This function is for backward compatibility with the original Rinbix
#' with a single gamma value. 
#' 
#' \code{snprank}
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.nature.com/variable/journal/v11/n8/full/gene201037a.html}
#'   {Genes & Immunity Paper}
#'   \item \url{http://insilico.utulsa.edu/index.php/snprank/}{SNPrank Software}
#' }
#' @family feature selection functions
#' @family GAIN functions
#' @family inbix synonym functions
#' @param G \code{matrix} genetic association interaction network.
#' @param alg.gamma \code{numeric} weighting, interactions (closer to 1) versus 
#' main effects (closer to 0).
#' @return sortedTable \code{data.frame} with variable, SNPrank, diagonal and 
#' degree columns.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' rinbixSnpranksDF <- snprank(rinbixRegain)
#' @export
snprank <- function(G, alg.gamma=0.85) {
  n <- nrow(G);
  variableNames <- colnames(G)
  Gdiag <- diag(G)
  Gtrace <- sum(Gdiag)
  #colsum <- colSums(G);
  diag(G) <- 0
  Gtrace <- Gtrace * n;
  colsumG <- colSums(G) + 0.001
  #rowSumG <- rowSums(G)
  rowsum_denom <- matrix(0, n, 1);
  for (i in 1:n) {
    localSum <- 0;
    for (j in 1:n) {
      factor <- ifelse(G[i, j] != 0, 1, 0);
      localSum <- localSum + (factor * colsumG[j]);
    }
    rowsum_denom[i] <- localSum;
  }
  gamma_vec <- rep(alg.gamma, n);
  gamma_matrix <- matrix(nrow = n, ncol = n, data = rep(gamma_vec, n))
  if (Gtrace) {
    b <- ((1.0 - gamma_vec) / n) + (Gdiag / Gtrace)  
  } else {
    b <- ((1.0 - gamma_vec) / n)
  }
  D <- matrix(nrow = n, ncol = n, data = c(0))
  diag(D) <- 1 / colsumG
  I <- diag(n)
  temp <- I - gamma_matrix * G %*% D
  r <- solve(temp, b)
  snpranks <- r / sum(r)
  saveTable <- data.frame(variable = variableNames, 
                          score = snpranks,
                          stringsAsFactors = FALSE)
  sortedTable <- saveTable[order(saveTable$score, decreasing = TRUE),]  
  sortedTable
}
