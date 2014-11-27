# inbix.R - Bill White - 1/30/14
#
# An R library that mimics many functions found in the Insilico 
# Bioinformatics (inbix) C++ library.

library(parallel)

#' A random dataset with 10 samples and 10 variables
#'
#' A dataset containing random variables var1 ... var10. 
#' Class/phenotypes are 5 0's and 5 1's.
#' The variables are as follows:
#'
#' @docType data
#' @keywords datasets
#' @name testdata10
#' @usage data(testdata10)
#' @format A data frame with 10 rows and 10 variables plus "Class" column
NULL

# ----------------------------------------------------------------------------
#' Create a random regression data set with binary class
#' 
#' \code{createRandomRegressionDataset} 
#' 
#' @param numRows number of rows (samples)
#' @param numCols number of columns (independent variables)
#' @return data frame
#' @export
createRandomRegressionDataset <- function(numRows, numCols) {
  dmatrix <- matrix(nrow=numRows, ncol=numCols, data=rnorm(numRows*numCols))
  dpheno <- c(rep(0, numRows/2), rep(1, numRows/2))
  dataset <- cbind(dmatrix, dpheno)
  colnames(dataset) <- c(paste("var", 1:numCols, sep=""), "Class")
  as.data.frame(dataset)
}

# ----------------------------------------------------------------------------
#' Read inbix numeric data set and phenotype file; return combined data frame
#' 
#' \code{readInbixNumericAsRegressionData} 
#' 
#' @param baseInbixName file base name
#' @return data frame with numeric data in first m columns and phenotye in m+1 column
#' @export
readInbixNumericAsRegressionData <- function(baseInbixName) {
  inbixNumericFile <- paste(baseInbixName, ".num", sep="")
  inbixNumericTable <- read.table(inbixNumericFile, header=T, sep="\t")
  inbixNumericTable <- inbixNumericTable[, 3:ncol(inbixNumericTable)]
  
  inbixPhenoFile <- paste(baseInbixName, ".pheno", sep="")
  inbixPhenoTable <- read.table(inbixPhenoFile, header=F, sep="\t")[, 3]
  
  regressionData <- cbind(inbixNumericTable, inbixPhenoTable)
  colnames(regressionData) <- c(colnames(inbixNumericTable), "Class")
  regressionData
}

# ----------------------------------------------------------------------------
#' Write regression data frame as  inbix numeric data set and phenotype file
#' 
#' \code{writeRegressionDataAsInbixNumeric} 
#' 
#' @param regressionData is a data frame with sample rows, gene columns 
#' plus phenotype column
#' @param baseInbixName file base name
#' @export
writeRegressionDataAsInbixNumeric <- function(regressionData, baseInbixName) {
  numSamples <- nrow(regressionData)
  numGenes <- ncol(regressionData)-1

  phenoCol <- ncol(regressionData)
  phenos <- regressionData[, phenoCol]
  subIds <- c(paste("subj", 1:numSamples, sep=""))
  phenosTable <- cbind(subIds, subIds, phenos)
  inbixPhenoFile <- paste(baseInbixName, ".pheno", sep="")
  write.table(phenosTable, inbixPhenoFile, quote=F, sep="\t", 
              col.names=F, row.names=F)
  
  inbixNumericData <- regressionData[, 1:numGenes]
  inbixNumericTable <- cbind(subIds, subIds, inbixNumericData)
  colnames(inbixNumericTable) <- c("FID", "IID", colnames(inbixNumericData))
  inbixNumericFile <- paste(baseInbixName, ".num", sep="")
  write.table(inbixNumericTable, inbixNumericFile, quote=F, sep="\t", 
              col.names=T, row.names=F)
}

# ----------------------------------------------------------------------------
#' Write simulated differential expression data matrix as inbix format
#' 
#' \code{writeSimulatedInbixDataset} 
#' 
#' @param D is a simulated differential coexpression data matrix
#' @param Dfileprefix is an output data file prefix
#' @export
# 
writeSimulatedInbixDataset <- function(D, Dfileprefix) {
  # ----------------------------------------------------------------------------
  # write the data set in inbix format for dcGAIN + SNPRank
  # save inbix phenotypes
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("case", 1:n1, sep=""), paste("ctrl", 1:n2, sep=""))
  phenos <- c(rep(1, n1), rep(0, n2))
  phenosTable <- cbind(subIds, subIds, phenos)
  datasimInbixPhenoFile <- paste(Dfileprefix, ".pheno", sep="")
  write.table(phenosTable, datasimInbixPhenoFile, quote=F, sep="\t", 
              col.names=F, row.names=F)
  # save inbix numeric file/data set
  inbixData <- t(D)
  dataTable <- cbind(subIds, subIds, inbixData)
  geneNames <- paste("gene", sprintf("%04d", 1:nrow(D)), sep="")
  colnames(dataTable) <- c("FID", "IID", geneNames)
  datasimInbixNumFile <- paste(Dfileprefix, ".num", sep="")
  write.table(dataTable, datasimInbixNumFile, quote=F, sep="\t", 
              col.names=T, row.names=F)
}

# ----------------------------------------------------------------------------
#' Fisher correlation transformation function R to Z distribution
#' 
#' \code{fiserRtoZ}
#' 
#' @param x a correlation value -1 to 1
#' @return transformed correlation value
#' @export
# 
fisherRtoZ <- function(x) { 
  .5 * log(abs((1 + x) / (1 - x))) 
}

# -----------------------------------------------------------------------------
#' Remove gene profiles with low absolute values
#' 
#' \code{genelowvalfilter} removes genes with values in the lowest 10%
#' 
#' @param dataMatrix a data frame with genes in rows and samples in columns
#' @return a list with the mask used and filtered data frame
#' @export
genelowvalfilter <- function(dataMatrix) {
  # Remove gene profiles with low absolute values in dataMatrix. Returns:
  # 1) a logical vector mask identifying gene expression profiles in dataMatrix
  #    that have absolute expression levels in the lowest 10% of the data set.
  # 2) a data matrix containing filtered expression profiles.
  threshold <- quantile(dataMatrix, c(0.1))
  mask <- apply(dataMatrix, 1, function(x) all(x < threshold))
  fdata <- dataMatrix[!mask, ]
  
  # return the row mask and filtered data
  list(mask=!mask, fdata=fdata)
}

# -----------------------------------------------------------------------------
#' Remove gene profiles with low expression variance
#' 
#' \code{genelowvarfilter} removes genes with variance below a 
#' variance percentile threshold
#' 
#' @param dataMatrix a data frame with genes in rows and samples in columns
#' @param percentile variance threshold below which genes will be removed
#' @return a list with the mask used and filtered data frame
#' @export
genevarfilter <- function(dataMatrix, percentile) {
  # calculates the variance for each gene expression profile in Data and returns Mask, 
  # which identifies the gene expression profiles with a variance less than the 10th 
  # percentile. Mask is a logical vector with one element for each row in Data. 
  # The elements of Mask corresponding to rows with a variance greater than the threshold 
  # have a value of 1, and those with a variance less than the threshold are 0.
  probeVariances <- apply(dataMatrix, 1, var)
  threshold <- quantile(probeVariances, c(percentile))
  mask <- apply(dataMatrix, 1, function(x) var(x) > threshold)
  fdata <- dataMatrix[mask, ]
  
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

# -----------------------------------------------------------------------------
#' Rank genes by univariate regression
#' 
#' \code{rankUnivariateRegression}
#' 
#' @param regressionData a data frame with gene in columns and samples in rows;
#' the last column should be labeled 'Class' and be 0 or 1 values
#' @return data frame with gene, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns
#' @export
rankUnivariateRegression <- function(regressionData) {
  # calculate the logistic regression coefficient for each variable in 
  # regressionData versus the "Class" column (last column)
  colNames <- colnames(regressionData)[1:(ncol(regressionData)-1)]
  numVars <- length(colNames)
  results <- NULL
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
  sortedResults <- results[order(results[,4]),]
}

# -----------------------------------------------------------------------------
#' Get the main effect of a variable using generalized linear regression - glm
#' 
#' \code{getMainEffect} 
#' 
#' @param data is the data frame with genes in columns and samples in rows
#' @param geneName is the column name of the gene used in the regression formula
#' @param depVarName is the phenotype column name used in the regression formula
#' @param regressionFamily is the regression family name for the glm
#' @return the regression beta coefficient 
#' @export
getMainEffect <- function(data, geneName, depVarName, regressionFamily) {
  regressionFormula <- as.formula(paste(depVarName, "~", paste("`", geneName, "`", sep=""), sep=" "))
  mainEffectModel <- glm(regressionFormula, family=regressionFamily, data=data)
  
  as.numeric(mainEffectModel$coef[2])
}

# -----------------------------------------------------------------------------
#' Parallel execution of the regression genetic association network algorithm
#' 
#' \code{regainParallel}
#' 
#' @param regressionData a data frame with gene in columns and samples in rows;
#' the last column should be labeled 'Class' and be 0 or 1 values
#' @param stdBetas flag to use standardized beta coefficients
#' @param absBetas flag to use absolute value of beta coefficients
#' @return regainMatrix a matrix of gene by gene regression coefficients
#' @export
regainParallel <- function(regressionData, stdBetas=FALSE, absBetas=FALSE) {
  transform <- ifelse(absBetas, "abs", "")
  rawBetas <- ifelse(stdBetas, FALSE, TRUE)
  mainEffects <- getMainEffects(regressionData, 
                                useBetas=rawBetas,
                                transformMethod=transform)
  regainMatrix <- getInteractionEffects(regressionData, 
                                        useBetas=rawBetas,
                                        transformMethod=transform)
  diag(regainMatrix) <- mainEffects
  colnames(regainMatrix) <- colnames(regressionData)[1:(ncol(regressionData)-1)]
  # replace NAs with zero
  regainMatrix[is.na(regainMatrix)] <- 0
  regainMatrix
}

# -----------------------------------------------------------------------------
#' Get main effects from generalized linear model regression (parallel)
#' 
#' \code{getMainEffects} 
#' 
#' @param data data is the data frame with genes in columns and samples in rows
#' @param regressionFamily glm regression family name
#' @param numCovariates number of included covariates
#' @param writeBetas flag indicating whther to write beta values to separate file
#' @param useBetas flag indicating betas rather than standardized betas used
#' @param transformMethod is an optional transform method
#' @param verbose flag to send verbose messages to stdout
#' @return mainEffectValues a vector of main effect values
#' @export
# -----------------------------------------------------------------------------
getMainEffects <- function(data, regressionFamily="binomial", numCovariates=0, 
                           writeBetas=FALSE, useBetas=FALSE, transformMethod="", 
                           verbose=FALSE) {
  allColumnNames <- colnames(data)
  numColumnNames <- length(allColumnNames)
  depVarName <- "Class"
  endVariableNames <- ncol(data) - numCovariates - 1
  variableNames <- allColumnNames[1:endVariableNames]
  numVariables <- length(variableNames)
  cat(paste("Computing GLM main effect models for each variable (", 
            numVariables,	"), in parallel", sep=""), "\n")
  results <- mclapply(variableNames, 
                      function(x) runMainEffectsTest(data, x, depVarName, 
                                                     regressionFamily, numCovariates))
  cat("Loading the reGAIN matrix diagonal with GLM main effect coefficients\n")
  #print(unlist(results))
  #stop("DEBUG")
  if(writeBetas) {
    betaInfo <- NULL
  }
  mainEffectValues <- vector(mode="numeric", length=numVariables)
  for(i in 1:numVariables) {
    regressionModel <- results[[i]]
    #print(regressionModel)
    #print(summary(regressionModel))
    glmConverged <- regressionModel$converged
    if(!glmConverged) {
      cat("WARNING: Regression model failed to converge for", variableNames[i], "\n")
    }
    mainCoeff <- regressionModel$coef
    mainStdErr <- regressionModel$se
    mainPval <- regressionModel$pval
    mainStat <- regressionModel$stdbeta
    mainEffectValue <- NA
    if(glmConverged) {
      if(useBetas) {
        mainEffectValue <- mainCoeff
      }
      else {
        mainEffectValue <- mainStat
      }
    }
    mainEffectValueTransformed <- mainEffectValue
    if(!is.na(mainEffectValue)) {
      if(transformMethod == "abs") {
        mainEffectValueTransformed <- abs(mainEffectValue)
      }
      if(transformMethod == "threshold") {
        mainEffectValueTransformed <- ifelse(mainEffectValue < 0, 0, mainEffectValue)
      }
    }
    #regainMatrix[i, i] <- mainEffectValueTransformed
    mainEffectValues[i] <- mainEffectValueTransformed
    if(writeBetas) {
      thisBetas <- c(mainCoeff, mainPval)
      if(numCovariates > 0) {
        for(cn in 1:numCovariates) {
          covBeta <- summary(regressionModel)$coefficients[2+cn, "Estimate"]
          covPval <- summary(regressionModel)$coefficients[2+cn, "Pr(>|z|)"]
          thisBetas <- c(thisBetas, covBeta, covPval)
        }
      }
      betaInfo <- rbind(betaInfo, thisBetas)
    }
  }
  if(writeBetas) {
    cat("Writing betas to mebetas.tab\n") 
    colNames <- c("Main Effect", "Pval")
    if(numCovariates > 0) {
      for(cn in 1:numCovariates) {
        colNames <- c(colNames, c(paste("COV", cn, sep=""), "Pval"))
      }
    }
    colnames(betaInfo) <- colNames
    rownames(betaInfo) <- variableNames
    write.table(betaInfo, file="mebetas.tab", quote=FALSE, sep="\t")
    rm(betaInfo)
  }
  
  # clean up memory
  rm(results)
  
  mainEffectValues
}

# -----------------------------------------------------------------------------
#' Get the main effect of a variable using generalized linear regression - glm
#' 
#' \code{runMainEffectsTest}
#' 
#' @param data data frame with genes in columns and samples in rows
#' @param variableName name of the variable to consider
#' @param depVarName name of the phenotype variable
#' @param regressionFamily glm regression family name
#' @param numCovariates number of included covariates
#' @return data frame with gene, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns
#' @export
runMainEffectsTest <- function(data, variableName, depVarName, regressionFamily, numCovariates) {
  if(numCovariates > 0) {
    covarsStart <- ncol(data) - numCovariates + 1
    covarNames <- colnames(data)[covarsStart:ncol(data)]
    
    covarsModelParts <- paste(covarNames, collapse=" + ")
    formulaString <- paste(depVarName, " ~ ", variableName, 
                           " + ", covarsModelParts, sep="")
    regressionFormula <- as.formula(formulaString)
  } else {
    regressionFormula <- as.formula(paste(depVarName, "~", 
                                          paste("`", variableName, "`", sep=""), sep=" "))
  }
  regressionModel <- glm(regressionFormula, family=regressionFamily, data=data)
  #print(regressionModel)
  
  #interceptCoeff <- summary(regressionModel)$coefficients[1, "Estimate"]
  mainCoeff <- summary(regressionModel)$coefficients[2, "Estimate"]
  mainStdErr <- summary(regressionModel)$coefficients[2, "Std. Error"]
  if(regressionFamily == "binomial") {
    #interceptPval <- summary(regressionModel)$coefficients[1, "Pr(>|z|)"]
    mainPval <- summary(regressionModel)$coefficients[2, "Pr(>|z|)"]
    mainStat <- summary(regressionModel)$coefficients[2, "z value"]
  }
  else {
    #interceptPval <- summary(regressionModel)$coefficients[1, "Pr(>|t|)"]
    mainPval <- summary(regressionModel)$coefficients[2, "Pr(>|t|)"]
    mainStat <- summary(regressionModel)$coefficients[2, "t value"]
  }

  data.frame(converged=regressionModel$converged, coef=mainCoeff, se=mainStdErr, 
             pval=mainPval, stdbeta=mainStat)
}

# -----------------------------------------------------------------------------
#' Get interaction effects from generalized linear model regression
#' 
#' \code{getInteractionEffects} 
#' 
#' @param data is the data frame with genes in columns and samples in rows
#' @param regressionFamily glm regression family name
#' @param numCovariates the number of included covariates
#' @param writeBetas flag indicating whther to write beta values to separate file
#' @param excludeMainEffects flag indicating whether to exclude main effect terms
#' @param useBetas flag indicating betas rather than standardized betas used
#' @param transformMethod is an optional transform method
#' @param verbose flag to send verbose messages to stdout
#' @return results a matrix of gene by gene regression coefficients
#' @export
getInteractionEffects <- function(data, regressionFamily="binomial", numCovariates=0,
                                  writeBetas=FALSE, excludeMainEffects=FALSE, useBetas=FALSE, 
                                  transformMethod="", verbose=FALSE) {
  allColumnNames <- colnames(data)
  numColumnNames <- length(allColumnNames)
  depVarName <- "Class"
  endVariableNames <- ncol(data) - numCovariates - 1
  variableNames <- allColumnNames[1:endVariableNames]
  numVariables <- length(variableNames)
  cat("From", numVariables ,"generating", choose(numVariables,2), 
      "variable index pairs for reGAIN upper triangular matrix\n")
  idxCombList <- combn(numVariables, 2, list)
  #print(idxCombList)
  cat("Computing GLM interaction models for each index pair, in parallel\n")
  lastIdx <- length(idxCombList)
  numSplits <- 10
  splitSize <- as.integer(lastIdx / numSplits) + 1
  results <- NULL
  for(i in 1:numSplits) {
    startIdx <- (i - 1) * splitSize + 1
    endIdx <- i *  splitSize
    if(endIdx > lastIdx) {
      endIdx <- lastIdx
    }
    if(startIdx < lastIdx) {
      cat("Running chunk", i, "split size:", splitSize, "start:", startIdx, 
          "end:", endIdx, "\n")
      results <- c(results, mclapply(idxCombList[startIdx:endIdx], 
                                     function(x) runInteractionEffectsTest(data, x, 
                                                                           depVarName, 
                                                                           regressionFamily, 
                                                                           numCovariates, 
                                                                           excludeMainEffects)))
    }
  }
  if(verbose) {
    print(results)
  }
  cat("Loading the reGAIN matrix upper and lower triangulars with",
      "GLM interaction coefficients\n")
  if(writeBetas) {
    betaInfo <- NULL
    betaRows <- NULL
  }
  interactionValues <- matrix(nrow=numVariables, ncol=numVariables)
  for(i in 1:length(idxCombList)) {
    thisComb <- idxCombList[[i]]
    variable1Idx <- thisComb[1]
    variable2Idx <- thisComb[2]
    interactionName <- paste(variableNames[variable1Idx], "x", 
                             variableNames[variable2Idx])
    
    regressionModel <- results[[i]]
    glmConverged <- regressionModel$converged
    if(!glmConverged) {
      cat("WARNING: Regression model failed to converge for", interactionName, "\n")
    }
    interactionCoeff <- regressionModel$coef
    interactionStdErr <- regressionModel$se
    interactionPval <- regressionModel$pval
    interactionStat <- regressionModel$stdbeta
    
    interactionValue <- NA
    if(glmConverged) {
      if(useBetas) {
        interactionValue <- interactionCoeff
      }
      else {
        interactionValue <- interactionStat
      }
    }
    interactionValueTransformed <- interactionValue
    if(!is.na(interactionValue)) {
      if(transformMethod == "abs") {
        interactionValueTransformed <- abs(interactionValue)
      }
      if(transformMethod == "threshold") {
        interactionValueTransformed <- ifelse(interactionValue < 0, 0, interactionValue)
      }
    }
    interactionValues[variable1Idx, variable2Idx] <- interactionValueTransformed
    interactionValues[variable2Idx, variable1Idx] <- interactionValueTransformed
    # 		cat("Interaction z value:", variable1Idx, variable2Idx, interactionStat, "\n")
    #     if(abs(interactionStat) > 4) {
    #       print(summary(regressionModel))
    #     }
    if(writeBetas) {
      allBetas <- NULL
      for(colIdx in 1:length(names(regressionModel$coefficients))) {
        varCoef <- summary(regressionModel)$coefficients[colIdx, "Estimate"]
        if(regressionFamily == "binomial") {
          varPval <- summary(regressionModel)$coefficients[colIdx, "Pr(>|z|)"]
        }
        else {
          varPval <- summary(regressionModel)$coefficients[colIdx, "Pr(>|t|)"]
        }
        #cat(colIdx, varCoef, varPval, "\n")
        allBetas <- c(allBetas, varCoef, varPval)
      }
      betaInfo <- rbind(betaInfo, allBetas)
      betaRows <- c(betaRows, interactionName)
    }
  }
  if(writeBetas) {
    cat("Writing betas to intbetas.tab\n") 
    # always
    betaCols <- c("Intercept", "Pval")
    # if main effects included
    if(!excludeMainEffects) {
      betaCols <- c(betaCols, paste(c("Main Effect 1", "Pval", "Main Effect 2", "Pval"), sep="\t"))
    }
    # if covariates
    if(numCovariates > 0) {
      for(cn in 1:numCovariates) {
        betaCols <- c(betaCols, c(paste("COV", cn, sep=""), "Pval"))
      }
    }
    # always
    betaCols <- c(betaCols, c("Interaction", "Pval"))
    # print(betaCols)
    colnames(betaInfo) <- betaCols
    rownames(betaInfo) <- betaRows
    write.table(betaInfo, file="intbetas.tab", quote=FALSE, sep="\t")
    rm(betaInfo)
  }
  
  # clean up memory
  rm(idxCombList)
  rm(results)
  
  interactionValues
}

# -----------------------------------------------------------------------------
#' Get the interaction effect of a pair of variables using generalized linear regression - glm
#' 
#' \code{runInteractionEffectsTest} 
#' 
#' @param data is the data frame with genes in columns and samples in rows
#' @param variableIndices column indices of variable pairs
#' @param depVarName name of the phenotype variable
#' @param regressionFamily glm regression family name
#' @param numCovariates number of covariates included
#' @param excludeMainEffects flag indicating whether to exclude main effect terms
#' @return data frame with gene, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns
#' @export
runInteractionEffectsTest <- function(data, variableIndices, depVarName, 
                                      regressionFamily, numCovariates, excludeMainEffects) {
  variable1Idx <- variableIndices[1]
  variable2Idx <- variableIndices[2]
  variableNames <- colnames(data)[1:(ncol(data)-1)]
  variable1Name <- variableNames[variable1Idx]  
  variable2Name <- variableNames[variable2Idx]
  if(excludeMainEffects) {
    interactionTerm <- paste("`", variable1Name, "`", ":", "`", 
                             variable2Name, "`", sep="")
  }
  else {
    interactionTerm <- paste("`", variable1Name, "`", "*", "`", 
                             variable2Name, "`", sep="")
  }
  if(numCovariates > 0) {
    covarsStart <- ncol(data) - numCovariates + 1
    covarNames <- colnames(data)[covarsStart:ncol(data)]
    covarsModelParts <- paste(covarNames, collapse=" + ")
    regressionFormula <- as.formula(paste(depVarName, "~", 
                                          paste(interactionTerm, " + ", covarsModelParts, sep=""), sep=" "))
  }
  else {
    regressionFormula <- as.formula(paste(depVarName, "~", interactionTerm, 
                                          sep=" "))
  }
  regressionModel <- glm(regressionFormula, family=regressionFamily, data=data)
  
  if(numCovariates > 0) {
    if(excludeMainEffects) {
      interactionTermIndex <- 1 + numCovariates + 1
    }
    else {
      interactionTermIndex <- 3 + numCovariates + 1
    }
  }
  else {
    if(excludeMainEffects) {
      interactionTermIndex <- 2
    }
    else {
      interactionTermIndex <- 4
    }
  }
  interactionCoeff <- summary(regressionModel)$coefficients[interactionTermIndex, "Estimate"]
  interactionStdErr <- summary(regressionModel)$coefficients[interactionTermIndex, "Std. Error"]
  if(regressionFamily == "binomial") {
    interactionPval <- summary(regressionModel)$coefficients[interactionTermIndex, "Pr(>|z|)"]
    interactionStat <- summary(regressionModel)$coefficients[interactionTermIndex, "z value"]
  }
  else {
    interactionPval <- summary(regressionModel)$coefficients[interactionTermIndex, "Pr(>|t|)"]
    interactionStat <- summary(regressionModel)$coefficients[interactionTermIndex, "t value"]
  }
  
  data.frame(converged=regressionModel$converged, coef=interactionCoeff, 
             stderr=interactionStdErr, pval=interactionPval, stdbeta=interactionStat)
}

# -----------------------------------------------------------------------------
#' Differential coexpression genetic association network algorithm
#' 
#' \code{dcgain} 
#' 
#' @param inbixData is a data frame with samples in rows, genes in columns
#' and phenotype in the last column
#' @return results a matrix of gene by gene differential coexpression values
#' @export
dcgain <- function(inbixData) {
  phenos <- inbixData[, ncol(inbixData)] + 1
  exprBySubj <- inbixData[, 1:(ncol(inbixData)-1)]
  exprByGene <- t(exprBySubj)
  varNames <- colnames(inbixData)[1:ncol(exprBySubj)]
  n1 <- nrow(exprBySubj) / 2
  n2 <- nrow(exprBySubj) / 2
  nVars <- ncol(exprBySubj)

  # determine group correlations
  cor_table_g1 <- cor(t(exprByGene[,phenos==1]))
  cor_table_g2 <- cor(t(exprByGene[,phenos==2]))
    
  # main effect diagonal
  results <- matrix(nrow=nVars, ncol=nVars)
  pvalues <- matrix(nrow=nVars, ncol=nVars)
  for(i in 1:nrow(exprByGene)) {
    g1_data <- exprByGene[i, phenos == 1]
    g2_data <- exprByGene[i, phenos == 2]
    g1_mean <- mean(g1_data)
    g2_mean <- mean(g2_data)
    # z-test
    z_i_1 <- 0.5 * log((abs((1 + g1_mean) / (1 - g1_mean))))
    z_i_2 <- 0.5 * log((abs((1 + g2_mean) / (1 - g2_mean))))
    Z_i <- abs(z_i_1 - z_i_2) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)))
    #Z_i <- abs(g1_mean - g2_mean) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)))
    results[i, i] <- Z_i
    pvalues[i, i] <- 1
    # t-test
#     t_result <- t.test(g1_data, g2_data)
#     t_i <- abs(t_result$statistic)
#     results[i, i] <- t_i
#     pvalues[i, i] <- t_result$p.value
    
    #cat(Z_i, t_i, "\n")
  }
  
  # ----------------------------------------------------------------------------
  # interactions
  for(i in 1:nVars) {
    for(j in 1:nVars) {
      if(j <= i) {
        next
      }
      r_ij_1 <- cor_table_g1[i, j]
      r_ij_2 <- cor_table_g2[i, j]
      z_ij_1 <- 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))))
      z_ij_2 <- 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))))
      Z_ij <- abs(z_ij_1 - z_ij_2) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)))
      pval <- 2 * pnorm(-abs(Z_ij))
      results[i, j] <- Z_ij
      results[j, i] <- Z_ij
      pvalues[i, j] <- pval
      pvalues[j, i] <- pval
    }
  }
  colnames(results) <- varNames
  colnames(pvalues) <- varNames
  
  list(scores=results, pvals=pvalues)
}

# -----------------------------------------------------------------------------
#' Differential modularity genetic association network algorithm
#' 
#' \code{dmgain} 
#' 
#' @param inbixData is a data frame with samples in rows, genes in columns
#' and phenotype in the last column
#' @return results a matrix of gene by gene differential modularity values
#' @export
dmgain <- function(inbixData) {
  phenos <- inbixData[, ncol(inbixData)] + 1
  exprBySubj <- inbixData[, 1:(ncol(inbixData)-1)]
  exprByGene <- t(exprBySubj)
  varNames <- colnames(inbixData)[1:ncol(exprBySubj)]
  n1 <- nrow(exprBySubj) / 2
  n2 <- nrow(exprBySubj) / 2
  nVars <- ncol(exprBySubj)
  
  # determine group correlations
  cor_table_g1 <- cor(t(exprByGene[,phenos == 1]))
  cor_table_g2 <- cor(t(exprByGene[,phenos == 2]))
  
  # main effect diagonal
  results <- matrix(nrow=nVars, ncol=nVars)
  pvalues <- matrix(nrow=nVars, ncol=nVars)
  for(i in 1:nrow(exprByGene)) {
    g1_data <- exprByGene[i, phenos == 1]
    g2_data <- exprByGene[i, phenos == 2]
    # z-test
    g1_mean <- mean(g1_data)
    g2_mean <- mean(g2_data)
    z_i_1 <- 0.5 * log((abs((1 + g1_mean) / (1 - g1_mean))))
    z_i_2 <- 0.5 * log((abs((1 + g2_mean) / (1 - g2_mean))))
    Z_i <- abs(z_i_1 - z_i_2) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)))
    results[i, i] <- Z_i
    pvalues[i, i] <- 1
    # t-test
#     t_result <- t.test(g1_data, g2_data)
#     t_i <- abs(t_result$statistic)
#     results[i, i] <- t_i
#     pvalues[i, i] <- t_result$p.value
  }
  
  # ----------------------------------------------------------------------------
  # interactions
  # added from bam by bcw 7/29/14
  k_1 <- rowSums(cor_table_g1)
  two_m_1 <- sum(k_1)
  k_2 <- rowSums(cor_table_g2)
  two_m_2 <- sum(k_2)
  for(i in 1:nVars) {
    for(j in 1:nVars) {
      if(j <= i) {
        next
      }
      r_ij_1 <- cor_table_g1[i, j]
      r_ij_2 <- cor_table_g2[i, j]
      # added from bam by bcw 7/29/14
      z_ij_1 <- cor_table_g1[i, j] - k_1[i] * k_1[j] / two_m_1
      z_ij_2 <- cor_table_g2[i, j] - k_2[i] * k_2[j] / two_m_2
      Z_ij <- abs(z_ij_1 - z_ij_2) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)))
      pval <- 2 * pnorm(-abs(Z_ij))
      results[i, j] <- Z_ij
      results[j, i] <- Z_ij
      pvalues[i, j] <- pval
      pvalues[j, i] <- pval
    }
  }
  colnames(results) <- varNames
  colnames(pvalues) <- varNames
  
  list(scores=results, pvals=pvalues)
}

# -----------------------------------------------------------------------------
#' Rank genes by SNPrank algorithm
#' 
#' \code{snprank}
#' 
#' @param G is a genetic association network matrix
#' @param gamma is a parameter weighting interactions versus main effects
#' @return sortedTable a data frame with gene, SNPrank, diagonal and degree columns
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

# -----------------------------------------------------------------------------
#' Detect modular structure in a network
#' 
#' \code{modularity}
#' 
#' @param G an adjacency matrix
#' @return a data frame of genes name and module assignment columns
#' @export
modularity <- function(G) {
  MODULARITY_THRESHOLD <- 0.000001
  n <- nrow(G);
  geneNames <- colnames(G)
  
  # create adjacency matrix by thresholding and/or conversion to binary
  # zero the diagonal
  diag(G) <- 0
  # create real symmetric modularity matrix B
  k <- colSums(G)
  m <- 0.5 * sum(k)
  B <- G - k %*% t(k) / (2.0 * m);
  
  # column indices
  firstModule <- seq(from=1, to=n)
  # stack for recursive subdivision of modules
  processStack <- list(firstModule)
  # list of vectors of column indices
  modules <- NULL
  # Q is a measure of "goodness" of the module decomposition
  Q <- 0
  
  # "recursive" loop to find submodules
  iteration <- 0
  while(length(processStack) > 0) {
    iteration <- iteration + 1
    # pop the stack for a list of column indices
    thisModule <- unlist(processStack[[length(processStack)]])
    #cat("this module:", thisModule, "\n")
    processStack[length(processStack)] <- NULL
    # create Bg, a submatrix of B based on current module indices
    newDim <- length(thisModule)
    Bg <- matrix(ncol=newDim, nrow=newDim, data=c(0))
    for(l1 in 1:newDim) {
      for(l2 in 1:newDim) {
        Bg[l1, l2] = B[thisModule[l1], thisModule[l2]];
      }
    }
    # adjust the diagonal
    rowsums <-rowSums(Bg)
    for(i in 1:length(rowsums)) {
      Bg[i, i] = Bg[i, i] - rowsums[i]
    }
    # get the best split of the modules based on eigenvector decomposition
    sub_modules <- modularityBestSplit(Bg, m)
    deltaQ <- sub_modules$Q
    s <- sub_modules$s_out
    # assign indices based on two groups
    s1 <- thisModule[s==-1]
    s2 <- thisModule[s==1]
    if((length(s1) == 0) || (length(s2) == 0)) {
      # stopping criteria for recursive splits
      modules[[length(modules)+1]] <- thisModule
      if(iteration == 1) {
        Q <- deltaQ
      }
    }
    else {
      if(deltaQ <= MODULARITY_THRESHOLD) {
        # stopping criteria for recursive splits
        modules[[length(modules)+1]] <- thisModule
      } else {
        # "recursive step": push the two groups onto the stack and repeat
        processStack[[length(processStack)+1]] <- s1
        processStack[[length(processStack)+1]] <- s2
        # update cummulative Q
        Q <- Q + deltaQ
      }
    }
  }
  
  # ----------------------------------------------------------------------------
  # output modules
  # ----------------------------------------------------------------------------
  cat("Q:", Q, "\n")
  cat("Number of modules:", length(modules), "\n")
  groupAssignments <- NULL
  for(i in 1:length(modules)) {
    modIdx <- modules[[i]]
    modNum <- length(modIdx)
    modGenes <- geneNames[modIdx]
    #cat("Module", i, ":", modGenes, "\n")
    modGroup <- rep(i, modNum)
    thisGroup <- cbind(modGenes, modGroup)
    groupAssignments <- rbind(groupAssignments, thisGroup)
  }
  colnames(groupAssignments) <- c("Gene", "Group")
  groupAssignments
}

# -----------------------------------------------------------------------------
#' Find the best split of a modularity matrix using eigenvector decomposition
#' 
#' \code{modularityBestSplit} 
#' 
#' @param B modularity matrix
#' @param m m
#' @return list with modularity value Q and best split vector
#' @export
modularityBestSplit <- function(B, m) {
  # function to split columns of matrix into two groups
  # get the maximum eigenvector
  eigResult <- eigen(B);
  eigval <- eigResult$values;
  eigvec <- eigResult$vectors
  # in R, this is the first vector
  maxeig_val <- eigval[1];
  maxeig_vec <- eigvec[,1];
  # use the sign of the eigenvector values to assign group status +/-1
  s_out <- ifelse(maxeig_vec < 0, -1, 1)
  # calculate Q for this split
  Q_mat <- t(s_out) %*% B %*% s_out
  Q <- Q_mat[1,1]
  Q <- Q * (1.0 / (m * 4.0))
  
  # return Q and list assignments
  list(Q=Q, s_out=s_out)
}
