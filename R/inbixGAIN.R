# -----------------------------------------------------------------------------
# inbixGAIN.R - Bill White - 10/10/15
#
# Rinbix package genetic association network (GAIN) functions.

# -----------------------------------------------------------------------------
#' Differential coexpression genetic association interaction network algorithm.
#' 
#' \code{dcgain} 
#' 
#' @param inbixData Data frame with samples in rows, genes in columns
#' and phenotype in the last column.
#' @param verbose Flag to send messages to stdout.
#' @return results Matrix of gene by gene differential coexpression values.
#' @family Genetic Association Interaction Network functions
#' @seealso \code{\link{dmgain}} for differential modularity.
#' @examples
#' data(testdata10)
#' rinbixDcgain <- dcgain(testdata10)
#' @export
dcgain <- function(inbixData, verbose=FALSE) {
  phenos <- inbixData[, ncol(inbixData)] + 1
  exprBySubj <- inbixData[, -ncol(inbixData)]
  exprByGene <- t(exprBySubj)
  varNames <- colnames(inbixData)[1:ncol(exprBySubj)]
  n1 <- length(which(phenos == 1))
  n2 <- length(which(phenos == 2))
  if(verbose) {
    print(dim(inbixData))
    cat("dcGAIN Group 1:", n1, "Group 2:", n2, "\n")
  }
  nVars <- ncol(exprBySubj)

  # determine group correlations
  expr_g1 <- exprByGene[, phenos == 1]
  expr_g2 <- exprByGene[, phenos == 2]
#   cat("DEBUG Rinbix package\n")
#   print(expr_g1[1:5, 1:5])
#   print(expr_g2[1:5, 1:5])
  cor_table_g1 <- cor(t(expr_g1))
  cor_table_g2 <- cor(t(expr_g2))
    
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
#' Differential modularity genetic association network algorithm.
#' 
#' \code{dmgain} 
#' 
#' @param inbixData Data frame with samples in rows, genes in columns
#' and phenotype in the last column.
#' @return results Matrix of gene by gene differential modularity values.
#' @examples
#' data(testdata10)
#' rinbixDmgain <- dmgain(testdata10)
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
#' Get interaction effects from generalized linear model regression.
#' 
#' \code{getInteractionEffects} 
#' 
#' @param data Data frame with genes in columns and samples in rows.
#' @param regressionFamily String glm regression family name.
#' @param numCovariates Number of included covariates.
#' @param writeBetas Flag indicating whether to write beta values to separate file.
#' @param excludeMainEffects Flag indicating whether to exclude main effect terms.
#' @param useBetas Flag indicating betas rather than standardized betas used.
#' @param transformMethod String optional transform method.
#' @param verbose Flag to send verbose messages to stdout.
#' @param numCores Number of processor cores to use in mclapply
#' @return results Matrix of gene by gene regression coefficients.
#' @export
getInteractionEffects <- function(data, regressionFamily="binomial", numCovariates=0,
                                  writeBetas=FALSE, excludeMainEffects=FALSE, useBetas=FALSE, 
                                  transformMethod="", verbose=FALSE, numCores=2) {
  allColumnNames <- colnames(data)
  numColumnNames <- length(allColumnNames)
  depVarName <- "Class"
  endVariableNames <- ncol(data) - numCovariates - 1
  variableNames <- allColumnNames[1:endVariableNames]
  numVariables <- length(variableNames)
  if(verbose) {
    cat("From", numVariables ,"generating", choose(numVariables,2), 
        "variable index pairs for reGAIN upper triangular matrix\n")
  }
  idxCombList <- combn(numVariables, 2, list)
  #print(idxCombList)
  if(verbose) { 
    cat("Computing GLM interaction models for each index pair, in parallel\n") 
  }
#  lastIdx <- length(idxCombList)
#  numSplits <- 10
#  splitSize <- as.integer(lastIdx / numSplits) + 1
  results <- NULL
  # for(i in 1:numSplits) {
  #   startIdx <- (i - 1) * splitSize + 1
  #   endIdx <- i *  splitSize
  #   if(endIdx > lastIdx) {
  #     endIdx <- lastIdx
  #   }
  #   if(startIdx < lastIdx) {
  #     # cat("Running chunk", i, "split size:", splitSize, "start:", startIdx, 
  #     #     "end:", endIdx, "\n")
  results <- c(results, 
               parallel::mclapply(idxCombList, mc.cores=numCores,
                                 function(x) runInteractionEffectsTest(data, x, 
                                                                       depVarName, 
                                                                       regressionFamily, 
                                                                       numCovariates, 
                                                                       excludeMainEffects)))
  #   }
  # }
  if(verbose) {
    #print(results)
    cat("Loading the reGAIN matrix upper and lower triangulars with",
        "GLM interaction coefficients\n")
  }
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
      if(verbose) {
        cat("WARNING: Regression model failed to converge for", interactionName, "\n")
      }
    }
    interactionCoeff <- regressionModel$coef
    interactionStdErr <- regressionModel$se
    interactionPval <- regressionModel$pval
    interactionStat <- regressionModel$stdbeta
    
    interactionValue <- NA
#    if(glmConverged) {
      if(useBetas) {
        interactionValue <- interactionCoeff
      }
      else {
        interactionValue <- interactionStat
      }
 #   } else {
  #    interactionValue <- 0
   # }
    if(interactionPval > 0.99) {
      if(verbose) {
        cat("WARNING: Interaction effect p-value > 0.99", interactionName, "\n")
      }
      #interactionValue <- 0
    }
    if(interactionPval < 2e-16) {
      if(verbose) {
        cat("WARNING: Interaction effect p-value < 2e-16", interactionName, "\n")
      }
      #interactionValue <- 0
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

    # -------------------------------------------------------------------------
    interactionValues[variable1Idx, variable2Idx] <- interactionValueTransformed
    interactionValues[variable2Idx, variable1Idx] <- interactionValueTransformed
    # -------------------------------------------------------------------------
    if(verbose) {
      cat("Interaction z value:", variable1Idx, variable2Idx, interactionStat, "\n")
      if(abs(interactionStat) > 4) {
        print(regressionModel)
        print(summary(regressionModel))
      }
    }

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
    if(verbose) cat("Writing betas to intbetas.tab\n") 
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
#' Get main effects from generalized linear model regression (parallel).
#' 
#' \code{getMainEffects} 
#' 
#' @param data Data frame with genes in columns and samples in rows.
#' @param regressionFamily String glm regression family name.
#' @param numCovariates Number of included covariates.
#' @param writeBetas Flag indicating whther to write beta values to separate file.
#' @param useBetas Flag indicating betas rather than standardized betas used.
#' @param transformMethod String optional transform method.
#' @param verbose Flag to send verbose messages to stdout.
#' @param numCores Number of processor cores to use in mclapply.
#' @return mainEffectValues Vector of main effect values.
#' @export
# -----------------------------------------------------------------------------
getMainEffects <- function(data, regressionFamily="binomial", numCovariates=0, 
                           writeBetas=FALSE, useBetas=FALSE, transformMethod="", 
                           verbose=FALSE, numCores=2) {
  allColumnNames <- colnames(data)
  numColumnNames <- length(allColumnNames)
  depVarName <- "Class"
  endVariableNames <- ncol(data) - numCovariates - 1
  variableNames <- allColumnNames[1:endVariableNames]
  numVariables <- length(variableNames)
  if(verbose) {
    cat(paste("Computing GLM main effect models for each variable (", 
              numVariables, "), in parallel", sep=""), "\n")
  }
  results <- parallel::mclapply(variableNames, mc.cores=numCores, 
                      function(x) runMainEffectsTest(data=data, 
                                                     variableName=x, 
                                                     depVarName=depVarName, 
                                                     regressionFamily=regressionFamily, 
                                                     numCovariates=numCovariates))
  if(verbose) {
    cat("Loading the reGAIN matrix diagonal with GLM main effect coefficients\n")
  }
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
      if(verbose) {
        cat("WARNING: Regression model failed to converge for", variableNames[i], "\n")
      }
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
    else {
      mainEffectValue <- 0
    }
    if(mainPval > 0.99) {
      if(verbose) {
        cat("WARNING: Main effect p-value > 0.99", variableNames[i], "\n")
      }
      #mainEffectValue <- 0
    }
    if(mainPval < 2e-16) {
      if(verbose) {
        cat("WARNING: Main effect p-value < 2e-16", variableNames[i], "\n")
      }
      #mainEffectValue <- 0
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
#' Parallel execution of the regression genetic association interaction network algorithm.
#' 
#' \code{regainParallel}
#' 
#' @param regressionData Data frame with gene in columns and samples in rows;
#' the last column should be labeled 'Class' and be 0 or 1 values.
#' @param stdBetas Flag to use standardized beta coefficients.
#' @param absBetas Flag to use absolute value of beta coefficients.
#' @param verbose Flag to send verbose messages to stdout.
#' @return regainMatrix Matrix of gene by gene regression coefficients.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' @export
regainParallel <- function(regressionData, stdBetas=FALSE, absBetas=FALSE, verbose=FALSE) {
  transform <- ifelse(absBetas, "abs", "")
  rawBetas <- ifelse(stdBetas, FALSE, TRUE)
  mainEffects <- getMainEffects(regressionData, 
                                useBetas=rawBetas,
                                transformMethod=transform,
                                verbose=verbose)
  regainMatrix <- getInteractionEffects(regressionData, 
                                        useBetas=rawBetas,
                                        transformMethod=transform,
                                        verbose=verbose)
  diag(regainMatrix) <- mainEffects
  colnames(regainMatrix) <- colnames(regressionData)[1:(ncol(regressionData)-1)]
  # replace NAs with zero
  regainMatrix[is.na(regainMatrix)] <- 0
  regainMatrix
}

# -----------------------------------------------------------------------------
#' Get the interaction effect of a pair of variables using generalized linear regression - glm.
#' 
#' \code{runInteractionEffectsTest} 
#' 
#' @param data Data frame with genes in columns and samples in rows.
#' @param variableIndices Vector of column indices of variable pairs.
#' @param depVarName String name of the phenotype variable.
#' @param regressionFamily String glm regression family name.
#' @param numCovariates Number of covariates included.
#' @param excludeMainEffects Flag indicating whether to exclude main effect terms.
#' @return Data frame with gene, convergence status, beta coefficient,
#' p-value, standard error and standardized beta columns.
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

  # ---------------------------------------------------------------------------
  # use glm2: Fits generalized linear models using the same model specification
  # as glm in the stats package, but with a modified default fitting method that 
  # provides greater stability for models that may fail to converge using glm
  # bcw - 10/18/15
  regressionModel <- glm2::glm2(regressionFormula, family=binomial(link="logit"), data=data)
  # ---------------------------------------------------------------------------
  
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
#' Get the main effect of a variable using generalized linear regression - glm.
#' 
#' \code{runMainEffectsTest}
#' 
#' @param data Data frame with genes in columns and samples in rows.
#' @param variableName String name of the variable to consider.
#' @param depVarName String name of the phenotype variable.
#' @param regressionFamily String glm regression family name.
#' @param numCovariates Number of included covariates.
#' @return Data frame with gene, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns.
#' @export
runMainEffectsTest <- function(data, variableName, depVarName, regressionFamily, numCovariates) {
  if(numCovariates > 0) {
    covarsStart <- ncol(data) - numCovariates
    covarNames <- colnames(data)[covarsStart:(ncol(data)-1)]
    covarModelParts <- NULL
    for(i in 1:length(covarNames)) {
      covarModelParts <- c(covarModelParts, paste("`", covarNames, "`", sep=""))
    }
    formulaString <- paste(depVarName, " ~ `", variableName, 
                           "` + ", paste(covarModelParts, collapse="+"), sep="")
    regressionFormula <- as.formula(formulaString)
  } else {
    regressionFormula <- as.formula(paste(depVarName, "~", 
                                          paste("`", variableName, "`", sep=""), sep=" "))
  }
  #print(formulaString)
  regressionModel <- glm2::glm2(regressionFormula, family=regressionFamily, data=data)
  
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
