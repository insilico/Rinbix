# -----------------------------------------------------------------------------
# inbixGAIN.R - Bill White - 10/10/15
#
# Rinbix package genetic association network (GAIN) functions.

# -----------------------------------------------------------------------------
#' Differential coexpression genetic association interaction network (dcGAIN) algorithm.
#' 
#' \code{dcgain} 
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.biodatamining.org/content/8/1/5}
#'   {Differential co-expression network centrality and machine learning feature selection for identifying susceptibility hubs in networks with scale-free structure}
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @keywords models array
#' @family GAIN functions
#' @family inbix synonym functions
#' @family Genetic Association Interaction Network functions
#' @seealso \code{\link{dmgain}} for differential modularity. 
#' \code{\link{dcgainInbix}} for inbix differential coexpression GAIN.
#' @param inbixData \code{data.frame} with samples in rows, genes in columns
#' and phenotype in the last column.
#' @param verbose \code{logical} to send messages to stdout.
#' @return results \code{matrix} of gene by gene differential coexpression values.
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
#' Differential modularity genetic association network algorithm (dmGAIN) algorithm.
#' 
#' \code{dmgain} 
#' 
#' @references 
#' \itemize{
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @keywords models array
#' @family GAIN functions
#' @family inbix synonym functions
#' \code{\link{dcgain}} for differential correlation GAIN. 
#' @param inbixData \code{data.frame} with samples in rows, genes in columns
#' and phenotype in the last column.
#' @return results \code{matrix} of gene by gene differential modularity values.
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
#' Get the interaction effect of a pair of variables using generalized linear regression.
#' 
#' \code{fitInteractionModel} 
#' 
#' @keywords models regression array
#' @family GAIN functions
#' @param data \code{data.frame} with genes in columns and samples in rows.
#' @param variableIndices \code{vector} of column indices of variable pairs.
#' @param depVarName \code{string} name of the phenotype variable.
#' @param regressionFamily \code{string} glm regression family name.
#' @param numCovariates \code{numeric}  of covariates included.
#' @param excludeMainEffects \code{logical} indicating whether to exclude main effect terms.
#' @return \code{data.frame} with gene, convergence status, beta coefficient,
#' p-value, standard error and standardized beta columns.
#' @keywords internal
fitInteractionModel <- function(data, variableIndices, depVarName, 
                                regressionFamily, numCovariates, excludeMainEffects) {
  variable1Idx <- variableIndices[1]
  variable2Idx <- variableIndices[2]
  variableNames <- colnames(data)[1:(ncol(data)-1)]
  variable1Name <- variableNames[variable1Idx]  
  variable2Name <- variableNames[variable2Idx]
  if(excludeMainEffects) {
    interactionTerm <- paste("`", variable1Name, "`", ":", "`", 
                             variable2Name, "`", sep="")
  } else {
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
#' Get the main effect of a variable using generalized linear regression .
#' 
#' \code{fitMainEffectModel}
#' 
#' @keywords models regression array internal
#' @family GAIN functions
#' @param data \code{data.frame} with genes in columns and samples in rows.
#' @param variableName \code{string} name of the variable to consider.
#' @param depVarName \code{string} name of the phenotype variable.
#' @param regressionFamily \code{string} glm regression family name.
#' @param numCovariates \code{numeric}  of included covariates.
#' @return \code{data.frame} with gene, convergence status, beta coefficient, 
#' p-value, standard error and standardized beta columns.
fitMainEffectModel <- function(data, variableName, depVarName, regressionFamily, 
                               numCovariates) {
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

# -----------------------------------------------------------------------------
#' Get interaction effects from generalized linear model regression.
#' 
#' \code{getInteractionEffects} 
#' 
#' @keywords models regression array internal
#' @family GAIN functions
#' @param data \code{data.frame} with genes in columns and samples in rows.
#' @param regressionFamily \code{string} glm regression family name.
#' @param numCovariates \code{numeric}  of included covariates.
#' @param writeBetas \code{logical} indicating whether to write beta values to separate file.
#' @param excludeMainEffects \code{logical} indicating whether to exclude main effect terms.
#' @param useBetas \code{logical} indicating betas rather than standardized betas used.
#' @param transformMethod \code{string} optional transform method.
#' @param numCores \code{numeric} number of processor cores to use in mclapply
#' @param verbose \code{logical} to send verbose messages to stdout.
#' @return results \code{matrix} of gene by gene regression coefficients.
getInteractionEffects <- function(data, regressionFamily="binomial", numCovariates=0,
                                  writeBetas=FALSE, excludeMainEffects=FALSE, useBetas=FALSE, 
                                  transformMethod="", numCores=2, verbose=FALSE) {
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
  if(verbose) { 
    cat("Computing GLM interaction models for each index pair, in parallel\n") 
  }
  #startTime <- proc.time()
  results <- parallel::mclapply(idxCombList, 
                                mc.cores=numCores,
                                function(x) 
                                  fitInteractionModel(data, x, 
                                                      depVarName, 
                                                      regressionFamily, 
                                                      numCovariates, 
                                                      excludeMainEffects))
  #endTime <- proc.time()
  #cat("getInteractionEffects mclapply GLMFIT runtime:", (endTime - startTime)[3], "\n")
  if(verbose) {
    cat("Loading the reGAIN matrix upper and lower triangulars with",
        "GLM interaction coefficients\n")
  }
  if(writeBetas) {
    betaInfo <- NULL
    betaRows <- NULL
  }
  #startTime <- proc.time()
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
    if(useBetas) {
      interactionValue <- interactionCoeff
    } else {
      interactionValue <- interactionStat
    }
    if(interactionPval > 0.99) {
      if(verbose) {
        cat("WARNING: Interaction effect p-value > 0.99", interactionName, "\n")
      }
    }
    if(interactionPval < 2e-16) {
      if(verbose) {
        cat("WARNING: Interaction effect p-value < 2e-16", interactionName, "\n")
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
        } else {
          varPval <- summary(regressionModel)$coefficients[colIdx, "Pr(>|t|)"]
        }
        allBetas <- c(allBetas, varCoef, varPval)
      }
      betaInfo <- rbind(betaInfo, allBetas)
      betaRows <- c(betaRows, interactionName)
    }
  }
  #endTime <- proc.time()
  #cat("getInteractionEffects mclapply PROCESSING runtime:", (endTime - startTime)[3], "\n")
  
  if(writeBetas) {
    if(verbose) cat("Writing interaction betas to intbetas.tab\n") 
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
#' @keywords models regression univar array internal
#' @family GAIN functions
#' @family inbix synonym functions
#' @seealso \code{\link{rankUnivariateRegression}}
#' @param data \code{data.frame} with genes in columns and samples in rows.
#' @param regressionFamily \code{string} glm regression family name.
#' @param numCovariates \code{numeric}  of included covariates.
#' @param writeBetas \code{logical} indicating whther to write beta values to separate file.
#' @param useBetas \code{logical} indicating betas rather than standardized betas used.
#' @param transformMethod \code{string} optional transform method.
#' @param numCores \code{numeric} number of processor cores to use in mclapply
#' @param verbose \code{logical} to send verbose messages to stdout.
#' @return mainEffectValues \code{vector} of main effect values.
# -----------------------------------------------------------------------------
getMainEffects <- function(data, regressionFamily="binomial", numCovariates=0, 
                           writeBetas=FALSE, useBetas=FALSE, transformMethod="", 
                           numCores=2, verbose=FALSE) {
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
  #startTime <- proc.time()
  results <- parallel::mclapply(variableNames, 
                                mc.cores=numCores, 
                                function(x) 
                                  fitMainEffectModel(data=data, 
                                                     variableName=x, 
                                                     depVarName=depVarName, 
                                                     regressionFamily=regressionFamily, 
                                                     numCovariates=numCovariates))
  #endTime <- proc.time()
  #cat("getMainEffects mclapply GLMFIT runtime:", (endTime - startTime)[3], "\n")
  if(verbose) {
    cat("Loading the reGAIN matrix diagonal with GLM main effect coefficients\n")
  }
  if(writeBetas) {
    betaInfo <- NULL
  }
  #startTime <- proc.time()
  mainEffectValues <- vector(mode="numeric", length=numVariables)
  for(i in 1:length(results)) {
    regressionModel <- results[[i]]
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
    }
    if(mainPval < 2e-16) {
      if(verbose) {
        cat("WARNING: Main effect p-value < 2e-16", variableNames[i], "\n")
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
    mainEffectValues[i] <- mainEffectValueTransformed
    if(writeBetas) {
      thisBetas <- c(mainCoeff, mainPval)
      if(numCovariates > 0) {
        for(cn in 1:numCovariates) {
          covBeta <- summary(regressionModel)$coefficients[2 + cn, "Estimate"]
          covPval <- summary(regressionModel)$coefficients[2 + cn, "Pr(>|z|)"]
          thisBetas <- c(thisBetas, covBeta, covPval)
        }
      }
      betaInfo <- rbind(betaInfo, thisBetas)
    }
  }
  #endTime <- proc.time()
  #cat("getMainEffects mclapply PROCESSING runtime:", (endTime - startTime)[3], "\n")
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
#' Front end to the reGAIN algorithm.
#' 
#' Tests number of variables/columns and runs either parallel R or parallel C++.
#' 
#' \code{regain}
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.nature.com/gene/journal/v11/n8/full/gene201037a.html}
#'   {Genes & Immunity Paper}
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @keywords models regression array
#' @family GAIN functions
#' @family inbix synonym functions
#' @param regressionData \code{data.frame} with genes in columns and samples in rows;
#' the last column should be labeled 'Class' coded 0 or 1.
#' @param stdBetas \code{logical} to use standardized beta coefficients.
#' @param absBetas \code{logical} to use absolute value of beta coefficients.
#' @param numCores \code{numeric} number of processor cores to use in mclapply
#' @param verbose \code{logical} to send verbose messages to stdout.
#' @return regainMatrix \code{matrix} of gene by gene regression coefficients.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regain(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' @export
regain <- function(regressionData, stdBetas=FALSE, absBetas=FALSE, numCores=2, 
                   verbose=FALSE) {
  numVars <- ncol(regressionData) - 1  
  regainMatrix <- NULL
  if(numVars <= 5000) {
    if(verbose) cat("[", numVars, "] Running parallel R reGAIN\n")
    regainMatrix <- regainParallel(regressionData, stdBetas, absBetas,
                                   numCores=numCores, verbose)
  } else {
    cat("WARNING: More than 5000 genes, attempting to run reGAIN in C++\n")
    if(verbose) cat("[", numVars, "] Running C++ inbix reGAIN\n")
    regainMatrix <- regainInbix(regressionData, stdBetas, absBetas)$reGAIN
  }
  regainMatrix
}

# -----------------------------------------------------------------------------
#' Parallel execution of the regression genetic association interaction network (reGAIN) algorithm.
#' 
#' \code{regainParallel}
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.nature.com/gene/journal/v11/n8/full/gene201037a.html}
#'   {Genes & Immunity Paper}
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @keywords models array regression
#' @family GAIN functions
#' @family inbix synonym functions
#' @param regressionData \code{data.frame} with gene in columns and samples in rows;
#' the last column should be labeled 'Class' and be 0 or 1 values.
#' @param stdBetas \code{logical} to use standardized beta coefficients.
#' @param absBetas \code{logical} to use absolute value of beta coefficients.
#' @param numCores \code{numeric} number of processor cores to use in mclapply
#' @param verbose \code{logical} to send verbose messages to stdout.
#' @return regainMatrix \code{matrix} of gene by gene regression coefficients.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' @export
regainParallel <- function(regressionData, stdBetas=FALSE, absBetas=FALSE, 
                           numCores=2, verbose=FALSE) {
  transform <- ifelse(absBetas, "abs", "")
  rawBetas <- ifelse(stdBetas, FALSE, TRUE)
  mainEffects <- getMainEffects(regressionData, 
                                useBetas=rawBetas,
                                transformMethod=transform,
                                numCores=numCores,
                                verbose=verbose)
  regainMatrix <- getInteractionEffects(regressionData, 
                                        useBetas=rawBetas,
                                        transformMethod=transform,
                                        numCores=numCores,
                                        verbose=verbose)
  diag(regainMatrix) <- mainEffects
  colnames(regainMatrix) <- colnames(regressionData)[1:(ncol(regressionData)-1)]
  # replace NAs with zero
  regainMatrix[is.na(regainMatrix)] <- 0
  regainMatrix
}
