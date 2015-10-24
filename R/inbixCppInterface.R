# ----------------------------------------------------------------------------
# inbixCppInterface.R - Bill White - 11/26/14
#
# Rinbix package functions to call C++ inbix.

# ----------------------------------------------------------------------------
#' Differential co-expression genetic association network - dcGAIN.
#' 
#' \code{dcgainInbix} 
#' 
#' @param regressionData Data frame with sample rows, gene columns 
#' plus phenotype column.
#' @param outPrefix String file output prefix.
#' @return List of gene scores and associated p-values.
#' @examples
#' data(testdata10)
#' rinbixCppDcgain <- dcgainInbix(testdata10)
#' @export
dcgainInbix <- function(regressionData, outPrefix="Rinbix") {
  inbixExists()
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, "Rinbix")
  # run inbix reGAIN
  inbixCmd <- paste("inbix --dcgain --dcgain-abs --numeric-file Rinbix.num --pheno Rinbix.pheno --1 --out", 
                    outPrefix)
  # inbixCmd <- paste("inbix --dcgain --numeric-file Rinbix.num --pheno Rinbix.pheno --1 --out", 
  #                   outPrefix)
  #cat("Running inbix command:", inbixCmd, "\n")
  regainStdout <- system(inbixCmd, intern=T)
  inbixDcgain <- read.table("Rinbix.dcgain", header=TRUE, sep="\t")
  inbixDcgainPvals <- read.table("Rinbix.pvals.dcgain", header=TRUE, sep="\t")
  
  # remove temporary files
  file.remove(c("Rinbix.dcgain", "Rinbix.pvals.dcgain", "Rinbix.pheno", "Rinbix.num", "Rinbix.log"))
  
  # return regain matrix
  list(scores=as.matrix(inbixDcgain), pvals=as.matrix(inbixDcgainPvals))
}

# ----------------------------------------------------------------------------
#' Check for the inbix executable in the PATH.
#' 
#' \code{inbixExists} 
#' 
#' @return TRUE if found, else exits with error
#' @export
inbixExists <- function() {
  if(Sys.which("inbix") == "") {
    stop("inbix is not in the PATH")
  }
  TRUE
}

# ----------------------------------------------------------------------------
#' Network module/community detection.
#' 
#' \code{modularityInbix} 
#' 
#' @param gainMatrix GAIN matrix.
#' @param outPrefix String file output prefix.
#' @return list of gene module assignments.
#' @examples
#' data(testdata10)
#' corMatrix <- cor(testdata10[, -ncol(testdata10)])
#' inbixModulesDF <- modularityInbix(corMatrix)
#' @export
modularityInbix <- function(gainMatrix, outPrefix="Rinbix") {
  inbixExists()
  # write gainMatrix to inbix-compatible file
  write.table(gainMatrix, file="Rinbix.gain", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
  inbixCmd <- paste("inbix --regain-file Rinbix.gain --modularity --out", outPrefix)
  #cat("Running inbix command:", inbixCmd, "\n")
  # run inbix
  modStdout <- system(inbixCmd, intern=T)
  #print(snprankStdout)
  # read results
  inbixModules <- read.table("Rinbix.modules", header=FALSE, sep="\t")
  colnames(inbixModules) <- c("Var", "Module")
  # remove temporary files
  file.remove(c("Rinbix.gain", "Rinbix.modules", "Rinbix.log"))
  # return module assignments
  inbixModules[order(inbixModules[,1]), ]
}

# ----------------------------------------------------------------------------
#' Permute a GAIN method.
#' 
#' \code{permuteGainInbix} 
#' 
#' @param regressionData Data frame with sample rows, gene columns
#' plus phenotype column.
#' @param method String regain or dcgain.
#' @param numPerms Numeric number of permutation runs.
#' @param pThresh Numeric p-value threshold for GAIN method.
#' @param threshold Numeric permute GAIN threshold.
#' @param outPrefix String file output prefix.
#' @return list of gene thresholds.
#' @export
permuteGainInbix <- function(regressionData, method="regain", numPerms=100, 
  pThresh=1, threshold=0.05, outPrefix="Rinbix") {
  inbixExists()
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, "Rinbix")
  # run inbix GAIN
  if(method == "regain") {
    methodOptions = "--regain-use-beta-values --regain-matrix-transform abs"
    #methodOptions <- "--regain-matrix-transform abs"
    #methodOptions <- ""
  }
  if(method == "dcgain") {
    #methodOptions <- ""
    methodOptions <- "--dcgain-abs"
  }
  methodOption <- paste("--permute-gain-method", method, methodOptions)
  numOption <- paste("--permute-gain-num", numPerms)
  thresholdOption <- paste("--permute-gain-threshold", threshold)
  pThresholdOption <- ""
  if(pThresh < 1) {
    pThresholdOption <- paste("--regain-pvalue-threshold", pThresh)
  }
  inbixCmd <- paste("inbix", methodOption, numOption, thresholdOption, pThresholdOption,
    "--numeric-file Rinbix.num --pheno Rinbix.pheno --1 --out", outPrefix)
  #cat("Running inbix command:", inbixCmd, "\n")
  regainStdout <- system(inbixCmd, intern=T)
  #print(regainStdout)
  inbixThresholds <- read.table("Rinbix_thresholds.txt", header=FALSE, sep="\t")
  
  # remove temporary files
  if(method == "regain") {
    file.remove(c("Rinbix_thresholds.txt", "Rinbix.num", "Rinbix.pheno", "Rinbix.block.sif", 
      "Rinbix.block.mebetas", "Rinbix.block.betas", "Rinbix.log", "Rinbix.perm"))
  }  
  if(method == "dcgain") {
    file.remove(c("Rinbix_thresholds.txt", "Rinbix.pheno", "Rinbix.num", "Rinbix.log", "Rinbix.perm"))
  }
  
  # return gene thresholds
  inbixThresholds
}

# ----------------------------------------------------------------------------
#' Read inbix numeric data set and phenotype file; return combined data frame.
#' 
#' \code{readInbixNumericAsRegressionData} 
#' 
#' @param baseInbixName String base filename.
#' @return data frame with numeric data in first m columns and phenotye in m+1 column.
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
#' Regression genetic association network - reGAIN.
#' 
#' \code{regainInbix} 
#' 
#' @param regressionData Data frame with sample rows, gene columns 
#' plus phenotype column.
#' @param stdBetas Flag use standardized beta coefficients.
#' @param absBetas Flag take absolute value of beta coefficients.
#' @param outPrefix String file output prefix.
#' @param pThreshold Numeric p-value threshold for GAIN method.
#' @return list of gene scores and associated p-values.
#' @examples
#' data(testdata10)
#' inbixRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' @export
regainInbix <- function(regressionData, stdBetas=TRUE, absBetas=TRUE, 
                        outPrefix="Rinbix", pThreshold=1) {
  inbixExists()
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, outPrefix)
  stdBetasCmd <- ""
  if(!stdBetas) {
    stdBetasCmd <- "--regain-use-beta-values"
  }
  absBetasCmd <- ""
  if(absBetas) {
    absBetasCmd <- "--regain-matrix-transform abs"  
  }
  pThresholdCmd <- ""
  if(pThreshold < 1) {
    pThresholdCmd <- paste("--regain-pvalue-threshold", pThreshold)
  }
  # run inbix reGAIN
  inbixCmd <- paste("inbix --regain --numeric-file Rinbix.num --pheno Rinbix.pheno --1 --out", 
                    outPrefix, stdBetasCmd, absBetasCmd, pThresholdCmd)
  #cat("Running inbix command:", inbixCmd, "\n")
  regainStdout <- system(inbixCmd, intern=T)
  inbixRegain <- read.table("Rinbix.block.regain", header=TRUE, sep="\t")
  
  # read warnings and errors if they exist
  warningsText <- ""
  failuresText <- ""
  if(file.exists("Rinbix.regression.warnings")) {
    warningsText <- readLines("Rinbix.regression.warnings")
    #cat(warningsText, "\n")
    file.remove("Rinbix.regression.warnings")
  }
  if(file.exists("Rinbix.regression.failures")) {
    failuresText <- readLines("Rinbix.regression.failures")
    #cat(failuresText, "\n")
    file.remove("Rinbix.regression.failures")
  }

  #file.copy(from="Rinbix.block.betas", to="foo.betas")

  # remove temporary files
  file.remove(c("Rinbix.block.regain", "Rinbix.num", "Rinbix.pheno", "Rinbix.block.sif", 
                "Rinbix.block.pvals.regain", "Rinbix.block.mebetas", "Rinbix.block.betas",
                "Rinbix.log"))
  
  # return regain matrix
  list(reGAIN=as.matrix(inbixRegain), warningsText=warningsText, failuresText=failuresText)
}

# ----------------------------------------------------------------------------
#' SNPrank network centrality gene ranker based on GAIN matrix.
#' 
#' \code{snprankInbix} 
#' 
#' @param gainMatrix GAIN matrix.
#' @param outPrefix String file output prefix.
#' @param gamma Numeric gamma damping parameter (see paper).
#' @return data frame of gene snpranks.
#' @examples
#' data(testdata10)
#' inbixRegain <- regainInbix(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' inbixSnprank <- snprankInbix(inbixRegain)
#' @export
snprankInbix <- function(gainMatrix, outPrefix="Rinbix", gamma=0.85) {
  inbixExists()
  # write gainMatrix to inbix-compatible file
  write.table(gainMatrix, file="Rinbix.gain", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
  gammaCmd <- ""
  if(gamma > 0) {
    gammaCmd <- paste("--rank-centrality-gamma", gamma)
  }
  # run inbix SNPrank
  #cat("gamma:", gamma, "gamma command:", gammaCmd, "\n")
  inbixCmd <- paste("inbix --regain-file Rinbix.gain --rank-by centrality --out", 
                    outPrefix, gammaCmd)
  #cat("Running inbix command:", inbixCmd, "\n")
  snprankStdout <- system(inbixCmd, intern=T)
  #print(snprankStdout)
  inbixSnpranks <- read.table("Rinbix.ranks", header=TRUE, sep="\t")
  
  # remove temporary files
  file.remove(c("Rinbix.gain", "Rinbix.ranks", "Rinbix.log"))
  
  # return snpranks
  data.frame(gene=inbixSnpranks$SNP, snprank=inbixSnpranks$SNPrank)
}

# ----------------------------------------------------------------------------
#' Write regression data frame as  inbix numeric data set and phenotype file.
#' 
#' \code{writeRegressionDataAsInbixNumeric} 
#' 
#' @param regressionData Data frame with sample rows, gene columns 
#' plus phenotype column.
#' @param baseInbixName String base filename to write.
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
