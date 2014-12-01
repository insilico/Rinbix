# inbixCppInterface.R - Bill White - 11/26/14
#
# Functions to support call C++ inbix library including:
# -reading inbix numeric format data and phenotypes
# -reading results of calling command-line inbix through system() calls:
#   -permuteGainInbix
#   -dcgainInbix
#   -regainInbix
#   -snprankInbix
#   -modularityInbix

# ----------------------------------------------------------------------------
#' Read inbix numeric data set and phenotypes
#' 
#' \code{readInbixNumericPheno} 
#' 
#' @param inbixNumericFilename, inbix numeric filename
#' @param inbixPhenoFilename, inbix phenotype filename
#' @return list with numeric matrix and phenotypes
#' @export
readInbixNumericPheno <- function(inbixNumericFilename, inbixPhenoFilename) {
  inbixNumericTable <- read.table(inbixNumericFilename, header=T, sep="\t")
  inbixPhenoTable <- read.table(inbixPhenoFilename, header=F, sep="\t")
  list(inbixNumeric=inbixNumericTable, inbixPheno=inbixPhenoTable)
}

# ----------------------------------------------------------------------------
#' Write regression data set in inbix format.
#' 
#' \code{writeRegressionInbixDataset}
#' 
#' @param regDs is a data frame with sample rows, gene columns 
#' plus phenotype column
#' @param filePrefix, filename prefix for numeric and phenotype files
#' @export
writeRegressionInbixDataset <- function(regDs, filePrefix) {
  numGenes <- ncol(regDs) - 1
  numSubjs <- nrow(regDs)
  numG1 <- numSubjs / 2
  numG2 <- numG1
  phenos <- regDs[, numGenes + 1]
  subIds <- c(paste("case", 1:numG1, sep=""), paste("ctrl", 1:numG2, sep=""))
  
  phenosTable <- cbind(subIds, subIds, phenos)
  datasimInbixPhenoFile <- paste(filePrefix, ".pheno", sep="")
  write.table(phenosTable, datasimInbixPhenoFile, quote=F, sep="\t",
              col.names=F, row.names=F)
  
  dataTable <- cbind(subIds, subIds, regDs[, 1:numGenes])
  #geneNames <- paste("gene", sprintf("%04d", 1:numGenes), sep="")
  geneNames <- colnames(regDs)[1:numGenes]
  colnames(dataTable) <- c("FID", "IID", geneNames)
  datasimInbixNumFile <- paste(filePrefix, ".num", sep="")
  write.table(dataTable, datasimInbixNumFile, quote=F, sep="\t",
              col.names=T, row.names=F)
}

# ----------------------------------------------------------------------------
#' Permute a GAIN method.
#' 
#' \code{permuteGainInbix} 
#' 
#' @param regressionData is a data frame with sample rows, gene columns 
#' plus phenotype column
#' @param method, regain or dcgain
#' @param numPerms, number of permutation runs
#' @param pThresh, p-value threshold for GAIN method
#' @param threshold, permute GAIN threshold
#' @param outPrefix, file output prefix
#' @return list of gene thresholds
#' @export
permuteGainInbix <- function(regressionData, method="regain", numPerms=100, 
  pThresh=1, threshold=0.05, outPrefix="Rinbix") {
  # write regressionData data frame to inbix files
  writeRegressionInbixDataset(regressionData, "Rinbix")
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
#' Differential co-expression genetic association network - dcGAIN
#' 
#' \code{dcgainInbix} 
#' 
#' @param regressionData is a data frame with sample rows, gene columns 
#' plus phenotype column
#' @param outPrefix, file output prefix
#' @return list of gene scores and associated p-values
#' @export
dcgainInbix <- function(regressionData, outPrefix="Rinbix") {
  # write regressionData data frame to inbix files
  writeRegressionInbixDataset(regressionData, "Rinbix")
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
#' Regression genetic association network - reGAIN.
#' 
#' \code{regainInbix} 
#' 
#' @param regressionData is a data frame with sample rows, gene columns 
#' plus phenotype column
#' @param stdBetas, use standardized beta coefficients
#' @param absBetas, take absolute value of beta coefficients
#' @param outPrefix, file output prefix
#' @param pThreshold, p-value threshold for GAIN method
#' @return list of gene scores and associated p-values
#' @export
regainInbix <- function(regressionData, stdBetas=TRUE, absBetas=TRUE, 
                        outPrefix="Rinbix", pThreshold=1) {
  # write regressionData data frame to inbix files
  writeRegressionInbixDataset(regressionData, outPrefix)
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
#' @param gainMatrix, GAIN matrix
#' @param outPrefix, file output prefix
#' @param gamma, gamma parameter (see paper)
#' @return data frame of gene snpranks
#' @export
snprankInbix <- function(gainMatrix, outPrefix="Rinbix", gamma=0.85) {
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
#' Network module/community detection.
#' 
#' \code{modularityInbix} 
#' 
#' @param gainMatrix, GAIN matrix
#' @param outPrefix, file output prefix
#' @return list of gene module assignments
#' @export
modularityInbix <- function(gainMatrix, outPrefix="Rinbix") {
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
