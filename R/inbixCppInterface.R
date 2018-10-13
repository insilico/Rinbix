# ----------------------------------------------------------------------------
# inbixCppInterface.R - Bill White - 11/26/14
#
# Rinbix package functions to call C++ inbix.
#
# TODO: Use Rcpp to implement the algorithm to avoid the standalone inbix
# compile and install. bcw 7/6/17

#' Differential co-expression genetic association network - dcGAIN.
#' 
#' \code{dcgainInbix} 
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.biodatamining.org/content/8/1/5}
#'   {Differential co-expression network centrality and machine learning feature selection for identifying susceptibility hubs in networks with scale-free structure}
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @family inbix interface functions
#' @family GAIN functions
#' @seealso \code{\link{dcgain}}
#' @param regressionData \code{data.frame} with sample rows, gene columns 
#' plus phenotype column.
#' @param outPrefix \code{string} file output prefix.
#' @return List of gene scores and associated p-values.
#' @examples
#' data(testdata10)
#' rinbixCppDcgain <- dcgainInbix(testdata10)
#' @export
dcgainInbix <- function(regressionData, outPrefix = "Rinbix", verbose = FALSE) {
  inbixExists()
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, outPrefix)
  # run inbix reGAIN
  inbixCmd <- paste("inbix --dcgain --dcgain-abs --numeric-file Rinbix.num --pheno Rinbix.pheno --1 --out", 
                    outPrefix)
  if(verbose) cat("Running inbix command:", inbixCmd, "\n")
  system(inbixCmd, intern = TRUE)
  inbixDcgain <- read.table("Rinbix.dcgain.tab", header = TRUE, sep = "\t")
  inbixDcgainPvals <- read.table("Rinbix.pvals.dcgain.tab", header = TRUE, sep = "\t")
  
  # remove temporary files
  file.remove(c("Rinbix.dcgain.tab", "Rinbix.pvals.dcgain.tab", "Rinbix.pheno", "Rinbix.num", "Rinbix.log"))
  
  # return regain matrix
  list(scores = as.matrix(inbixDcgain), pvals = as.matrix(inbixDcgainPvals))
}

#' Check for the inbix executable in the PATH.
#' 
#' \code{inbixExists} 
#' 
#' @return \code{logical} TRUE if found, else exits with error
#' @export
inbixExists <- function() {
  if (Sys.which("inbix") == "") {
    stop("inbix is not in the PATH")
  }
  TRUE
}

#' Network module/community detection.
#' 
#' \code{modularityInbix} 
#' 
#' @references 
#' \itemize{
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @family inbix interface functions
#' @seealso \code{\link{modularity}}
#' @param gainMatrix \code{matrix} GAIN matrix.
#' @param outPrefix \code{string} file output prefix.
#' @return \code{list} of gene module assignments.
#' @examples
#' data(testdata10)
#' corMatrix <- cor(testdata10[, -ncol(testdata10)])
#' inbixModulesDF <- modularityInbix(corMatrix)
#' @export
modularityInbix <- function(gainMatrix, outPrefix = "Rinbix") {
  inbixExists()
  # write gainMatrix to inbix-compatible file
  write.table(gainMatrix, file = "Rinbix.gain", quote = FALSE, 
              col.names = TRUE, row.names = FALSE, sep = "\t")
  inbixCmd <- paste("inbix --regain-file Rinbix.gain --modularity --out", outPrefix)
  #cat("Running inbix command:", inbixCmd, "\n")
  # run inbix
  system(inbixCmd, intern = TRUE)
  #print(snprankStdout)
  # read results
  inbixModules <- read.table("Rinbix.modules", header = FALSE, sep = "\t")
  colnames(inbixModules) <- c("Var", "Module")
  # remove temporary files
  file.remove(c("Rinbix.gain", "Rinbix.modules", "Rinbix.log"))
  # return module assignments
  inbixModules[order(inbixModules[,1]), ]
}

#' Permute a GAIN method.
#' 
#' \code{permuteGainInbix} 
#' 
#' @references 
#' \itemize{
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @family inbix interface functions
#' @param regressionData \code{data.frame} with sample rows, gene columns
#' plus phenotype column.
#' @param method \code{string} regain or dcgain.
#' @param numPerms \code{numeric} number of permutation runs.
#' @param pThresh \code{numeric} p-value threshold for GAIN method.
#' @param threshold \code{numeric} permute GAIN threshold.
#' @param outPrefix \code{string} file output prefix.
#' @return \code{list} of gene thresholds.
#' @examples 
#' data(testdata10)
#' inbixRegainThresholds <- permuteGainInbix(testdata10)
#' @export
permuteGainInbix <- function(regressionData, method = "regain", numPerms = 100, 
                             pThresh = 1, threshold = 0.05, outPrefix = "Rinbix") {
  inbixExists()
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, "Rinbix")
  # run inbix GAIN
  if (method == "regain") {
    methodOptions <- "--regain-use-beta-values --regain-matrix-transform abs"
    #methodOptions <- "--regain-matrix-transform abs"
    #methodOptions <- ""
  }
  if (method == "dcgain") {
    #methodOptions <- ""
    methodOptions <- "--dcgain-abs"
  }
  methodOption <- paste("--permute-gain-method", method, methodOptions)
  numOption <- paste("--permute-gain-num", numPerms)
  thresholdOption <- paste("--permute-gain-threshold", threshold)
  pThresholdOption <- ""
  if (pThresh < 1) {
    pThresholdOption <- paste("--regain-pvalue-threshold", pThresh)
  }
  inbixCmd <- paste("inbix", methodOption, numOption, thresholdOption, pThresholdOption,
    "--numeric-file Rinbix.num --pheno Rinbix.pheno --1 --out", outPrefix)
  #cat("Running inbix command:", inbixCmd, "\n")
  system(inbixCmd, intern = TRUE)
  #print(regainStdout)
  inbixThresholds <- read.table("Rinbix_thresholds.txt", header = FALSE, sep = "\t")
  
  # remove temporary files
  if (method == "regain") {
    file.remove(c("Rinbix_thresholds.txt", "Rinbix.num", "Rinbix.pheno", "Rinbix.block.sif", 
      "Rinbix.block.mebetas", "Rinbix.block.betas", "Rinbix.log", "Rinbix.perm"))
  }  
  if (method == "dcgain") {
    file.remove(c("Rinbix_thresholds.txt", "Rinbix.pheno", "Rinbix.num", "Rinbix.log", "Rinbix.perm"))
  }
  
  # return gene thresholds
  inbixThresholds
}

#' Call C++ inbix Evaporative Cooling Privacy
#'
#' Assumes the inbix executable is in the PATH.
#'
#' @param train.ds A data frame with training data and class labels
#' @param holdout.ds A data frame with holdout data and class labels
#' @param validation.ds A data frame with validation data and class labels
#' @param label A character vector of the class variable column name
#' @param is.simulated Is the data simulated (or real?)
#' @param bias A numeric for effect size in simulated signal variables
#' @param update.freq A integer for the number of steps before update
#' @param start.temp A numeric for EC starting temperature
#' @param final.temp A numeric for EC final temperature
#' @param tau.param A numeric for tau to control reduction schedule
#' @param rf.ntree An integer the number of trees in the random forest
#' @param rf.mtry An integer the number of variables sampled at each random forest node split
#' @param save.file A character vector for the results filename or NULL to skip
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @note inbix must be instaled in the path
#' @return A list containing:
#' \describe{
#'   \item{algo.acc}{data frame of results, a row for each update}
#'   \item{ggplot.data}{melted results data frame for plotting}
#'   \item{correct}{number of variables detected correctly in each data set}
#'   \item{elapsed}{total elapsed time}
#' }
#' @examples
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' sim.data <- privateEC::createSimulation(num.samples = num.samples,
#'                                         num.variables = num.variables,
#'                                         pct.train = 1 / 3,
#'                                         pct.holdout = 1 / 3,
#'                                         pct.validation = 1 / 3,
#'                                         pct.signals = pct.signals,
#'                                         sim.type = "mainEffect",
#'                                         verbose = FALSE)
#' pec.results <- privateECinbix(train.ds = sim.data$train,
#'                               holdout.ds = sim.data$holdout,
#'                               validation.ds = sim.data$validation,
#'                               verbose = FALSE)
#' @family classification
#' @export
privateECinbix <- function(train.ds=NULL, holdout.ds=NULL, validation.ds=NULL,
                           label="phenos",
                           is.simulated=TRUE,
                           bias=0.4,
                           update.freq=50,
                           start.temp=0.1,
                           final.temp=10 ^ (-5),
                           tau.param=100,
                           rf.ntree=500,
                           rf.mtry=NULL,
                           save.file=NULL,
                           verbose=FALSE) {
  # check for inbix in the PATH or stop with an error
  if (Sys.which("inbix") == "") {
    stop("inbix is not installed or not in the PATH")
  }
  # check input parameters
  train.good <- is.null(train.ds)
  holdout.good <- is.null(holdout.ds)
  validation.good <- is.null(validation.ds)
  sum.good <- sum(train.good, holdout.good, validation.good)
  if (sum.good >= 2) {
    stop("At least train and holdout data sets must be provided")
  }
  d <- ncol(train.ds) - 1
  param.mtry <- rf.mtry
  if (is.null(rf.mtry)) {
    param.mtry <- floor(sqrt(d))
  }
  if (param.mtry < 1 | param.mtry > d) {
    stop(paste("mtry parameter", param.mtry, "out of range 1 < mtry < d"))
  }
  ptm <- proc.time()
  unique.sim.prefix <- tempfile("pec")
  correct.detect.inbix <- vector(mode = "numeric")
  # write simple tab-separated vales (tsv) files
  if (verbose) cat("Writing ", unique.sim.prefix,
                   ".sim.(train|holdout|validation).tab files\n")
  utils::write.table(train.ds, paste(unique.sim.prefix, ".sim.train.tab",
                                     sep = ""), sep = "\t",
                     row.names = FALSE, col.names = TRUE, quote = FALSE)
  utils::write.table(holdout.ds, paste(unique.sim.prefix, ".sim.holdout.tab",
                                       sep = ""), sep = "\t",
                     row.names = FALSE, col.names = TRUE, quote = FALSE)
  if (sum.good > 2) {
    utils::write.table(validation.ds, paste(unique.sim.prefix, ".sim.test.tab",
                                            sep = ""), sep = "\t",
                       row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  if (verbose) cat("Running inbix privacy EC on saved simulation data sets\n")
  system.command <- paste("OMP_NUM_THREADS=4 inbix ",
                          "--ec-privacy ",
                          "--depvarname ", label, " ",
                          "--ec-privacy-train-file ", unique.sim.prefix, ".sim.train.tab ",
                          "--ec-privacy-holdout-file ", unique.sim.prefix, ".sim.holdout.tab ",
                          "--ec-privacy-start-temp ", start.temp, " ",
                          "--ec-privacy-final-temp ", final.temp, " ",
                          "--ec-privacy-tau ", tau.param, " ",
                          "--ec-privacy-update-frequency ", update.freq, " ",
                          "--ntree ", rf.ntree, " ",
                          "--mtry ", param.mtry, " ",
                          "--out ", unique.sim.prefix, ".privateec", sep = "")
  if (sum.good > 2) {
    c(system.command,
      paste("--ec-privacy-test-file ", unique.sim.prefix, ".sim.test.tab ",
            sep = ""))
  }
  if (verbose) cat(system.command, "\n")
  system(system.command, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  # if(verbose) cat("Reading inbix privacy EC attributes used results\n")
  # attrs.file <- paste(unique.sim.prefix, ".privateec.selected.attributes.tab", sep="")
  # if(verbose) cat(attrs.file, "\n")
  # attrs.table <- utils::read.table(attrs.file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  
  if (verbose) cat("Reading inbix privacy EC algorithm run details\n")
  iters.file <- paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep = "")
  if (verbose) cat(iters.file, "\n")
  iters.table <- utils::read.table(iters.file, sep = "\t", header = TRUE,
                                   stringsAsFactors = FALSE)
  if (is.simulated) {
    correct.detect.inbix <- as.integer(iters.table$Correct)
  }
  fxplots <- data.frame(vars.remain = as.integer(iters.table$Keep),
                        ftrain = iters.table$TrainAcc,
                        fholdout = iters.table$HoldoutAcc,
                        fvalidation = iters.table$TestAcc,
                        alg = 5)
  melted.fx <- reshape2::melt(fxplots, id = c("vars.remain", "alg"))
  if (verbose) cat("Cleaning up inbix private EC results files\n")
  inbix.temp.files <- c(paste(unique.sim.prefix, ".sim.train.tab", sep = ""),
                        paste(unique.sim.prefix, ".sim.holdout.tab", sep = ""),
                        paste(unique.sim.prefix, ".privateec.privacy.iterations.tab", sep = ""),
                        paste(unique.sim.prefix, ".privateec.log", sep = ""))
  if (sum.good > 2) {
    inbix.temp.files <- c(inbix.temp.files, 
                          paste(unique.sim.prefix, ".sim.test.tab", sep = ""))
  }
  #file.remove(c(inbix.conv.files, rf.files, attrs.file, iters.file))
  file.remove(inbix.temp.files)
  if (!is.null(save.file)) {
    save.file <- paste(unique.sim.prefix, ".Rdata", sep = "")
    if (verbose) {
      cat("saving results to ", save.file, "\n")
    }
    save(fxplots, melted.fx, correct.detect.inbix, bias,
         train.ds, holdout.ds, validation.ds, file = save.file)
  }
  if (verbose) cat("privacyECinbix elapsed time:", (proc.time() - ptm)[3], "\n")
  
  list(algo.acc = fxplots,
       ggplot.data = melted.fx,
       correct = correct.detect.inbix,
       elasped = (proc.time() - ptm)[3])
}

#' Read inbix numeric data set and phenotype file; return combined data frame.
#' 
#' \code{readInbixNumericAsRegressionData} 
#' 
#' @references 
#' \itemize{
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @family inbix interface functions
#' @param baseInbixName \code{string} base filename.
#' @return \code{data.frame} with numeric data in first m columns and phenotye in m+1 column.
#' @export
readInbixNumericAsRegressionData <- function(baseInbixName) {
  inbixNumericFile <- paste(baseInbixName, ".num", sep = "")
  inbixNumericTable <- read.table(inbixNumericFile, header = TRUE, sep = "\t", 
                                  stringsAsFactors = FALSE)
  inbixNumericTable <- inbixNumericTable[, 3:ncol(inbixNumericTable)]
  inbixPhenoFile <- paste(baseInbixName, ".pheno", sep = "")
  inbixPhenoTable <- read.table(inbixPhenoFile, header = FALSE, sep = "\t", 
                                stringsAsFactors =  FALSE)[, 3]
  regressionData <- data.frame(inbixNumericTable, inbixPhenoTable)
  colnames(regressionData) <- c(colnames(inbixNumericTable), "Class")
  regressionData
}

#' Regression genetic association interaction network - reGAIN.
#' 
#' \code{regainInbix} 
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.nature.com/gene/journal/v11/n8/full/gene201037a.html}
#'   {Genes & Immunity Paper}
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @family inbix interface functions
#' @family GAIN functions
#' @seealso \code{\link{regainParallel}}
#' @param regressionData \code{data.frame} with sample rows, gene columns 
#' plus phenotype column.
#' @param stdBetas \code{logical} use standardized beta coefficients (coef/stderr).
#' @param absBetas \code{logical} take absolute value of beta coefficients.
#' @param outPrefix \code{string} file output prefix.
#' @param pThreshold \code{numeric} p-value threshold for GAIN method.
#' @param verbose \code{logical} to send verbose messages to stdout.
#' @return \code{list} of gene scores and associated p-values.
#' @examples
#' data(testdata10)
#' inbixRegain <- regainInbix(testdata10, stdBetas = TRUE, absBetas = TRUE)
#' @export
regainInbix <- function(regressionData, stdBetas = TRUE, absBetas = TRUE, 
                        outPrefix = "Rinbix", pThreshold = 1, verbose = FALSE) {
  inbixExists()
  # write regressionData data frame to inbix files
  writeRegressionDataAsInbixNumeric(regressionData, outPrefix)
  stdBetasCmd <- ""
  if (!stdBetas) {
    stdBetasCmd <- "--regain-use-beta-values"
  }
  absBetasCmd <- ""
  if (absBetas) {
    absBetasCmd <- "--regain-matrix-transform abs"  
  }
  pThresholdCmd <- ""
  if (pThreshold < 1) {
    pThresholdCmd <- paste("--regain-pvalue-threshold", pThreshold)
  }
  verboseCmd <- ""
  if (verbose) {
    verboseCmd <- "--verbose"
  }
  # run inbix reGAIN
  numericFile <- paste(outPrefix, "num ", sep = ".")
  phenoFile <- paste(outPrefix, "pheno", sep = ".")
  inbixCmd <- paste("inbix --regain-minimal", 
                    "--numeric-file", numericFile, 
                    "--pheno", phenoFile, "--1",
                    "--out", outPrefix, 
                    stdBetasCmd, absBetasCmd, pThresholdCmd, verboseCmd)
  if (verbose) cat("Running inbix command:", inbixCmd, "\n")
  system(inbixCmd, intern = TRUE)
  if (verbose) print(list.files(path = ".", pattern = paste(outPrefix, "*")))
  regainMatrixFile <- paste(outPrefix, ".regain.min.tab", sep = "")
  if (file.exists(regainMatrixFile)) {
    inbixRegain <- read.table(regainMatrixFile, header = TRUE, sep = "\t")
  } else {
    print(list.files(path = ".", pattern = paste(outPrefix, "*")))
    stop(paste(regainMatrixFile, "not generated by inbix reGAIN algorithm!"))
  }

  # read warnings and errors if they exist
  warningsText <- ""
  failuresText <- ""
  warningsFile <- paste(outPrefix, ".regression.warnings", sep = "")
  failuresFile <- paste(outPrefix, ".regression.failures", sep = "")
  if (file.exists(warningsFile)) {
    warningsText <- readLines(warningsFile)
    if (verbose) cat(warningsText, "\n")
    file.remove(warningsFile)
  }
  if (file.exists(failuresFile)) {
    failuresText <- readLines(failuresFile)
    if (verbose) cat(failuresText, "\n")
    file.remove(failuresFile)
  }

  # remove temporary files
  tempFileSuffixes <- c("regain.min.tab",
                        "regain.min.pvals.tab",
                        "regain.min.sif",
                        "runinfo.tab",
                        "num",
                        "pheno",
                        "log")
  file.remove(paste(outPrefix, tempFileSuffixes, sep = "."))
  
  # return regain matrix along with any warnings and/or errors
  list(reGAIN = as.matrix(inbixRegain), 
       warningsText = warningsText, 
       failuresText = failuresText)
}

#' Rank by inbix ReliefSeq.
#' 
#' \code{rankReliefSeqInbix} 
#' 
#' @references 
#' \itemize{
#'   \item \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081527}
#'   {ReliefSeq: A variable-Wise Adaptive-K Nearest-Neighbor Feature Selection Tool for Finding Gene-Gene Interactions and Main Effects in mRNA-Seq Gene Expression Data}
#'   \item \url{http://insilico.utulsa.edu/index.php/reliefseq/}{Relief-Seq Software}
#' }
#' @family inbix interface functions
#' @family feature selection functions
#' @param labelledDataFrame \code{data.frame} of predictors and final class column.
#' @param outPrefix \code{string} output file prefix for temporary files.
#' @param k \code{numeric} k-nearest neighbors in the Relief-F algorithm.
#' @return \code{data.frame} with ReliefSeq results: variable, score.
#' @examples
#' data(testdata100ME4)
#' rankReliefSeqResults <- rankReliefSeq(testdata100ME4)
#' @export
rankReliefSeqInbix <- function(labelledDataFrame, outPrefix="Rinbix", k=10) {
  inbixExists()
  # write labelledDataFrame data frame to inbix files
  writeRegressionDataAsInbixNumeric(labelledDataFrame, outPrefix)
  # run reliefseq command
  #relieffCmd <- paste("reliefseq -n Rinbix.num -a Rinbix.pheno -k", k , "-o", outPrefix)
  # replaced relifseq with inbix version 6/24/17
  relieffCmd <- paste("inbix --numeric-file Rinbix.num --pheno Rinbix.pheno --1", 
                      "--relieff --algorithm-mode reliefseq",
                      "--k-nearest-neighbors", k,
                      "--out", outPrefix)
  #cat("Running relieff command:", relieffCmd, "\n")
  system(relieffCmd, intern = TRUE)
  #relieffStdout <- system(relieffCmd, intern=TRUE)
  #cat("stdout:", relieffStdout, "\n")
  relieffRankings <- read.table("Rinbix.relieff.tab", header = FALSE, sep = "\t")
  file.remove(c("Rinbix.num", "Rinbix.pheno", "Rinbix.relieff.tab"))
  # cmd=relieffCmd, stdout=relieffStdout
  data.frame(variable = relieffRankings[, 2], 
             score = relieffRankings[, 1],
             stringsAsFactors = FALSE)
}

#' Write inbix numeric and phenotype files (PLINK format) for a list
#' of simulated data sets (from privateEC, perhaps?)
#'
#' @param data.sets A list of train, holdout and validation data frames
#' @param base.sim.prefix A character vector for the input and saved file prefixes
#' @param verbose A flag for sending verbose output to stdout
#' @return NULL
saveSimAsInbixNative <- function(data.sets=NULL,
                                 base.sim.prefix,
                                 verbose=FALSE) {
  if (is.null(data.sets)) {
    stop("privateEC: No data sets provided as first argument")
  }
  X_train <- data.sets$train
  X_holdo <- data.sets$holdout
  X_test <- data.sets$validation
  
  # train
  train.expr.matrix <- X_train[, 1:(ncol(X_train) - 1)]
  var.names <- colnames(train.expr.matrix)
  train.num.subj <- nrow(train.expr.matrix)
  train.subj.names <- paste("subj", 1:train.num.subj, sep = "")
  train.phenotype <- ifelse(X_train[, ncol(X_train)] == -1, 0, 1)
  train.inbix <- cbind(train.subj.names, train.subj.names, train.expr.matrix)
  colnames(train.inbix) <- c("FID", "IID", var.names)
  utils::write.table(train.inbix, file=paste(base.sim.prefix, ".train.sim.num", sep=""),
                     quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  train.inbix.pheno <- cbind(train.subj.names, train.subj.names, train.phenotype)
  utils::write.table(train.inbix.pheno, file=paste(base.sim.prefix, ".train.sim.pheno", sep=""),
                     quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  # holdout
  holdo.expr.matrix <- X_holdo[, 1:(ncol(X_holdo)-1)]
  var.names <- colnames(holdo.expr.matrix)
  holdo.num.subj <- nrow(holdo.expr.matrix)
  holdo.subj.names <- paste("subj", 1:holdo.num.subj, sep="")
  holdo.phenotype <- ifelse(X_holdo[, ncol(X_holdo)] == -1, 0, 1)
  holdo.inbix <- cbind(holdo.subj.names, holdo.subj.names, holdo.expr.matrix)
  colnames(holdo.inbix) <- c("FID", "IID", var.names)
  utils::write.table(holdo.inbix, file=paste(base.sim.prefix, ".holdo.sim.num", sep=""),
                     quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  holdo.inbix.pheno <- cbind(holdo.subj.names, holdo.subj.names, holdo.phenotype)
  utils::write.table(holdo.inbix.pheno, file=paste(base.sim.prefix, ".holdo.sim.pheno", sep=""),
                     quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
  # validation
  validation.expr.matrix <- X_test[, 1:(ncol(X_test)-1)]
  var.names <- colnames(validation.expr.matrix)
  validation.num.subj <- nrow(validation.expr.matrix)
  validation.subj.names <- paste("subj", 1:validation.num.subj, sep="")
  validation.phenotype <- ifelse(X_test[, ncol(X_test)] == -1, 0, 1)
  validation.inbix <- cbind(validation.subj.names, validation.subj.names, validation.expr.matrix)
  colnames(validation.inbix) <- c("FID", "IID", var.names)
  utils::write.table(validation.inbix, file=paste(base.sim.prefix, ".validation.sim.num", sep=""),
                     quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
  validation.inbix.pheno <- cbind(validation.subj.names, validation.subj.names, validation.phenotype)
  utils::write.table(validation.inbix.pheno, file=paste(base.sim.prefix, ".validation.sim.pheno", sep=""),
                     quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

#' SNPrank network centrality gene ranker based on GAIN matrix.
#' 
#' \code{snprankInbix} 
#' 
#' @references 
#' \itemize{
#'   \item \url{http://www.nature.com/gene/journal/v11/n8/full/gene201037a.html}
#'   {Genes & Immunity Paper}
#'   \item \url{http://insilico.utulsa.edu/index.php/snprank/}{SNPrank Software}
#'   \item \url{https://github.com/hexhead/inbix}{C++ inbix on github}
#' }
#' @family inbix interface functions
#' @family feature selection functions
#' @keywords eigenvectors
#' @seealso \code{\link{snprank}}
#' @param gainMatrix \code{matrix} GAIN matrix.
#' @param outPrefix \code{string} file output prefix.
#' @param gammaParam \code{numeric} gamma damping parameter (see paper).
#' @return \code{data.frame} of gene snpranks.
#' @examples
#' data(testdata10)
#' inbixRegain <- regainInbix(testdata10, stdBetas = TRUE, absBetas = TRUE)
#' inbixSnpranksDF <- snprankInbix(inbixRegain$reGAIN)
#' @export
snprankInbix <- function(gainMatrix, outPrefix = "Rinbix", gammaParam = 0.85) {
  inbixExists()
  # write gainMatrix to inbix-compatible file
  write.table(gainMatrix, file = "Rinbix.gain", quote = FALSE, col.names = TRUE, 
              row.names = FALSE, sep = "\t")
  gammaCmd <- ""
  if (gammaParam > 0) {
    gammaCmd <- paste("--rank-centrality-gamma", gammaParam)
  }
  # run inbix SNPrank
  #cat("gamma:", gamma, "gamma command:", gammaCmd, "\n")
  inbixCmd <- paste("inbix --regain-file Rinbix.gain --rank-by centrality --out", 
                    outPrefix, gammaCmd)
  #cat("Running inbix command:", inbixCmd, "\n")
  #snprankStdout <- system(inbixCmd, intern = T)
  system(inbixCmd, intern = TRUE)
  #print(snprankStdout)
  inbixSnpranks <- read.table("Rinbix.ranks", header = TRUE, sep = "\t")
  
  # remove temporary files
  file.remove(c("Rinbix.gain", "Rinbix.ranks", "Rinbix.log"))
  
  # return snpranks
  data.frame(variable = inbixSnpranks$SNP, score = inbixSnpranks$SNPrank)
}

#' Write regression data frame as  inbix numeric data set and phenotype file.
#' 
#' \code{writeRegressionDataAsInbixNumeric} 
#' 
#' @family inbix interface functions
#' @param regressionData \code{data.frame} with sample rows, gene columns 
#' plus phenotype column.
#' @param baseInbixName \code{string} base filename to write.
#' @export
writeRegressionDataAsInbixNumeric <- function(regressionData, baseInbixName) {
  numSamples <- nrow(regressionData)
  numGenes <- ncol(regressionData) - 1
  phenoCol <- ncol(regressionData)
  phenos <- regressionData[, phenoCol]
  subIds <- c(paste("subj", 1:numSamples, sep = ""))
  phenosTable <- data.frame(subIds, subIds, phenos)
  inbixPhenoFile <- paste(baseInbixName, ".pheno", sep = "")
  write.table(phenosTable, inbixPhenoFile, quote = FALSE, sep = "\t", 
              col.names = FALSE, row.names = FALSE)
  inbixNumericData <- regressionData[, 1:numGenes]
  inbixNumericTable <- data.frame(subIds, subIds, inbixNumericData)
  colnames(inbixNumericTable) <- c("FID", "IID", colnames(inbixNumericData))
  inbixNumericFile <- paste(baseInbixName, ".num", sep = "")
  write.table(inbixNumericTable, inbixNumericFile, quote = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE)
}
