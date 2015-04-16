# inbixPredictRnaseq.R - Bill White - 4/10/15
#
# Functions to support Insilico Bioinformatics (inbix) RNASeq 
# classification pipeline:
#   * predictRnaseq
#   * preprocess
#   * filterGenes
#   * classify

# ----------------------------------------------------------------------------
#' Predict a response based on RNASeq gene expression
#' 
#' \code{predictRnaseq} 
#' 
#' @param rnaseqCountsTrain, RNASeq counts matrix for training
#' @param groupLabelsTrain, group labels corresponding to the columns in rnaseqCountsTrain
#' @param rnaseqCountsTest, RNASeq counts matrix for testing
#' @param groupLabelsTest, group labels corresponding to the columns in rnaseqCountsTest
#' @param preprocessMethod, pre-processing method: none, scale, log2, log2scale
#' @param filterMethod, filtering method to reduce the dimensions for classification
#' @param topN, top number of genes to keep in the filter step
#' @param classifierMethod, classifier used to build a model and predict new data
#' @return list with ranked genes, classification metrics
#' @export
predictRnaseq <- function(rnaseqCountsTrain, groupLabelsTrain,
                          rnaseqCountsTest, groupLabelsTest,
                          preprocessMethod="none", filterMethod="none", 
                          topN=10, classifyMethod="none") {
  cat("\tCalling preprocess(", preprocessMethod, ")\n")
  preprocessResult <- preprocess(preprocessMethod, rnaseqCountsTrain, rnaseqCountsTest)
  cat("Before preprocessing:\n")
  print(rnaseqCountsTrain[1:5, 1:5])
  cat("After preprocessing:\n")
  print(preprocessResult$train[1:5, 1:5])
  
  cat("\tCalling filter(", filterMethod, topN, ")\n")
  filterResult <- filterGenes(filterMethod, preprocessResult$train, groupLabelsTrain,
                              preprocessResult$test, groupLabelsTest, topN)
  cat("Genes before:", ncol(preprocessResult$train), "after:", ncol(filterResult$train), "\n")
  
  cat("\tCalling classify(", classifyMethod, ")\n")
  classifyStats <- classify(classifyMethod, filterResult$train, filterResult$test)
  
  list(stats=classifyStats)
}

# ----------------------------------------------------------------------------
#' Pre-processing step of the predictRnaseq function.
#' 
#' \code{preprocess} 
#' 
#' @param method, pre-processing method: none, scale, log2, log2scale
#' @param countsTrain, RNASeq counts matrix for training
#' @param countsTest, RNASeq counts matrix for testing
#' @return list with preprocessed training and testing matrices
#' @export
preprocess <- function(method="none", countsTrain, countsTest) {
  returnTrain <- countsTrain
  returnTest <- countsTest
  if(method == "scale") {
    cat("\t\tscale\n")
    returnTrain <- scale(countsTrain, center=F)
    returnTest <- scale(countsTest, center=F)
  }
  if(method == "log2") {
    cat("\t\tlog2\n")
    returnTrain <- log2(countsTrain)
    returnTest <- log2(countsTest)
  }
  if(method == "log2scale") {
    cat("\t\tlog2scale\n")
    returnTrain <- log2(countsTrain)
    returnTest <- log2(countsTest)
    returnTrain <- scale(returnTrain)
    returnTest <- scale(returnTest)
  }
  # make return data integer counts
  returnTrain <- apply(returnTrain, c(1, 2), function(x) { (as.integer(x)) })
  returnTest <- apply(returnTest, c(1, 2), function(x) { (as.integer(x)) })
  
  list(train=returnTrain, test=returnTest)
}

# ----------------------------------------------------------------------------
#' Filtering step of the predictRnaseq function.
#' 
#' \code{filter} 
#' 
#' @param method, filtering method: none, relieff, edger, deseq2, randomforests
#' @param dataTrain, RNASeq counts matrix for training
#' @param labelsTrain, group labels for training
#' @param dataTest, RNASeq counts matrix for testing
#' @param labelsTest, group labels for testing
#' @param nTopGenes, number of top genes to remain after filtering
#' @return list with filtered training and testing data sets
#' @export
filterGenes <- function(method="none", dataTrain, labelsTrain, dataTest, labelsTest, nTopGenes) {
  if(method == "edger") {
    require(edgeR)
    cat("\t\tedgeR\n")
    y <- DGEList(counts=t(dataTrain), group=labelsTrain)
    y <- estimateCommonDisp(y)
    dgeResult <- exactTest(y)
    topGenes <- rownames(dgeResult$table[order(dgeResult$table$PValue), ])[1:nTopGenes]
  }
  if(method == "deseq2") {
    require(DESeq2)
    cat("\t\tDESeq2\n")
    dds <- DESeqDataSetFromMatrix(countData=t(predictorsTrain),
                                  colData=data.frame(response=factor(responseTrain)),
                                  design=~response)
    dds <- DESeq(dds)
    res <- results(dds)
    topGenes <- rownames(res[order(res$pvalue), ])[1:nTopGenes]
  }
  if(method == "relieff") {
    require(CORElearn)
    cat("\t\tRelief-F\n")
    classData <- as.data.frame(cbind(dataTrain, labelsTrain))
    colnames(classData) <- c(colnames(dataTrain), "Class")
    relieff <- attrEval(Class ~ ., classData, estimator="Relief")
    topGenes <- names(sort(relieff, decreasing=T))[1:nTopGenes]
  }
  if(method == "randomforests") {
    require(CORElearn)
    cat("\t\tRandom Forests\n")
    classData <- as.data.frame(cbind(dataTrain, labelsTrain))
    colnames(classData) <- c(colnames(dataTrain), "Class")
    modelRF <- CoreModel(Class ~ ., classData, model="rf",
                         selectionEstimator="MDL", minNodeWeightRF=5,
                         rfNoTrees=100, maxThreads=1)
    rf <- rfAttrEval(modelRF)
    topGenes <- names(sort(rf, decreasing=T))[1:nTopGenes]
  }
  print(topGenes)
  list(train=dataTrain[, topGenes], test=dataTest[, topGenes])
}

# ----------------------------------------------------------------------------
#' Classifying step of the predictRnaseq function.
#' 
#' \code{classify} 
#' 
#' @param method, classifier method: none, knn, svm, logreg
#' @param dataTrain, RNASeq counts matrix for training
#' @param dataTest, RNASeq counts matrix for testing
#' @return list with classifier stats for method
#' @export
classify <- function(method="none", dataTrain, labelsTrain, dataTest, labelsTest) {
  list()
}

