# inbixPredictRnaseq.R - Bill White - 4/10/15
#
# Functions to support Insilico Bioinformatics (inbix) RNASeq 
# classification pipeline:
#   * predictRnaseq
#   * preprocess
#   * filterGenes
#   * classify
#   * getClassificationStats

# ----------------------------------------------------------------------------
#' Predict a response based on RNASeq gene expression.
#' 
#' \code{predictRnaseq} 
#' 
#' @param rnaseqCountsTrain RNASeq counts matrix for training.
#' @param groupLabelsTrain Group labels corresponding to the columns in rnaseqCountsTrain.
#' @param rnaseqCountsTest RNASeq counts matrix for testing.
#' @param groupLabelsTest Group labels corresponding to the columns in rnaseqCountsTest.
#' @param preprocessMethod Pre-processing method: none, scale, log2, log2scale.
#' @param filterMethod Filtering method to reduce the dimensions for classification.
#' @param topN Top number of genes to keep in the filter step.
#' @param classifierMethod Classifier used to build a model and predict new data.
#' @param verbose Verbose flag for more output set TRUE
#' @return List with ranked genes, classification metrics.
#' @export
predictRnaseq <- function(rnaseqCountsTrain, groupLabelsTrain,
                          rnaseqCountsTest, groupLabelsTest,
                          preprocessMethod="none", filterMethod="none", 
                          topN=10, classifierMethod="none", verbose=FALSE) {
  if(verbose) { cat("\tCalling preprocess(", preprocessMethod, ")\n") }
  preprocessResult <- preprocess(preprocessMethod, rnaseqCountsTrain, 
                                 rnaseqCountsTest, verbose)
  if(verbose) {
    cat("Before preprocessing:\n")
    print(rnaseqCountsTrain[1:5, 1:5])
    cat("After preprocessing:\n")
    print(preprocessResult$train[1:5, 1:5])
    
    cat("\tCalling filter(", filterMethod, topN, ")\n")
  }
  filterResult <- filterGenes(filterMethod, preprocessResult$train, groupLabelsTrain,
                              preprocessResult$test, groupLabelsTest, topN, verbose)
  if(verbose) { 
    cat("Genes before filtering:", ncol(preprocessResult$train), 
        "after:", ncol(filterResult$train), "\n")
  }
  
  if(verbose) { cat("\tCalling classify(", classifierMethod, ")\n") }
  classifyStats <- classify(classifierMethod, filterResult$train, groupLabelsTrain,
                            filterResult$test, groupLabelsTest, verbose)
  
  list(stats=classifyStats)
}

# ----------------------------------------------------------------------------
#' Pre-processing step of the predictRnaseq function.
#' 
#' \code{preprocess} 
#' 
#' @param method Pre-processing method: none, scale, log2, log2scale.
#' @param countsTrain RNASeq counts matrix for training.
#' @param countsTest RNASeq counts matrix for testing.
#' @param verbose Verbose flag for more output set TRUE
#' @return List with preprocessed training and testing matrices.
#' @export
preprocess <- function(method="none", countsTrain, countsTest, verbose=FALSE) {
  returnTrain <- countsTrain
  returnTest <- countsTest
  if(method == "scale") {
    if(verbose) { cat("\t\tscale\n") }
    returnTrain <- scale(countsTrain, center=F)
    returnTest <- scale(countsTest, center=F)
  }
  if(method == "log2") {
    if(verbose) { cat("\t\tlog2\n") }
    returnTrain <- log2(countsTrain)
    returnTest <- log2(countsTest)
  }
  if(method == "log2scale") {
    if(verbose) { cat("\t\tlog2scale\n") }
    returnTrain <- log2(countsTrain)
    returnTest <- log2(countsTest)
    returnTrain <- scale(returnTrain)
    returnTest <- scale(returnTest)
  }
  # make return values integer (counts)
  returnTrain <- apply(returnTrain, c(1, 2), function(x) { (as.integer(x)) })
  returnTest <- apply(returnTest, c(1, 2), function(x) { (as.integer(x)) })
  
  list(train=returnTrain, test=returnTest)
}

# ----------------------------------------------------------------------------
#' Filtering step of the predictRnaseq function.
#' 
#' \code{filterGenes} 
#' 
#' @param method Filtering method: none, relieff, edger, deseq2, randomforests.
#' @param dataTrain RNASeq counts matrix for training.
#' @param labelsTrain Group labels for training.
#' @param dataTest RNASeq counts matrix for testing.
#' @param labelsTest Group labels for testing.
#' @param nTopGenes Number of top genes to remain after filtering.
#' @param verbose Verbose flag for more output set TRUE
#' @return List with filtered training and testing data sets.
#' @export
filterGenes <- function(method="none", dataTrain, labelsTrain, dataTest, 
                        labelsTest, nTopGenes, verbose=FALSE) {
  if(method == "edger") {
    if(verbose) { cat("\t\tedgeR\n") }
    y <- edgeR::DGEList(counts=t(dataTrain), group=factor(labelsTrain))
    y <- edgeR::estimateCommonDisp(y)
    dgeResult <- edgeR::exactTest(y)
    topGenes <- rownames(dgeResult$table[order(dgeResult$table$PValue), ])[1:nTopGenes]
  }
  if(method == "deseq2") {
    if(verbose) { cat("\t\tDESeq2\n") }
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=t(dataTrain),
                                          colData=data.frame(response=factor(labelsTrain)),
                                          design=~response)
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    topGenes <- rownames(res[order(res$pvalue), ])[1:nTopGenes]
  }
  if(method == "relieff") {
    if(verbose) { cat("\t\tRelief-F\n") }
    classData <- as.data.frame(cbind(dataTrain, labelsTrain))
    colnames(classData) <- c(colnames(dataTrain), "Class")
    classData$Class <- factor(classData$Class)
    relieff <- CORElearn::attrEval(Class ~ ., classData, estimator="Relief")
    topGenes <- names(sort(relieff, decreasing=T))[1:nTopGenes]
  }
  if(method == "randomforests") {
    if(verbose) { cat("\t\tRandom Forests\n") }
    classData <- as.data.frame(cbind(dataTrain, labelsTrain))
    colnames(classData) <- c(colnames(dataTrain), "Class")
    classData$Class <- factor(classData$Class)
    modelRF <- CORElearn::CoreModel(Class ~ ., classData, model="rf")
    rf <- CORElearn::rfAttrEval(modelRF)
    topGenes <- names(sort(rf, decreasing=T))[1:nTopGenes]
  }
  if(verbose) {
    print("Top genes from filter:")
    print(topGenes)
  }
  list(train=dataTrain[, topGenes], test=dataTest[, topGenes])
}

# ----------------------------------------------------------------------------
#' Classifying step of the predictRnaseq function.
#' 
#' \code{classify} 
#' 
#' @param method Classifier method: none, knn, svm, logreg.
#' @param dataTrain RNASeq counts matrix for training.
#' @param labelsTrain Group labels for training.
#' @param dataTest RNASeq counts matrix for testing.
#' @param labelsTest Group labels for testing.
#' @param verbose Verbose flag for more output set TRUE
#' @return List with classifier stats for method.
#' @export
classify <- function(method="none", dataTrain, labelsTrain, dataTest, 
                     labelsTest, verbose=FALSE) {
  if(method == "svm") {
    if(verbose) { cat("\t\tSVM - Linear\n") }
    yTrain <- factor(ifelse(labelsTrain == 1, -1, 1))
    yTest <- factor(ifelse(labelsTest == 1, -1, 1))
    fitControl <- caret::trainControl(method="cv", number=5)
    fit <- caret::train(dataTrain, yTrain, method="svmLinear", 
                        verbose=TRUE, trControl=fitControl, metric="Accuracy")
    predictionsTrain <- predict(fit, newdata=dataTrain)
    classStatsTrain <- getClassificationStats(yTrain, predictionsTrain, 
                                              classLevels=c(-1, 1))
    predictionsTest <- predict(fit, newdata=dataTest)
    classStatsTest <- getClassificationStats(yTest, predictionsTest, 
                                             classLevels=c(-1, 1))
  }
  if(verbose) { 
    cat("\t\tCV10 ", 
        "Train", 
        "SENS:", classStatsTrain$TPR, 
        "SPEC:", classStatsTrain$SPC, 
        "ACC:", classStatsTrain$ACC, 
        "Test:", 
        "SENS:", classStatsTest$TPR, 
        "SPEC:", classStatsTest$SPC, 
        "ACC:", classStatsTest$ACC, 
        "\n")
  }
  list(statsTrain=classStatsTrain, statsTest=classStatsTest)
}

# -----------------------------------------------------------------------------
#' Compute classification stats from true and predicted binary states.
#'
#' \code{getClassificationStats} 
#' 
#' @param trueClassification Vector of true class labels.
#' @param predictedClassification Vector of predicted class labels.
#' @param classLevels Levels of the binary class.
#' @return Data frame of classification stats.
#' @export
getClassificationStats <- function(trueClassification, 
                                   predictedClassification,
                                   classLevels=c(0,1)) {
  # compute confusion matrix
  confusionMatrix <- table(factor(predictedClassification, levels=classLevels),
                           factor(trueClassification, levels=classLevels))
#   print(trueClassification)
#   print(predictedClassification)
#   print(confusionMatrix)
  
  TP <- confusionMatrix[1, 1]
  FP <- confusionMatrix[1, 2]
  FN <- confusionMatrix[2, 1]
  TN <- confusionMatrix[2, 2]
  
  # calculate classification metrics from contingency table
  TPR <- TP / (TP + FN) # TPR/recall/hit rate/sensitivity
  SPC <- TN / (TN + FP) # TNR/specificity/SPC
  PPV <- TP / (TP + FP) # precision/PPV
  NPV <- TN / (TN + FN) # negative predictive value/NPV
  FPR <- FP / (FP + TN) # fall-out/FPR/false positive rate
  FDR <- FP / (FP + TP) # false discovery rate/FDR
  FNR <- FN / (FN + TP) # false negative rate/FNR/miss rate
  
  # accuracy of the model
  ACC <- (TP + TN) / (TP + FP + TN + FN)
  BAC <- (TPR + SPC) / 2
  F1 <- (2 * TP) / (2 * TP + FP + FN)
  MCC <-
    ((TP * TN) - (FP * FN)) /
    sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  
  # package and return all the computed results
  data.frame(TP=TP, FP=FP, TN=TN, FN=FN, TPR=TPR, SPC=SPC, PPV=PPV, NPV=NPV, 
             FPR=FPR, FDR=FDR, FNR=FNR, ACC=ACC, BAC=BAC, F1=F1, MCC=MCC)
}
