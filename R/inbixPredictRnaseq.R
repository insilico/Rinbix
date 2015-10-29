# ----------------------------------------------------------------------------
# inbixPredictRnaseq.R - Bill White - 4/10/15
#
# Functions to support Insilico Bioinformatics (inbix) RNASeq 
# classification pipeline:
#   * predictRnaseq
#   * preprocessRnaseq
#   * filterRnaseq
#   * classifyRnaseq
  
# ----------------------------------------------------------------------------
#' Predict a response based on RNASeq gene expression.
#' 
#' \code{predictRnaseq} 
#' 
#' @param rnaseqCountsTrain Matrix RNASeq counts.
#' @param groupLabelsTrain Vector group labels.
#' @param rnaseqCountsTest Matrix RNASeq counts.
#' @param groupLabelsTest Vector group labels.
#' @param preprocessMethod String pre-processing method: none, scale, log2, log2scale.
#' @param filterMethod String filtering method to reduce the dimensions for classification.
#' @param topN Numeric top number of genes to keep in the filter step.
#' @param classifierMethod String classifier used to build a model and predict new data.
#' @param verbose Flag verbose for more output set TRUE.
#' @return List with ranked genes, classification metrics.
#' @examples
#' data(simrnaseq)
#' predictResult <- predictRnaseq(rnaseqCountsTrain=predictorsTrain, 
#'                                groupLabelsTrain=responseTrain, 
#'                                rnaseqCountsTest=predictorsTest, 
#'                                groupLabelsTest=responseTest, 
#'                                preprocessMethod="none", 
#'                                filterMethod="randomforests", 
#'                                topN=2, 
#'                                classifierMethod="svm",
#'                                verbose=FALSE)
#' @export
predictRnaseq <- function(rnaseqCountsTrain=NULL,
                          groupLabelsTrain=NULL,
                          rnaseqCountsTest=NULL,
                          groupLabelsTest=NULL,
                          preprocessMethod="none",
                          filterMethod="none", 
                          topN=1,
                          classifierMethod="svm",
                          verbose=FALSE) {
  if(verbose) { cat("\tCalling preprocess(", preprocessMethod, ")\n") }
  preprocessResult <- preprocessRnaseq(preprocessMethod, rnaseqCountsTrain, 
                                 rnaseqCountsTest, verbose)
  if(verbose) {
    cat("Before preprocessing:\n")
    print(rnaseqCountsTrain[1:5, 1:5])
    cat("After preprocessing:\n")
    print(preprocessResult$train[1:5, 1:5])
    cat("\tCalling filter(", filterMethod, topN, ")\n")
  }
  filterResult <- filterRnaseq(filterMethod, preprocessResult$train, groupLabelsTrain,
                               preprocessResult$test, groupLabelsTest, topN, verbose)
  if(verbose) { 
    cat("Genes before filtering:", ncol(preprocessResult$train), 
        "after:", ncol(filterResult$train), "\n")
  }
  
  if(verbose) { cat("\tCalling classifier(", classifierMethod, ")\n") }
  classifyStats <- classifyRnaseq(classifierMethod, filterResult$train, groupLabelsTrain,
                            filterResult$test, groupLabelsTest, verbose)
  
  list(stats=classifyStats)
}

# ----------------------------------------------------------------------------
#' Pre-processing step of the predictRnaseq function.
#' 
#' \code{preprocessRnaseq} 
#' 
#' @param method Pre-processing method: none, scale, log2, log2scale.
#' @param countsTrain RNASeq counts matrix for training.
#' @param countsTest RNASeq counts matrix for testing.
#' @param verbose Verbose flag for more output set TRUE
#' @return List with preprocessed training and testing matrices.
#' @examples
#' data(simrnaseq)
#' preprocessResult <- preprocessRnaseq(method="none", 
#'                                predictorsTrain, 
#'                                predictorsTest, 
#'                                verbose=FALSE)
#' @export
preprocessRnaseq <- function(method="none", countsTrain, countsTest, verbose=FALSE) {
  returnTrain <- countsTrain
  returnTest <- countsTest
  if(method == "scale") {
    if(verbose) { cat("\t\tscale\n") }
    returnTrain <- scale(countsTrain, center=F)
    returnTest <- scale(countsTest, center=F)
  }
  if(method == "log2") {
    if(verbose) { cat("\t\tlog2\n") }
    returnTrain <- log2(countsTrain + 1)
    returnTest <- log2(countsTest + 1)
  }
  if(method == "log2scale") {
    if(verbose) { cat("\t\tlog2scale\n") }
    returnTrain <- log2(countsTrain + 1)
    returnTest <- log2(countsTest + 1)
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
#' \code{filterRnaseq} 
#' 
#' @param method String filtering method: none, relieff, edger, deseq2, randomforests.
#' @param dataTrain Matrix RNASeq counts.
#' @param labelsTrain Vector group labels.
#' @param dataTest Matrix RNASeq counts.
#' @param labelsTest Vector group labels.
#' @param nTopGenes Numeric number of top genes to remain after filtering.
#' @param verbose Flag verbose for more output set TRUE.
#' @return List with filtered training and testing data sets.
#' @examples
#' data(simrnaseq)
#' filteredGenes <- filterRnaseq(method="none", 
#'                              predictorsTrain, 
#'                              responseTrain, 
#'                              predictorsTest, 
#'                              responseTest,
#'                              nTopGenes=10, 
#'                              verbose=FALSE)
#' @export
filterRnaseq <- function(method="none", dataTrain, labelsTrain, dataTest, 
                        labelsTest, nTopGenes, verbose=FALSE) {
  topGenes <- colnames(dataTrain)
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
#' \code{classifyRnaseq} 
#' 
#' @param method String filtering method: none, relieff, edger, deseq2, randomforests.
#' @param dataTrain Matrix RNASeq counts.
#' @param labelsTrain Vector group labels.
#' @param dataTest Matrix RNASeq counts.
#' @param labelsTest Vector group labels.
#' @param verbose Flag verbose for more output set TRUE.
#' @return List with classifier stats for method.
#' @examples
#' data(simrnaseq)
#' classifyStats <- classifyRnaseq(method="none", 
#'                           predictorsTrain, 
#'                           responseTrain, 
#'                           predictorsTest, 
#'                           responseTest,
#'                           verbose=FALSE)
#' @export
classifyRnaseq <- function(method="none", dataTrain, labelsTrain, dataTest, 
                     labelsTest, verbose=FALSE) {
  classStatsTrain <- NULL
  classStatsTest <- NULL
  if(method == "svm") {
    if(verbose) { cat("\t\tSVM - Linear\n") }
    yTrain <- factor(ifelse(labelsTrain == 1, -1, 1))
    yTest <- factor(ifelse(labelsTest == 1, -1, 1))
    fitControl <- caret::trainControl(method="cv", number=5)
    fit <- caret::train(dataTrain, yTrain, method="svmLinear", 
                        verbose=TRUE, trControl=fitControl, metric="Accuracy")
    predictionsTrain <- predict(fit, newdata=dataTrain)
    classStatsTrain <- classifyPredictedValues(yTrain, predictionsTrain, 
                                               classLevels=c(-1, 1))
    predictionsTest <- predict(fit, newdata=dataTest)
    classStatsTest <- classifyPredictedValues(yTest, predictionsTest, 
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
