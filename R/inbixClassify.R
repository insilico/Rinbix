# ------------------------------------------------------------------------------------
# inbixClassify.R - Bill White - 10/10/15
#
# Rinbix package machine learning classification functions.

# ------------------------------------------------------------------------------------
#' Compute and return classifier stats from a confusion matrix.
#' 
#' \code{classifyConfusionMatrix} 
#' 
#' @family classification functions
#' @param confusionMatrix \code{matrix} of classification counts.
#' @return \code{data.frame} of classifier stats.
#' @examples
#' testValues <- c(0,0,0,0,0,1,1,1,1,1)
#' trueValues <- c(0,1,0,1,0,1,0,1,0,0)
#' classifierStats <- classifyConfusionMatrix(table(testValues, trueValues))
#' @export
classifyConfusionMatrix <- function(confusionMatrix) {
  TN <- confusionMatrix[1, 1]
  FN <- confusionMatrix[1, 2]
  FP <- confusionMatrix[2, 1]
  TP <- confusionMatrix[2, 2]
  
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

# -----------------------------------------------------------------------------
#' Classify one variable using Weka Logistic
#' 
#' \code{classifyOneVarWeka} 
#' 
#' @family classification functions
#' @param varName \code{string} variable name.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @return \code{list} of classifier stats.
#' @examples
#' data(testdata100ME4)
#' testdata100ME4$Class <- factor(testdata100ME4$Class)
#' classifierStats <- classifyOneVarWeka("gene0005", testdata100ME4, testdata100ME4)
#' @export
classifyOneVarWeka <- function(varName, trainData, testData) {
  fit_formula <- paste("Class ~", varName)
  weka_model <- RWeka::Logistic(as.formula(fit_formula), data=trainData)
  weka_model_acc <- summary(weka_model)$details["pctCorrect"]
  weka_model_eval <- RWeka::evaluate_Weka_classifier(weka_model, newdata=testData, 
                                                     complexity=FALSE, class=TRUE)
  weka_model_eval_acc <- weka_model_eval$details["pctCorrect"]
  
  list(train_acc=weka_model_acc, test_acc=weka_model_eval_acc)
}

# -----------------------------------------------------------------------------
#' Classify a pair variable using Weka Logistic.
#' 
#' \code{classifyPairWeka} 
#' 
#' @family classification functions
#' @param var1Name \code{string} variable name.
#' @param var2Name \code{string} variable name.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @return \code{list} of classifier stats.
#' @examples
#' data(testdata100ME4)
#' testdata100ME4$Class <- factor(testdata100ME4$Class)
#' classifierStats <- classifyPairWeka("gene0005", "gene0010", 
#'                                     testdata100ME4, testdata100ME4)
#' @export
classifyPairWeka <- function(var1Name, var2Name, trainData, testData) {
  fit_formula <- paste("Class ~", var1Name, "+" , var2Name)
  weka_model <- RWeka::Logistic(as.formula(fit_formula), data=trainData)
  weka_model_acc <- summary(weka_model)$details["pctCorrect"]
  weka_model_eval <- RWeka::evaluate_Weka_classifier(weka_model, newdata=testData, 
                                              complexity=FALSE, class=TRUE)
  weka_model_eval_acc <- weka_model_eval$details["pctCorrect"]
  weka_model_eval_hr <- RWeka::evaluate_Weka_classifier(weka_model, newdata=testData, 
                                                 complexity=FALSE, class=TRUE)
  list(train_acc=weka_model_acc, test_acc=weka_model_eval_acc)
}

# ------------------------------------------------------------------------------------
#' Compute and return classifier stats from classifier predicted values versus 
#' true clasifications.
#' 
#' \code{classifyPredictedValues} 
#' 
#' @family classification functions
#' @param someClassification \code{vector} of predicted values.
#' @param trueClassification \code{vector} of true values.
#' @param classLevels \code{vector} of class levels.
#' @return \code{data.frame} of classifier stats.
#' @examples
#' testValues <- c(0,0,0,0,0,1,1,1,1,1)
#' trueValues <- c(0,1,0,1,0,1,0,1,0,0)
#' classifierStats <- classifyPredictedValues(testValues, trueValues)
#' @export
classifyPredictedValues <- function(someClassification, 
                                    trueClassification,
                                    classLevels=c(0,1)) {
  confusionMatrix <- table(factor(someClassification, levels=classLevels),
                           factor(trueClassification, levels=classLevels))
  classifyConfusionMatrix(confusionMatrix)
}

# ------------------------------------------------------------------------------------
#' Classify and predict a train and test data set pair for a cross validation fold.
#' 
#' \code{classifyPredictFold} 
#' 
#' @param foldIdx \code{numeric} fold index.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @param my_seed \code{numeric} random generator seed.
#' @param samp_method \code{string} over, under or none sampling.
#' @param wrapper \code{string} feature slection algorithm.
#' @param top_n \code{numeric} top n features to select.
#' @return \code{data.frame} of classification and prediction statistics.
#' @export
classifyPredictFold <- function(foldIdx, trainData, testData, 
                                my_seed=5627, samp_method="orig",
                                wrapper="none", top_n=ncol(trainData)-1) {
  grp_table_trn <- table(trainData$Class)
  grp_table_tst <- table(testData$Class)
  gene_expr <- trainData[, -ncol(trainData)]
  num_genes <- ncol(gene_expr)
  gene_names <- colnames(gene_expr)
  #print(grp_table_trn)
  # rank the variables before training
  ranked_vars <- seq(from=num_genes, to=1)
  names(ranked_vars) <- gene_names
  top_vars <- colnames(gene_expr)
  if(wrapper == "relieff") {
    ranked_vars <- CORElearn::attrEval(Class ~ ., trainData, estimator="Relief")
  }
  if(wrapper == "rf") {
    rf_model <- CORElearn::CoreModel(Class ~ ., trainData, model="rf")
    ranked_vars <- CORElearn::rfAttrEval(rf_model)
  }
  if(wrapper == "glmnet") {
    ranked_vars <- rankGlmnet(trainData)$names
  }
  if(wrapper == "ttest") {
    t_test_pvals <- vector(mode="numeric", length=num_genes)
    names(t_test_pvals) <- gene_names
    for(gene_idx in 1:num_genes) {
      t_test_pvals[gene_idx] <-  t.test(trainData[, gene_idx] ~ trainData$Class)$p.value
    }
    ranked_vars <- t_test_pvals
  }
  if(wrapper == "snprank") {
    if(num_genes < 3) {
      cat("WARNING: SNPrank selection wrapper requires at least 3 variables\n")
      cat("WARNING: Continuing with unranked variables\n")
    } else {
      # run reGAIN
      regain_results <- regainParallel(trainData)
      # run SNPrank on reGAIN matrix
      regain_matrix <- regain_results
      colnames(regain_matrix) <- colnames(trainData[, -ncol(trainData)])
      snprank_results <- snprank(regain_matrix)
      ranked_vars <- snprank_results$snprank
      names(ranked_vars) <- as.character(snprank_results$gene)
    }
  }
  num_ranked_vars <- length(ranked_vars)
  if(num_ranked_vars < top_n) {
    #     cat("WARNING glmnet selected less than specified top N:", top_n)
    #     cat(" setting top N to length glnmnet selection:", num_ranked_vars, "\n")
    top_n <- num_ranked_vars
  }
  if(num_ranked_vars < 1) {
    fold_train_acc <- 0
    fold_test_acc <- 0
  } else {
    if(wrapper == "ttest") {
      top_vars <- names(sort(ranked_vars, decreasing=F)[1:top_n])
    } else {
      top_vars <- names(sort(ranked_vars, decreasing=T)[1:top_n])
    }
#     cat("Fold:", foldIdx, "\t", wrapper, " classifier selected top", top_n, 
#         "variables: ", top_vars, "\n")
    # Weka J48, C4.5-like decision tree
    # training
    trainData <- trainData[, c(top_vars, "Class")]
    trainData$Class <- as.factor(trainData$Class)
    # For J48 options use: control= Weka_control(M=1,U=TRUE)
    # use WOW(J48) to see all options
    fit <- RWeka::J48(Class~., data=trainData, control=RWeka::Weka_control(M=1, U=FALSE))
    #fit <- J48(Class~., data=trainData)
    eval_train <- RWeka::evaluate_Weka_classifier(fit, newdata=trainData, numFolds=0, class=TRUE)
    eval_stats <- classifyConfusionMatrix(eval_train$confusionMatrix)
    fold_train_acc <- eval_stats$F1
    # testing
    testData$Class <- as.factor(testData$Class)
    eval_test <- RWeka::evaluate_Weka_classifier(fit, newdata=testData, numFolds=0, class=TRUE)
    eval_stats <- classifyConfusionMatrix(eval_test$confusionMatrix)
    fold_test_acc <- eval_stats$F1

    # caret models
    # ctrl <- trainControl(method="cv", number=5)
    # ctrl <- trainControl(method="none")
    # fit <- train(Class~., data=trainData, method="LMT", trConrol=ctrl)
    # pred <- predict(fit, trainData[, top_vars])
    # pred_metrics <- confusionMatrix(pred, trainData$Class)
    # fold_train_acc <- pred_metrics$byClass["Balanced Accuracy"]
    # # testing - predict using final model - balanced accuracy
    # pred <- predict(fit, testData[, top_vars])
    # pred_metrics <- confusionMatrix(pred, testData$Class)
    # fold_test_acc <- pred_metrics$byClass["Balanced Accuracy"]
  }
  # return summary info
  total_subjects <- grp_table_trn[1] + grp_table_trn[2]
  cat(foldIdx, 
      samp_method, 
      grp_table_trn[1], round(grp_table_trn[1] / total_subjects, 6),
      grp_table_trn[2], round(grp_table_trn[2] / total_subjects, 6),
      grp_table_tst[1], round(grp_table_tst[1] / total_subjects, 6),
      grp_table_tst[2], round(grp_table_tst[2] / total_subjects, 6),
      wrapper, 
      top_n,
      round(fold_train_acc, 6), 
      round(fold_test_acc, 6),
      paste(top_vars, collapse=", "),
      "\n", sep="\t")
  data.frame(fold=foldIdx,
             method=samp_method,
             trn0=grp_table_trn[1],
             trn1=grp_table_trn[2],
             tst0=grp_table_tst[1],
             tst1=grp_table_tst[2],
             trn=fold_train_acc,
             tst=fold_test_acc,
             wrapper=wrapper,
             topn=top_n,
             topvars=paste(top_vars, collapse=", "))
}

# ------------------------------------------------------------------------------------
#' Cross validated classify and predict outer loop.
#' 
#' \code{crossValidate} 
#' 
#' @family classification functions
#' @param classData \code{data.frame} of subject by variable plus 'Class' column.
#' @param k_folds \code{numeric} number of cross validation folds.
#' @param repeat_cv \code{numeric} number of times to repeat k-folds cross validation.
#' @param my_seed \code{numeric} random number generator seed.
#' @param samp_method \code{string} over, under or none sampling.
#' @param wrapper \code{string} feature slection algorithm.
#' @param top_n \code{numeric} top n features to select.
#' @return \code{list} of all results and averaged results.
#' @examples
#' data(testdata100ME4)
#' cv_res <- crossValidate(testdata100ME4, k_folds=10, repeat_cv=1, my_seed=5627, 
#'                         samp_method="orig", wrapper="none", top_n=ncol(testdata100ME4)-1)
#' @export
crossValidate <- function(classData, k_folds=10, repeat_cv=1, my_seed=5627, 
  samp_method="orig", wrapper="none", top_n=ncol(classData)-1) {
  results <- NULL
  # repeated CV
  for(repeat_idx in 1:repeat_cv) {
    cat("Repeat:", repeat_idx, "\n")
    # k-fold cross validation with feature selection CV wrapper and various down/up sampling 
    cat("Creating ", k_folds, " folds for cross validation\n", sep="")
    split_idx <- caret::createFolds(classData$Class, k=k_folds, list=TRUE, returnTrain=FALSE)
    # CV split loops
    cat("---------------------------------------------------------------------------\n")
    cat("Fold\tRSmpl\ttrn0\ttrn1\ttrn0%\ttrn1%\ttst0\ttst1\ttst0%\ttst1%\tWrapper\tTop N\tTrain Acc\tTest Acc\tTop Vars\n")
    for(fold_idx in 1:length(split_idx)) {
      fold_instance_idx <- split_idx[[fold_idx]]
      imbal_train <- classData[-fold_instance_idx, ]
      imbal_test  <- classData[fold_instance_idx, ]
      this_result <- classifyPredictFold(fold_idx, imbal_train, imbal_test, 
                                           my_seed, samp_method, wrapper, top_n=top_n)
      results <- rbind(results, this_result)
    }
  }
  rownames(results) <- paste("run", 1:nrow(results), sep="")
  
  # ------------------------------------------------------------------------------------
  # average errors
  metric_avg <- colMeans(results[results$method == samp_method & 
                                   results$wrapper == wrapper, c("trn", "tst")])
  
  list(results=results, avgs=metric_avg)
}

# -----------------------------------------------------------------------------
#' Compute main effect GLM main effect for a variable name.
#' 
#' \code{glmMainEffect} 
#' 
#' @keywords models regression array
#' @family classification functions
#' @param varName \code{string} variable name.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @return \code{data.frame} of model fit information.
#' @examples
#' data(testdata100ME4)
#' testdata100ME4$Class <- factor(testdata100ME4$Class)
#' classifierStats <- glmMainEffect("gene0005", testdata100ME4, testdata100ME4)
#' @export
glmMainEffect <- function(varName, trainData, testData) {
  regression_formula <- paste("Class ~", varName, sep="")
  fit <- glm(as.formula(regression_formula), data=trainData, family="binomial")
  fit_ys <- predict(fit, newdata=testData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- testData[, ncol(testData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  #print(classification_stats)
  
  maineffect_term_idx <- 2
  maineffect_coeff <- summary(fit)$coefficients[maineffect_term_idx, "Estimate"]
  maineffect_stdcoeff <- summary(fit)$coefficients[maineffect_term_idx, "z value"]
  maineffect_stderr <- summary(fit)$coefficients[maineffect_term_idx, "Std. Error"]
  maineffect_pval <- summary(fit)$coefficients[maineffect_term_idx, "Pr(>|z|)"]
  data.frame(vara=varName, 
             converged=fit$converged, 
             coef=maineffect_coeff, 
             stdcoef=maineffect_stdcoeff,
             stderr=maineffect_stderr, 
             pval=maineffect_pval,
             accuracy=classification_stats$ACC)
}

# -----------------------------------------------------------------------------
#' Compute glm on variable indices.
#' 
#' \code{glmVarList} 
#' 
#' @keywords models regression array
#' @family classification functions
#' @param varIdx \code{string} variable name.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @return \code{data.frame} of classification results.
#' @examples
#' data(testdata100ME4)
#' testdata100ME4$Class <- factor(testdata100ME4$Class)
#' classifierStats <- glmVarList(c(5, 10), testdata100ME4, testdata100ME4)
#' @export
glmVarList <- function(varIdx, trainData, testData) {
  var_names <- colnames(trainData[, 1:(ncol(trainData)-1)])

  regression_formula <- paste("Class ~", paste(var_names[varIdx], sep="+"), sep="")
  fit <- glm(as.formula(regression_formula), data=trainData, family="binomial")

  fit_ys <- predict(fit, newdata=trainData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- trainData[, ncol(trainData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  train_acc <- classification_stats$ACC

  fit_ys <- predict(fit, newdata=testData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- testData[, ncol(testData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  test_acc <- classification_stats$ACC
  
  data.frame(converged=fit$converged, train.acc=train_acc, test.acc=test_acc)
}

# -----------------------------------------------------------------------------
#' Compute glm with a pair of variables including their interaction.
#' 
#' \code{glmWithInteractionTerm} 
#' 
#' @keywords models regression array
#' @family classification functions
#' @param var1Name \code{string} variable name.
#' @param var2Name \code{string} variable name.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @return \code{data.frame} of classification results.
#' @examples
#' data(testdata100ME4)
#' testdata100ME4$Class <- factor(testdata100ME4$Class)
#' classifierStats <- glmWithInteractionTerm("gene0005", "gene0010", testdata100ME4, testdata100ME4)
#' @export
glmWithInteractionTerm <- function(var1Name, var2Name, trainData, testData) {
  regression_formula <- paste("Class ~", var1Name, " + ", var2Name, " + ",
                              var1Name, "*", var2Name, sep="")
  fit <- glm(as.formula(regression_formula), data=trainData, family="binomial")
  fit_ys <- predict(fit, newdata=trainData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- trainData[, ncol(trainData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  train_acc <- classification_stats$ACC

  fit_ys <- predict(fit, newdata=testData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- testData[, ncol(testData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  test_acc <- classification_stats$ACC
  
  # get the interaction term model results and keep appending rows 
  # to a results data frame
  int_term_idx <- 4
  interaction_coeff <- summary(fit)$coefficients[int_term_idx, "Estimate"]
  interaction_stdcoeff <- summary(fit)$coefficients[int_term_idx, "z value"]
  interaction_stderr <- summary(fit)$coefficients[int_term_idx, "Std. Error"]
  interaction_pval <- summary(fit)$coefficients[int_term_idx, "Pr(>|z|)"]
  data.frame(vara=var1Name, 
             varb=var2Name,
             converged=fit$converged, 
             coef=interaction_coeff, 
             stdcoef=interaction_stdcoeff,
             stderr=interaction_stderr, 
             pval=interaction_pval,
             accuracy=classification_stats$ACC,
             train.acc=train_acc,
             test.acc=test_acc)
}

# -----------------------------------------------------------------------------
#' Compute glm with a pair of variables including their interaction and squared interaction.
#' 
#' \code{glmWithSquaredTerms} 
#' 
#' @keywords models regression array
#' @family classification functions
#' @param var1Name \code{string} variable name.
#' @param var2Name \code{string} variable name.
#' @param trainData \code{data.frame} with class column.
#' @param testData \code{data.frame} with class column.
#' @return \code{data.frame} of classification results.
#' @examples
#' data(testdata100ME4)
#' testdata100ME4$Class <- factor(testdata100ME4$Class)
#' classifierStats <- glmWithSquaredTerms("gene0005", "gene0010", testdata100ME4, testdata100ME4)
#' @export
glmWithSquaredTerms <- function(var1Name, var2Name, trainData, testData) {
  regression_formula <- paste("Class ~", var1Name, " + ", var2Name, " + ",
                              var1Name, "*", var2Name, 
                              " + I(",  var1Name, "^2)",
                              " + I(",  var2Name, "^2)", sep="")
  fit <- glm(as.formula(regression_formula), data=trainData, family="binomial")
  fit_ys <- predict(fit, newdata=trainData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- trainData[, ncol(trainData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  train_acc <- classification_stats$ACC

  fit_ys <- predict(fit, newdata=testData, type="response")
  fit_phenos <- ifelse(fit_ys > 0.5, 1, 0)    
  true_phenos <- testData[, ncol(testData)]
  classification_stats <- classifyPredictedValues(true_phenos, fit_phenos)
  test_acc <- classification_stats$ACC
  
  # get the interaction term model results and keep appending rows 
  # to a results data frame
  int_term_idx <- 4
  interaction_coeff <- summary(fit)$coefficients[int_term_idx, "Estimate"]
  interaction_stdcoeff <- summary(fit)$coefficients[int_term_idx, "z value"]
  interaction_stderr <- summary(fit)$coefficients[int_term_idx, "Std. Error"]
  interaction_pval <- summary(fit)$coefficients[int_term_idx, "Pr(>|z|)"]
  data.frame(vara=var1Name, 
             varb=var2Name,
             converged=fit$converged, 
             coef=interaction_coeff, 
             stdcoef=interaction_stdcoeff,
             stderr=interaction_stderr, 
             pval=interaction_pval,
             accuracy=classification_stats$ACC,
             train.acc=train_acc,
             test.acc=test_acc)
}
