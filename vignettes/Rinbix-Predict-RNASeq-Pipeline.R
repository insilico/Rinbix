## ------------------------------------------------------------------------
library(Rinbix)
data(simrnaseq)

## ------------------------------------------------------------------------
set.seed(1965)
predictResult <- predictRnaseq(rnaseqCountsTrain=predictorsTrain, 
                               groupLabelsTrain=responseTrain, 
                               rnaseqCountsTest=predictorsTest, 
                               groupLabelsTest=responseTest, 
                               preprocessMethod="none", 
                               filterMethod="randomforests", 
                               topN=10, 
                               classifierMethod="svm",
                               verbose=FALSE)
predictResult

## ------------------------------------------------------------------------
set.seed(1965)
preprocessResult <- preprocessRnaseq(method="log2", 
                                     predictorsTrain, 
                                     predictorsTest, 
                                     verbose=FALSE)
cat("Before preprocessing:\n")
predictorsTrain[1:5, 1:5]
cat("After preprocessing:\n")
preprocessResult$train[1:5, 1:5]

filterResult <- filterRnaseq(method="randomforests", 
                             preprocessResult$train, 
                             responseTrain, 
                             preprocessResult$test, 
                             responseTest,
                             nTopGenes=10, 
                             verbose=FALSE)

cat("Genes before filtering:", ncol(preprocessResult$train), 
    "after:", ncol(filterResult$train), "\n")

classifyStats <- classifyRnaseq(method="svm", 
                                filterResult$train, 
                                responseTrain, 
                                filterResult$test, 
                                responseTest,
                                verbose=FALSE)
classifyStats

