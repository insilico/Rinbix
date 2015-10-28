library(Rinbix)
context("Classify")

test_that("classifyConfusionMatrix", {
	testValues <- matrix(c(0,0,0,0,0,1,1,1,1,1))
	trueValues <- matrix(c(0,1,0,1,0,1,0,1,0,0))
	classifierStats <- classifyConfusionMatrix(table(testValues, trueValues))
  expect_equal(classifierStats$ACC, 0.5)
})

test_that("classifyOneVarWeka", {
	require(RWeka)
	data(testdata100ME4)
	testdata100ME4$Class <- factor(testdata100ME4$Class)
	classifierStats <- classifyOneVarWeka("gene0005", testdata100ME4, testdata100ME4)
	#expect_equal(classifierStats$train_acc, 100)
  expect_equal(classifierStats$train_acc, classifierStats$test_acc)
})

test_that("classifyPairWeka", {
	require(RWeka)
	data(testdata100ME4)
	testdata100ME4$Class <- factor(testdata100ME4$Class)
	classifierStats <- classifyPairWeka("gene0005", "gene0010", testdata100ME4, testdata100ME4)
	#expect_equal(classifierStats$train_acc, 100)
  expect_equal(classifierStats$train_acc, classifierStats$test_acc)
})

test_that("classifyPredictedValues", {
	testValues <- matrix(c(0,0,0,0,0,1,1,1,1,1))
	trueValues <- matrix(c(0,1,0,1,0,1,0,1,0,0))
	classifierStats <- classifyPredictedValues(testValues, trueValues)
  expect_equal(classifierStats$ACC, 0.5)
})

# test_that("crossValidate", {
# 	data(testdata100ME4)
# 	cv_res <- crossValidate(testdata100ME4, k_folds=10, repeat_cv=1, my_seed=5627, 
#   	        					    samp_method="orig", wrapper="none", top_n=ncol(testdata100ME4)-1)
#   expect_equal(TRUE, TRUE)
# })

test_that("glmMainEffect", {
	data(testdata100ME4)
	testdata100ME4$Class <- factor(testdata100ME4$Class)
	classifierStats <- glmMainEffect("gene0005", testdata100ME4, testdata100ME4)
  expect_equal(classifierStats$accuracy, 1)
})

test_that("glmVarList", {
	data(testdata100ME4)
	testdata100ME4$Class <- factor(testdata100ME4$Class)
	classifierStats <- glmVarList(c(5, 10), testdata100ME4, testdata100ME4)
  expect_equal(classifierStats$train.acc, 1)
})

test_that("glmWithInteractionTerm", {
	data(testdata100ME4)
	testdata100ME4$Class <- factor(testdata100ME4$Class)
	classifierStats <- glmWithInteractionTerm("gene0005", "gene0010", testdata100ME4, testdata100ME4)
  expect_equal(classifierStats$accuracy, 1)
})

test_that("glmWithSquaredTerms", {
	data(testdata100ME4)
	testdata100ME4$Class <- factor(testdata100ME4$Class)
	classifierStats <- glmWithSquaredTerms("gene0005", "gene0010", testdata100ME4, testdata100ME4)
  expect_equal(classifierStats$accuracy, 1)
})
