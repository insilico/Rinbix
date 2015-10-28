library(Rinbix)
context("Data Simulations")

test_that("createRandomRegressionDataset", {
  ds <- createRandomRegressionDataset(100, 100)
	expect_equal(nrow(ds), 100)
  expect_equal(ncol(ds), 101)
  expect_equal(nlevels(factor(ds[, ncol(ds)])), 2)
})

test_that("createRandomMatrix", {
  ds <- createRandomMatrix(100, 100, 7, 0.05)$regressionData
  expect_equal(nrow(ds), 100)
  expect_equal(ncol(ds), 101)
  expect_equal(nlevels(factor(ds[, ncol(ds)])), 2)
})

test_that("createDiffCoexpMatrix", {
  data("scaleFreeNetwork")
  dsobj <- createDiffCoexpMatrix(M=100, 
                               N=100, 
                               meanExpression=7, 
                               A=scaleFreeNetwork, 
                               randSdNoise=0.05, 
                               sdNoise=1.5, 
                               mGenesToPerturb=3,
                               sampleIndicesMainEffects = c(5, 10, 15),
                               sampleIndicesInteraction = c(5, 10, 15),
                               mainEffectMode=1,
                               mainEffect=4,
                               verbose=F)
  ds <- dsobj$regressionData  
  expect_equal(nrow(ds), 100)
  expect_equal(ncol(ds), 101)
  expect_equal(nlevels(factor(ds[, ncol(ds)])), 2)
})

test_that("createDiffCoexpMatrixNoME", {
  data("scaleFreeNetwork")
  dsobj <- createDiffCoexpMatrixNoME(M=100, 
                               N=100, 
                               meanExpression=7, 
                               A=scaleFreeNetwork, 
                               randSdNoise=0.05, 
                               sdNoise=1.5, 
                               sampleIndicesInteraction = c(5, 10, 15))
  ds <- dsobj$regressionData  
  expect_equal(nrow(ds), 100)
  expect_equal(ncol(ds), 101)
  expect_equal(nlevels(factor(ds[, ncol(ds)])), 2)
})

test_that("createDiffCoexpMatrixNull", {
  data("scaleFreeNetwork")
  dsobj <- createDiffCoexpMatrixNull(M=100, 
                               N=100, 
                               meanExpression=7, 
                               A=scaleFreeNetwork, 
                               randSdNoise=0.05, 
                               sdNoise=1.5)
  ds <- dsobj$regressionData  
  expect_equal(nrow(ds), 100)
  expect_equal(ncol(ds), 101)
  expect_equal(nlevels(factor(ds[, ncol(ds)])), 2)
})

test_that("createMainEffectsMatrix", {
  dsobj <- createMainEffectsMatrix(M=100, 
                               N=100, 
                               meanExpression=7, 
                               randSdNoise=0.05, 
                               sampleIndicesMainEffects = c(5, 10),
                               mainEffect=4, 
                               doScale=FALSE, 
                               doLog=FALSE)
  ds <- dsobj$regressionData  
  expect_equal(nrow(ds), 100)
  expect_equal(ncol(ds), 101)
  expect_equal(nlevels(factor(ds[, ncol(ds)])), 2)
})

test_that("getFoldChange", {
  dsobj <- createMainEffectsMatrix(M=100, 
                               N=100, 
                               meanExpression=7, 
                               randSdNoise=0.05, 
                               sampleIndicesMainEffects = c(5, 10),
                               mainEffect=4, 
                               doScale=FALSE, 
                               doLog=FALSE)
  ds <- dsobj$regressionData  
  fc <- getFoldChange(t(ds), 5, 50, 100)
  expect_equal(object=fc, expected=4, tolerance=0.1)
})

test_that("simCorrMatrix", {
  mat <- simCorrMatrix(n=400, num_clust=20, max_noise_corr=0.5, lower_true_corr=0.5)
  expect_equal(nrow(mat), 400)
  expect_equal(ncol(mat), 400)
})

test_that("simCorrMatrixNonUniform", {
  mat <- simCorrMatrixNonUniform(n=400, num_clust=20, max_noise_corr=0.5, lower_true_corr=0.5)
  expect_equal(nrow(mat), 400)
  expect_equal(ncol(mat), 400)
})
