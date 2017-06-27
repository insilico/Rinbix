library(Rinbix)
data(testdata10)
context("Clustering and Network Modularity")

test_that("Rinbix modularity from reGAIN stdBetas=TRUE, absBetas=TRUE", {
	require(fossil)
  inbixModulesDF <- read.table("testdata10-abs-zval.modules", header=F)
  inbixModulesGrpsKey <- as.integer(substr(inbixModulesDF$V1, 4, length(inbixModulesDF$V1)))
  inbixModulesGrps <- inbixModulesDF[order(inbixModulesGrpsKey), ]
  colnames(inbixModulesGrps) <- c("Gene", "Group")
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
  rinbixModulesDF <- Rinbix::modularity(rinbixRegain)
  moduleListGrps <- as.data.frame(rinbixModulesDF$groups)
  moduleListGrps$Gene <- as.character(moduleListGrps$Gene)
  moduleListGrps$Group <- as.integer(moduleListGrps$Group)
  moduleListGrpsKey <- as.integer(substr(moduleListGrps$Gene, 4, length(moduleListGrps$Gene)))
  moduleListGrps <- moduleListGrps[order(moduleListGrpsKey), ]
  expect_equal(fossil::rand.index(inbixModulesGrps$Group, inbixModulesGrps$Group), 1)
})

test_that("Rinbix modularity from reGAIN stdBetas=FALSE, absBetas=TRUE", {
	require(fossil)
  inbixModulesDF <- read.table("testdata10-abs-beta.modules", header=F)
  inbixModulesGrpsKey <- as.integer(substr(inbixModulesDF$V1, 4, length(inbixModulesDF$V1)))
  inbixModulesGrps <- inbixModulesDF[order(inbixModulesGrpsKey), ]
  colnames(inbixModulesGrps) <- c("Gene", "Group")
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=TRUE)
  rinbixModulesDF <- Rinbix::modularity(rinbixRegain)
  moduleListGrps <- as.data.frame(rinbixModulesDF$groups)
  moduleListGrps$Gene <- as.character(moduleListGrps$Gene)
  moduleListGrps$Group <- as.integer(moduleListGrps$Group)
  moduleListGrpsKey <- as.integer(substr(moduleListGrps$Gene, 4, length(moduleListGrps$Gene)))
  moduleListGrps <- moduleListGrps[order(moduleListGrpsKey), ]
  expect_equal(fossil::rand.index(inbixModulesGrps$Group, inbixModulesGrps$Group), 1)
})

test_that("Rinbix modularity from reGAIN stdBetas=TRUE, absBetas=FALSE", {
	require(fossil)
  inbixModulesDF <- read.table("testdata10-noabs-zval.modules", header=F)
  inbixModulesGrpsKey <- as.integer(substr(inbixModulesDF$V1, 4, length(inbixModulesDF$V1)))
  inbixModulesGrps <- inbixModulesDF[order(inbixModulesGrpsKey), ]
  colnames(inbixModulesGrps) <- c("Gene", "Group")
  rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=FALSE)
  rinbixModulesDF <- Rinbix::modularity(rinbixRegain)
  moduleListGrps <- as.data.frame(rinbixModulesDF$groups)
  moduleListGrps$Gene <- as.character(moduleListGrps$Gene)
  moduleListGrps$Group <- as.integer(moduleListGrps$Group)
  moduleListGrpsKey <- as.integer(substr(moduleListGrps$Gene, 4, length(moduleListGrps$Gene)))
  moduleListGrps <- moduleListGrps[order(moduleListGrpsKey), ]
  expect_equal(fossil::rand.index(inbixModulesGrps$Group, inbixModulesGrps$Group), 1)
})

test_that("Rinbix modularity from reGAIN stdBetas=FALSE, absBetas=FALSE", {
	require(fossil)
  inbixModulesDF <- read.table("testdata10-noabs-beta.modules", header=F)
  inbixModulesGrpsKey <- as.integer(substr(inbixModulesDF$V1, 4, length(inbixModulesDF$V1)))
  inbixModulesGrps <- inbixModulesDF[order(inbixModulesGrpsKey), ]
  colnames(inbixModulesGrps) <- c("Gene", "Group")
  rinbixRegain <- regainParallel(testdata10, stdBetas=FALSE, absBetas=FALSE)
  rinbixModulesDF <- Rinbix::modularity(rinbixRegain)
  moduleListGrps <- as.data.frame(rinbixModulesDF$groups)
  moduleListGrps$Gene <- as.character(moduleListGrps$Gene)
  moduleListGrps$Group <- as.integer(moduleListGrps$Group)
  moduleListGrpsKey <- as.integer(substr(moduleListGrps$Gene, 4, length(moduleListGrps$Gene)))
  moduleListGrps <- moduleListGrps[order(moduleListGrpsKey), ]
  expect_equal(fossil::rand.index(inbixModulesGrps$Group, inbixModulesGrps$Group), 1)
})

test_that("rip-M", {
  simMatrix <- simCorrMatrix(n=400, 
                             num_clust=20, 
                             max_noise_corr=0.8, 
                             lower_true_corr=0.2) 
  modListRipm <- ripM(simMatrix, 
                      thresholdType="hard", 
                      thresholdValue=0.8, 
                      startMergeOrder=2, 
                      maxMergeOrder=4, 
                      minModuleSize=10, 
                      maxModuleSize=50, 
                      useAbs=TRUE, 
                      useWeighted=TRUE,
                      verbose=FALSE)
  expect_equal(TRUE, TRUE)
})

test_that("rip-M File I/O", {
  simMatrix <- simCorrMatrix(n=400, num_clust=20, max_noise_corr=0.8, lower_true_corr=0.2) 
  modListRipm <- ripM(simMatrix, 
                      thresholdType="hard", 
                      thresholdValue=0.8, 
                      startMergeOrder=2, maxMergeOrder=4, 
                      minModuleSize=10, maxModuleSize=50, 
                      useAbs=TRUE, useWeighted=TRUE,
                      verbose=FALSE)
  saveRipmResultsToRdata(modListRipm, "foobar")
  modListChk <- loadRipmResultsFromRdata("foobar_results.rds")
  file.remove(c("foobar_results.rds"))
  expect_equal(modListChk, modListRipm)
})
