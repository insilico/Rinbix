# inbixData.R - Bill White - 10/16/15

#' A random dataset with 10 samples and 10 variables.
#'
#' A dataset containing random variables var1 ... var10. 
#' Class/phenotypes are 5 0's and 5 1's.
#'
#' @docType data
#' @keywords datasets
#' @name testdata10
#' @usage data(testdata10)
#' @format Data frame with 10 rows and 10 variables plus "Class" column.
NULL

#' A main effects data set with 100 samples and 100 variables.
#'
#' A dataset containing 100 variables var1 ... var100. 
#' Class/phenotypes are case-control coded 0/1. 
#' Indices 1 and 5 are the simulated main effects variables.
#'
#' @docType data
#' @keywords datasets
#' @name testdata100ME4
#' @usage data(testdata100ME4)
#' @format Data frame with 100 rows and 100 variables plus "Class" column.
#' @source testdata100ME4 <- createMainEffectsMatrix(M=100, N=100, meanExpression=7, 
#' randSdNoise=0.05, sampleIndicesMainEffects=c(1, 5), mainEffect=4) 
NULL

#' Simulated RNASeq data set with 40 subjects and 100 genes.
#'
#' A dataset containing  with 40 subjects (20 cases, 20 controls), 100 genes, 
#' 10 differentially expressed. Training and testing data with 40 responses.
#'
#' @docType data
#' @keywords datasets
#' @name simrnaseq
#' @usage data(simrnaseq)
NULL

#' Simulated RNASeq expression with 40 subjects and 100 genes for training.
#'
#' A dataset containing 40 subjects (20 cases, 20 controls), 100 genes, 
#' 10 differentially expressed.
#'
#' @format Data frames with 40 rows and 100 columns
#' @source \code{createDiffCoexpMatrix}
#' @name predictorsTrain
NULL

#' Simulated RNASeq expression with 40 subjects and 100 genes for testing.
#'
#' A dataset containing 40 subjects (20 cases, 20 controls), 100 genes, 
#' 10 differentially expressed.
#'
#' @format Data frames with 40 rows and 100 columns
#' @source \code{createDiffCoexpMatrix}
#' @name predictorsTest
NULL

#' Simulated RNASeq phenotypes with 40 subjects for training.
#'
#' A vector containing 40 subjects' phenotypes coded 1/2.
#'
#' @format Vector of 40 phenotypes coded 1/2.
#' @source \code{createDiffCoexpMatrix}
#' @name responseTrain
NULL

#' Simulated RNASeq phenotypes with 40 subjects for training.
#'
#' A vector containing 40 subjects' phenotypes coded 1/2.
#'
#' @format Vector of 40 phenotypes coded 1/2.
#' @source \code{createDiffCoexpMatrix}
#' @name responseTest
NULL

#' Erdos-Reyni Random Network adjacency matrix.
#'
#' A random network adjacency matrix used in data simulation.
#'
#' @docType data
#' @keywords datasets
#' @name erdosRenyiNetwork
#' @usage data("erdosRenyiNetwork")
#' @format Matrix 100x100.
NULL

#' Scale free network adjacency matrix.
#'
#' A scale-free network adjacency matrix used in data simulation.
#'
#' @docType data
#' @keywords datasets
#' @name scaleFreeNetwork
#' @usage data("scaleFreeNetwork")
#' @format Matrix 100x100.
NULL
