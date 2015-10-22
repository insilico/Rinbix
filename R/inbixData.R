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

#' Simulated RNASeq data set with 40 subjects and 100 genes.
#'
#' A dataset containing  with 40 subjects (20 cases, 20 controls), 100 genes, 
#' 10 differentially expressed. Training and testing data with 40 responses.
#'
#' @docType data
#' @keywords datasets
#' @name simrnaseq
#' @usage data(simrnaseq)
#' @format Two data frames with 40 rows and 100 columns. Two vectors of 40 responses.
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
