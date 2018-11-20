#' Rank variables by EpistasisRank algorithm.
#' 
#' This function replaces the original Rinbix snprank().
#'
#' \code{EpistasisRank}
#' 
#' @family feature selection functions
#' @param G \code{matrix} genetic association interaction network.
#' @param Gamma_vec \code{numeric} gamma vector, either a constant value or a vector
#' @return sortedTable \code{data.frame} with variable, EpistasisRank, 
#' diagonal and degree columns.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' rinbixER_DF <- EpistasisRank(rinbixRegain)
#' @export
EpistasisRank <- function(G = NULL, Gamma_vec = 0.85)
{
  n <- nrow(G)
  geneNames <- colnames(G)
  Gdiag <- diag(G)
  Gtrace <- sum(Gdiag)
  #colsum <- colSums(G)
  diag(G) <- 0
  Gtrace <- Gtrace * n
  colsumG <- colSums(G)
  #rowSumG <- rowSums(G)
  rowsum_denom <- matrix(0, n, 1)
  for (i in 1:n) {
    localSum <- 0
    for (j in 1:n) {
      factor <- ifelse(G[i, j] != 0, 1, 0)
      localSum <- localSum + (factor * colsumG[j])
    }
    rowsum_denom[i] <- localSum
  }
  if (length(Gamma_vec) == 1) {
    gamma_vec <- rep(Gamma_vec, n)
  } else {
    gamma_vec <- Gamma_vec
  }
  gamma_matrix <- matrix(nrow = n, ncol = n, data = rep(gamma_vec, n))
  if (Gtrace) {
    b <- ((1 - gamma_vec) / n) + (Gdiag / Gtrace)
  }
  else {
    b <- ((1 - gamma_vec) / n)
  }
  D <- matrix(nrow = n, ncol = n, data = c(0))
  diag(D) <- 1 / colsumG
  I <- diag(n)
  temp <- I - gamma_matrix * G %*% D
  r <- solve(temp, b)
  ERank <- r / sum(r)
  saveTable <- data.frame(gene = geneNames, ER = ERank)
  saveTable[order(saveTable$ER, decreasing = TRUE), ]
}

#' Rank variables by EpistasisKatz algorithm.
#' 
#' This function can be use as original Katz algorithm or as EpistasisKatz
#' that incorporates prior knowledge.
#'
#' \code{EpistasisKatz}
#' 
#' @family feature selection functions
#' @param A \code{matrix} network of features either in adjacency or interaction format.
#' @param alpha \code{numeric} a vector with numeric values.
#' @param beta \code{numeric} either a vector of constant values or prior knowledge.
#' @return features ranking \code{vector} with features name.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' alpha <- 1 / mean(colSums(rinbixRegain))
#' beta <- diag(rinbixRegain)
#' rinbixEK_DF <- EpistasisKatz(rinbixRegain, alpha, beta)
#' @export
EpistasisKatz <- function(A = NULL, alpha = NULL, beta = NULL) {
  n <- nrow(A)
  I <- diag(1, n); 
  c <- solve(I - alpha * A,beta)  
  return(as.vector(c))
}

#' Rank variables by PageRank algorithm.
#' 
#' This function is the original PageRank algorithm that incorporates prior knowledge.
#'
#' \code{PageRank}
#' 
#' @family feature selection functions
#' @param A \code{matrix} network of features in adjacency format.
#' @param damping \code{numeric} damping factor, either a constant number or a vector of prior knowledge.
#' @return features ranking \code{vector} with features name.
#' @examples
#' data(testdata10)
#' adj <- ifelse((cor(testdata10))>0, 1, 0)
#' damping <- 0.85
#' myPageRank <- PageRank(adj, damping)
#' @export
PageRank <- function(A = NULL, damping = 0.85){
  #make sure A is matrix type
  n <-nrow(A)
  rand.jump.col <- (1-damping)/n*matrix(rep(1, n),nrow=n)
  gDegs <- sapply(1:n, function(x) sum(A[x,]))
  Dk_cols<-matrix(rep(gDegs,n), nrow=n, byrow=TRUE)
  Amarkov <- A/Dk_cols
  my.sol <- solve(diag(n) - damping*Amarkov, rand.jump.col)
  return(my.sol)
}