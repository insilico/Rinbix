#' Rank variables by EpistasisRank algorithm.
#' 
#' This function replaces the original Rinbix snprank().
#'
#' \code{EpistasisRank}
#' 
#' @family feature selection functions
#' @param G \code{matrix} genetic association interaction network.
#' @param Gamma_vec \code{numeric} gamma vector
#' @return sortedTable \code{data.frame} with variable, EpistasisRank, 
#' diagonal and degree columns.
#' @examples
#' data(testdata10)
#' rinbixRegain <- regainParallel(testdata10, stdBetas=TRUE, absBetas=TRUE)
#' rinbixSnpranksDF <- EpistasisRank(rinbixRegain)
#' @export
EpistasisRank <- function(G, Gamma_vec)
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
  snpranks <- r / sum(r)
  saveTable <- data.frame(gene = geneNames, snprank = snpranks)
  saveTable[order(saveTable$snprank, decreasing = TRUE), ]
}