# ----------------------------------------------------------------------------
# inbixSim.R - Bill White - 10/10/15
#
# Rinbix package data simulation functions.

# ----------------------------------------------------------------------------
#' Create a random regression data set with binary class.
#' 
#' \code{createRandomRegressionDataset} 
#' 
#' @param numRows Numeric number of rows (samples).
#' @param numCols Numeric number of columns (independent variables),
#' @return data frame with class column.
#' @examples
#' ds <- createRandomRegressionDataset(100, 100)
#' @export
createRandomRegressionDataset <- function(numRows, numCols) {
  dmatrix <- matrix(nrow=numRows, ncol=numCols, data=rnorm(numRows*numCols))
  dpheno <- c(rep(0, numRows/2), rep(1, numRows/2))
  dataset <- cbind(dmatrix, dpheno)
  colnames(dataset) <- c(paste("var", 1:numCols, sep=""), "Class")
  as.data.frame(dataset)
}

# ----------------------------------------------------------------------------
#' Create a random matrix for differential co-expression simulation.
#' 
#' \code{createRandomMatrix} 
#' 
#' @param M Numeric number of genes.
#' @param N Numeric number of subjects.
#' @param meanExpression Numeric mean gene expression.
#' @param randSdNoise Numeric standard deviation of random normal (rnorm) noise.
#' @return list with subject by gene data frame with class column.
#' @examples
#' ds <- createRandomMatrix(100, 100, 7, 0.05)$regressionData
#' @export
createRandomMatrix <- function(M, N, meanExpression, randSdNoise) {
  # create a random data matrix
  D <- matrix(nrow=M, ncol=N, data=rnorm(M*N, mean=meanExpression, sd=randSdNoise))
  # return a regression ready data frame
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("ctrl", 1:n1, sep=""), paste("case", 1:n2, sep=""))
  phenos <- c(rep(0, n1), rep(1, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("gene", sprintf("%04d", 1:M), sep=""), "Class")
  rownames(newD) <- subIds
  list(regressionData=newD)
}

# ----------------------------------------------------------------------------
#' Create a differentially coexpressed data set with interactions and main effects,
#' 
#' \code{createDiffCoexpMatrix} 
#' 
#' @param M genes
#' @param N subjects
#' @param meanExpression mean gene expression
#' @param A network model adjacency matrix
#' @param randSdNoise standard deviation of random normal (rnorm) noise
#' @param sdNoise standard deviation of noise in differential correlation
#' @param mGenesToPerturb number of genes to perturb
#' @param sampleIndicesInteraction interaction gene indices
#' @param sampleIndicesMainEffects main effects indices
#' @param mainEffectMode 1=use interaction indices, else use main effects indices
#' @param mainEffect desired fold change
#' @param verbose verbose output to stdout
#' @return list with subject by gene data frame with class column and fold changes.
#' @examples
#' data("scaleFreeNetwork")
#' dsobj <- createDiffCoexpMatrix(M=100, 
#'                              N=100, 
#'                              meanExpression=7, 
#'                              A=scaleFreeNetwork, 
#'                              randSdNoise=0.05, 
#'                              sdNoise=1.5, 
#'                              mGenesToPerturb=3,
#'                              sampleIndicesMainEffects = c(5, 10, 15),
#'                              sampleIndicesInteraction = c(5, 10, 15),
#'                              mainEffectMode=1,
#'                              mainEffect=4,
#'                              verbose=FALSE)
#' ds <- dsobj$regressionData  
#' @export
createDiffCoexpMatrix <- function(M, N, meanExpression, A, 
                                  randSdNoise, sdNoise, mGenesToPerturb, 
                                  sampleIndicesInteraction, sampleIndicesMainEffects, 
                                  mainEffectMode, mainEffect, verbose=FALSE) {
  # create a random data matrix
  D <- matrix(nrow=M, ncol=N, data=rnorm(M*N, mean=meanExpression, sd=randSdNoise))
  # add co-expression
  already_modified <- rep(0, M)
  already_modified[1] <- 1 
  for(i in 1:(M-1)) {
    for(j in (i+1):M) {
      #cat("Condidering A: row", i, "column", j, "\n")
      if((A[i, j] == 1) && (!already_modified[j])) {
        #cat("Making row", j, "from row", i, "\n")
        D[j, ] <- D[i, ] + rnorm(N, mean=0, sd=as.numeric(sdNoise))
        already_modified[j] <- 1
      } else {
        if(already_modified[j]==1 && !already_modified[i]) {
          # if j is already modified, we want to modify i, 
          # unless i is already modified then do nothing 
          D[i,] <- D[j,] + rnorm(N, mean=0, sd=as.numeric(sdNoise))
        }
      }
    }
  }
  
  # perturb to get differential coexpression
  n1 <- N / 2;
  foldChangesIntrBefore <- NULL
  foldChangesMainBefore <- NULL
  foldChangesIntrAfter <- NULL
  foldChangesMainAfter <- NULL
  for(i in 1:mGenesToPerturb) { 
    # get the group 2 gene expression and randomly order
    geneIdxInteraction <- sampleIndicesInteraction[i]
    if(mainEffectMode == 1) {
      geneIdxMainEffects <- geneIdxInteraction
    } else {
      geneIdxMainEffects <- sampleIndicesMainEffects[i]
    }
    # fold change before adding main effect
    g0 <- D[geneIdxMainEffects, (n1 + 1):N]
    g1 <- D[geneIdxMainEffects, 1:n1]
    fcInteractionBefore <- getFoldChange(D, geneIdxInteraction, n1, N)
    fcMainEffectsBefore <- getFoldChange(D, geneIdxMainEffects, n1, N)
    foldChangesIntrBefore <- c(foldChangesIntrBefore, fcInteractionBefore)
    foldChangesMainBefore <- c(foldChangesMainBefore, fcMainEffectsBefore)
    
    # calculate the amount to add to each value in the affected group
    amt_to_add_per <- 0
    if(mainEffect > 1) {
      mu_g0 <- mean(g0)
      mu_g1 <- mean(g1)
      target_fc <- mainEffect
      new_mu_g1 <- target_fc * mu_g0
      new_sum_g1 <- new_mu_g1 * n1
      sum_to_add <- new_sum_g1 - sum(g1)
      amt_to_add_per <- sum_to_add / n1
    }
    if(verbose) {
      cat("To achieve a main effects fold change of", mainEffect, 
          "adding", amt_to_add_per, "to each value\n")
    }
    
    # differential coexpression + main effect
    x <- D[geneIdxInteraction, 1:n1]
    x <- x[order(runif(length(x)))]
    # add a main effect to either the same as interaction gene or another gene
    # main effect is the desired fold change, so calculate the amount to add to 
    # create the desired fold change
    if(mainEffectMode == 1) {
      D[geneIdxInteraction, 1:n1] <- x + amt_to_add_per
    } else {
      D[geneIdxInteraction, 1:n1] <- x
      z <- D[geneIdxMainEffects, 1:n1]
      D[geneIdxMainEffects, 1:n1] <-  z + amt_to_add_per
    }
    # fold change after adding main effect
    fcInteractionAfter <- getFoldChange(D, geneIdxInteraction, n1, N)
    fcMainEffectsAfter <- getFoldChange(D, geneIdxMainEffects, n1, N)
    foldChangesIntrAfter <- c(foldChangesIntrAfter, fcInteractionAfter)
    foldChangesMainAfter <- c(foldChangesMainAfter, fcMainEffectsAfter)
    if(verbose) {
      cat("--------------------------------------------------\n")
      cat("main effects gene index: ", geneIdxMainEffects, 
          ", FC before: ", fcMainEffectsBefore, 
          ", after: ", fcMainEffectsAfter, "\n", sep="")
      cat("interaction gene index: ", geneIdxInteraction, 
          ", FC before: ", fcInteractionBefore, 
          ", after: ", fcInteractionAfter, "\n", sep="")
    }
  }
  
  # return a regression ready data frame
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("case", 1:n1, sep=""), paste("ctrl", 1:n2, sep=""))
  phenos <- c(rep(1, n1), rep(0, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("gene", sprintf("%04d", 1:M), sep=""), "Class")
  rownames(newD) <- subIds
  list(regressionData=newD,
       fcIntrBefore=foldChangesIntrBefore,
       fcMainBefore=foldChangesMainBefore,
       fcIntrAfter=foldChangesIntrAfter,
       fcMainAfter=foldChangesMainAfter
  )
}

# ----------------------------------------------------------------------------
#' Create a differentially coexpressed data set without main effects.
#' 
#' \code{createDiffCoexpMatrixNoME} 
#' 
#' @param M Numeric number of genes.
#' @param N Numeric number of subjects.
#' @param meanExpression Numeric mean gene expression.
#' @param A Matrix network adjacency matrix.
#' @param randSdNoise Numeric standard deviation of random normal (rnorm) noise.
#' @param sdNoise Numeric standard deviation of noise in differential correlation.
#' @param sampleIndicesInteraction Vector interaction gene indices.
#' @return list with subject by gene data frame with class column.
#' @examples
#' data("scaleFreeNetwork")
#' dsobj <- createDiffCoexpMatrixNoME(M=100, 
#'                              N=100, 
#'                              meanExpression=7, 
#'                              A=scaleFreeNetwork, 
#'                              randSdNoise=0.05, 
#'                              sdNoise=1.5, 
#'                              sampleIndicesInteraction = c(5, 10, 15))
#' ds <- dsobj$regressionData  
#' @export
createDiffCoexpMatrixNoME <- function(M, N, meanExpression, A, randSdNoise, 
  sdNoise, sampleIndicesInteraction) {
  # create a random data matrix
  D <- matrix(nrow=M, ncol=N, data=rnorm(M*N, mean=meanExpression, sd=randSdNoise))
  
  # add co-expression
  already_modified <- rep(0, M)
  already_modified[1] <- 1 
  for(i in 1:(M-1)) {
    for(j in (i+1):M) {
      #cat("Condidering A: row", i, "column", j, "\n")
      if((A[i, j] == 1) && (!already_modified[j])) {
        #cat("Making row", j, "from row", i, "\n")
        D[j, ] <- D[i, ] + rnorm(N, mean=0, sd=as.numeric(sdNoise))
        already_modified[j] <- 1
      } else {
        if(already_modified[j]==1 && !already_modified[i]) {
          # if j is already modified, we want to modify i, 
          # unless i is already modified then do nothing 
          D[i,] <- D[j,] + rnorm(N, mean=0, sd=sdNoise)
        }
      }
    }
  }
  
  # perturb to get differential coexpression
  n1 <- N / 2;
  mGenesToPerturb <- length(sampleIndicesInteraction)
  for(i in 1:mGenesToPerturb) { 
    geneIdxInteraction <- sampleIndicesInteraction[i]

    g0 <- D[sampleIndicesInteraction, (n1 + 1):N]
    g1 <- D[sampleIndicesInteraction, 1:n1]
    
    # get the group 2 gene expression and randomly order for differential coexpression
    x <- D[geneIdxInteraction, 1:n1]
    x <- x[order(runif(length(x)))]
    D[geneIdxInteraction, 1:n1] <- x
  }

  # return a regression ready data frame
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("case", 1:n1, sep=""), paste("ctrl", 1:n2, sep=""))
  phenos <- c(rep(1, n1), rep(0, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("gene", sprintf("%04d", 1:M), sep=""), "Class")
  rownames(newD) <- subIds
  list(regressionData=newD)
}

# ----------------------------------------------------------------------------
#' Create a differentially coexpressed data set without any association.
#' 
#' \code{createDiffCoexpMatrixNull} 
#' 
#' @param M Numeric number of genes.
#' @param N Numeric number of subjects.
#' @param meanExpression Numeric mean gene expression.
#' @param A Matrix network adjacency matrix.
#' @param randSdNoise Numeric standard deviation of random normal (rnorm) noise.
#' @param sdNoise Numeric standard deviation of noise in differential correlation.
#' @return list with subject by gene data frame with class column.
#' @examples
#' data("scaleFreeNetwork")
#' dsobj <- createDiffCoexpMatrixNull(M=100, 
#'                              N=100, 
#'                              meanExpression=7, 
#'                              A=scaleFreeNetwork, 
#'                              randSdNoise=0.05, 
#'                              sdNoise=1.5)
#' ds <- dsobj$regressionData  
#' @export
createDiffCoexpMatrixNull <- function(M, N, meanExpression, A, randSdNoise, sdNoise) {
  # create a random data matrix
  D <- matrix(nrow=M, ncol=N, data=rnorm(M*N, mean=meanExpression, sd=randSdNoise))
  
  # add co-expression
  already_modified <- rep(0, M)
  already_modified[1] <- 1 
  for(i in 1:(M-1)) {
    for(j in (i+1):M) {
      #cat("Condidering A: row", i, "column", j, "\n")
      if((A[i, j] == 1) && (!already_modified[j])) {
        #cat("Making row", j, "from row", i, "\n")
        D[j, ] <- D[i, ] + rnorm(N, mean=0, sd=as.numeric(sdNoise))
        already_modified[j] <- 1
      } else {
        if(already_modified[j]==1 && !already_modified[i]) {
          # if j is already modified, we want to modify i, 
          # unless i is already modified then do nothing 
          D[i,] <- D[j,] + rnorm(N, mean=0, sd=sdNoise)
        }
      }
    }
  }
  
  # return a regression ready data frame
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2
  subIds <- c(paste("case", 1:n1, sep=""), paste("ctrl", 1:n2, sep=""))
  phenos <- c(rep(1, n1), rep(0, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("gene", sprintf("%04d", 1:M), sep=""), "Class")
  rownames(newD) <- subIds
  list(regressionData=newD)
}

# ----------------------------------------------------------------------------
#' Create a main effects-only data set.
#' 
#' \code{createMainEffectsMatrix} 
#' 
#' @param M Numeric number of genes.
#' @param N Numeric number of subjects.
#' @param meanExpression Numeric mean gene expression.
#' @param randSdNoise Numeric standard deviation of random normal (rnorm) noise.
#' @param sampleIndicesMainEffects main effects indices
#' @param mainEffect Numeric desired fold change.
#' @param doScale Flag scale the resulting matrix.
#' @param doLog Flag log2 transform resulting matrix.
#' @return list with subject by gene data frame with class column and fold changes.
#' @examples
#' dsobj <- createMainEffectsMatrix(M=100, 
#'                              N=100, 
#'                              meanExpression=7, 
#'                              randSdNoise=0.05, 
#'                              sampleIndicesMainEffects = c(5, 10),
#'                              mainEffect=4, 
#'                              doScale=FALSE, 
#'                              doLog=FALSE)
#' ds <- dsobj$regressionData  
#' @export
createMainEffectsMatrix <- function(M, N, meanExpression, randSdNoise, 
  sampleIndicesMainEffects, mainEffect, doScale=FALSE, doLog=FALSE) {

  # create a random data matrix
  D <- matrix(nrow=M, ncol=N, data=rnorm(M*N, mean=meanExpression, sd=randSdNoise))
  dimN <- ncol(D)
  n1 <- dimN / 2
  n2 <- dimN / 2

  # fold change before adding main effect
  foldChangesMainBefore <- NULL
  foldChangesMainAfter <- NULL
  for(i in 1:length(sampleIndicesMainEffects)) { 
    geneIdxMainEffects <- sampleIndicesMainEffects[i]

    # fold change before adding main effect
    g1 <- D[geneIdxMainEffects, 1:n1]
    g2 <- D[geneIdxMainEffects, (n1 + 1):N]
    fcMainEffectsBefore <- getFoldChangeME(D, geneIdxMainEffects, n1, N)
    foldChangesMainBefore <- c(foldChangesMainBefore, fcMainEffectsBefore)
    
    # calculate the amount to add to each value in the affected group
    amt_to_add_per <- 0
    mu_g1 <- mean(g1)
    mu_g2 <- mean(g2)
    target_fc <- mainEffect
    new_mu_g2 <- target_fc * mu_g1
    new_sum_g2 <- (new_mu_g2 * n2) - sum(g2)
    amt_to_add_per <- new_sum_g2 / n2
    #cat("To achieve a main effects fold change of", mainEffect, 
    #    "adding", amt_to_add_per, "to each value\n")
    
    # differential coexpression + main effect
    x <- D[geneIdxMainEffects, (n1+1):N]
    D[geneIdxMainEffects, (n1+1):N] <-  x + amt_to_add_per
    
    # fold change after adding main effect
    fcMainEffectsAfter <- getFoldChangeME(D, geneIdxMainEffects, n1, N)
    foldChangesMainAfter <- c(foldChangesMainAfter, fcMainEffectsAfter)
    # cat("main effects gene index: ", geneIdxMainEffects, 
    #     ", FC before: ", fcMainEffectsBefore, 
    #     ", after: ", fcMainEffectsAfter, "\n", sep="")
    # cat("--------------------------------------------------\n")
  }

  # this handles the case of large sd noise creating negative expression
  D[D <= 0] <- 0.1

  if(doLog) {
    D <- log(D)
  }

  if(doScale) {
    D <- scale(D)
  }

  # return a regression ready data frame
  subIds <- c(paste("ctrl", 1:n1, sep=""), paste("case", 1:n2, sep=""))
  phenos <- c(rep(0, n1), rep(1, n2))
  newD <- cbind(t(D), phenos)
  colnames(newD) <- c(paste("gene", sprintf("%04d", 1:M), sep=""), "Class")
  rownames(newD) <- subIds
  list(regressionData=newD,
       fcMainBefore=fcMainEffectsBefore,
       fcMainAfter=fcMainEffectsAfter)
}

# ----------------------------------------------------------------------------
#' Compute interaction fold change between groups for gene at idx.
#' 
#' \code{getFoldChange} 
#' 
#' @param D Matrix data matrix.
#' @param idx Numeric index of variable to get fold change.
#' @param n1 Numeric size of group 1.
#' @param N Numeric number of subjects.
#' @return numeric fold change.
#' @examples
#' dsobj <- createMainEffectsMatrix(M=100, 
#'                              N=100, 
#'                              meanExpression=7, 
#'                              randSdNoise=0.05, 
#'                              sampleIndicesMainEffects = c(5, 10),
#'                              mainEffect=4, 
#'                              doScale=FALSE, 
#'                              doLog=FALSE)
#' ds <- dsobj$regressionData  
#' fc <- getFoldChange(t(ds), 5, 50, 100)
#' @export
getFoldChange <- function(D, idx, n1, N) {
  # fold change before changes
  g0 <- D[idx, 1:n1]
  g1 <- D[idx, (n1 + 1):N]
  fc <- mean(g1) / mean(g0)
}
getFoldChangeME <- function(D, idx, n1, N) {
  getFoldChange(D, idx, n1, N)
}

# -----------------------------------------------------------------------------
#' Create a simulated correlation matrix with block diagonal clusters.
#' 
#' \code{simCorrMatrix} 
#' 
#' @param n Numeric dim of matrix, number of genes.
#' @param num_clust Numeric number of clusters.
#' @param max_noise_corr Numeric background noise correlation.
#' @param lower_true_corr Numeric minimum strength of correlation within clusters
#' lower limit for a "true" correlation.
#' @return n x n correlation matrix.
#' @examples
#' mat <- simCorrMatrix(n=400, num_clust=20, max_noise_corr=0.5, lower_true_corr=0.5)
#' @export
simCorrMatrix <- function(n=4000, num_clust=20, max_noise_corr=0.2, lower_true_corr=0.6) {
  # creates same-size clusters, choose num_clust so that it divides evenly into n
  if(n %% num_clust != 0) {
    stop("sim_corr_matrix: Number of clusters doesn't evenly divide number of genes")
  } else { ######## the rest of code is in this else statement
    sim_matrix <- matrix(runif(n*n,0,max_noise_corr),nrow=n)  # background noise
    ######################## create matrix #####################
    clust_id_names <- integer(n) # initialize
    block_dim <- n / num_clust # size of the each cluster
    clust_id <- 1
    clust_id_names <- integer(n)
    for(i in seq(0, n - block_dim, block_dim)) {
      block_diag<-matrix(runif(block_dim * block_dim, lower_true_corr, 1), nrow=block_dim)
      sim_matrix[(i + 1):(block_dim + i), (i + 1):(block_dim + i)] <- block_diag
      clust_id_names[(i + 1):(block_dim + i)] <- rep(clust_id, block_dim)
#       clust_id_names[(i + 1):(block_dim + i)] <- paste("var", 1:block_dim, 
#                                                        ".cls", clust_id, 
#                                                        ".blk", block_dim, sep="")
      clust_id <- clust_id + 1
    }
    # print matrix
    # format(sim_matrix, scientific=FALSE, digits=2)
  }
  sim_matrix[lower.tri(sim_matrix)] = t(sim_matrix)[lower.tri(sim_matrix)]
  rownames(sim_matrix) <- clust_id_names
  colnames(sim_matrix) <- clust_id_names
  
  sim_matrix
}

# -----------------------------------------------------------------------------
#' Create a simulated correlation matrix with non-uniform block diagonal clusters.
#' 
#' \code{simCorrMatrixNonUniform} 
#' 
#' @param n Numeric dim of matrix, number of genes.
#' @param num_clust Numeric number of clusters.
#' @param max_noise_corr Numeric background noise correlation.
#' @param lower_true_corr Numeric minimum strength of correlation within clusters
#' lower limit for a "true" correlation.
#' @return n x n correlation matrix
#' @examples
#' mat <- simCorrMatrixNonUniform(n=400, num_clust=20, max_noise_corr=0.5, lower_true_corr=0.5)
#' @export
simCorrMatrixNonUniform <- function(n=400, num_clust=20, max_noise_corr=0.2, 
                                    lower_true_corr=0.6) {
  ## create a simulated correlation matrix with block diagonal clusters.
  # right now creates same-size clusters, choose num_clust so that it divides evenly into n
  ############### input ##############################
  #n <- 400 # dim of matrix, number of genes
  # make two same-size clusters
  #num_clust <- 20
  if(n %% num_clust != 0) {
    cat("ERROR in sim_corr_matrix_non_uniform: Number of clusters doesn't evenly divide number of genes.\n")
  } else { ######## the rest of code is in this else statement
    # background noise correlation    
    #max_noise_corr <-.8  # .2 #more noise
    sim_matrix <- matrix(runif(n * n, 0, max_noise_corr), nrow=n)  # background noise
    # minimum strength of correlation within clusters
    #lower_true_corr <- .2  # lower limit for a "true" correlation
    
    ######################## create matrix #####################
    unif_block_dim <- n / num_clust # start size of the each cluster
    clust_id <- 1
    clust_id_names <- integer(n) # initialize
    
    # make the number of clusters stay the same but merge and split
    
    # keep track of size of each block
    # case of n=400, num_clust = 20
    #     size 80         size 40
    #     index 80        index 120         
    cluster_index_vec <- c(0,unif_block_dim*4,unif_block_dim*6,  # merge 4 times then 2 times = 6 times
                           #  the rest (below) are uniform size 20, except for the last 3, which I want to split smaller
                           #  the number 3 is chosen because I merged 6 earliers, so I need to split 6/2 times to get
                           #  the original number of clusters
                           # index 140 to 340
                           seq(unif_block_dim*7,unif_block_dim*(num_clust-3),unif_block_dim),
                           # now we split the last three modules in half to get 6
                           # index 350, 360, ... 390
                           seq(unif_block_dim*(num_clust-3)+unif_block_dim/2,n-unif_block_dim/2,unif_block_dim/2))
    # vector of the sizes of the each cluster
    block_dim_vec <- c(unif_block_dim*4,unif_block_dim*2,rep(unif_block_dim,11),rep(unif_block_dim/2,6))
    
    for (i in 1:(length(cluster_index_vec))){
      #cat(i, cluster_index_vec[i],block_dim_vec[i],"\n")
      block_diag<-matrix(runif(block_dim_vec[i]*block_dim_vec[i],lower_true_corr,1),nrow=block_dim_vec[i])
      sim_matrix[(cluster_index_vec[i]+1):(cluster_index_vec[i]+block_dim_vec[i]),(cluster_index_vec[i]+1):(cluster_index_vec[i]+block_dim_vec[i])] <- block_diag
      clust_id_names[(cluster_index_vec[i]+1):(cluster_index_vec[i]+block_dim_vec[i])] <- rep(i,block_dim_vec[i])
    }
    rownames(sim_matrix) <- clust_id_names
    colnames(sim_matrix) <- clust_id_names
    
    # make the lower triangle equal to the upper triangle so matrix is symmetric
    sim_matrix[lower.tri(sim_matrix)] = t(sim_matrix)[lower.tri(sim_matrix)]
    
    # print matrix
    #format(sim_matrix,scientific = FALSE, digits = 2)
  }
  sim_matrix
}
