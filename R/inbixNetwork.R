# ----------------------------------------------------------------------------
# inbixNetwork.R - Bill White - 10/10/15
#
# Rinbix package network/graph functions.

# ----------------------------------------------------------------------------
#' Get stats about the IGraph g.
#' 
#' \code{getIGraphStats}
#' 
#' @family network functions
#' @param g IGraph graph object.
#' @return \code{data.frame} of graph statistics.
#' @examples
#' require(igraph)
#' g <- erdos.renyi.game(1000, 1/1000)
#' igstats <- getIGraphStats(g)
#' @export
getIGraphStats <- function(g) {
  return(
    data.frame(
      nnodes=igraph::vcount(g), 
      maxDegree=max(igraph::degree(g)),
      minDegree=min(igraph::degree(g)),
      nedges=igraph::ecount(g), 
      radius=igraph::radius(g), 
      density=igraph::graph.density(g), 
      transitivity=igraph::transitivity(g), 
      avgPathLen=igraph::average.path.length(g), 
      assortativeDegree=igraph::assortativity.degree(g)))
}

# ----------------------------------------------------------------------------
#' Prepare and adjacency matrix for network analysis.
#' 
#' \code{prepareAdjacencyMatrix}
#' 
#' @family network functions
#' @param adjacencyMatrix \code{matrix} adjacency matrix
#' @param thresholdType \code{string} threshold type: "hard" or "soft".
#' @param thresholdValue \code{numeric} hard threshold correlation value or soft threshold power.
#' @param useAbs \code{logical} take absolute value of the correlation matrix.
#' @param useWeighted \code{logical} use weighted adjacency matrix versus binary.
#' @param verbose \code{logical} to send runtime messages to stdout.
#' @return \code{matrix} transformed for network analysis.
#' @export
prepareAdjacencyMatrix <- function(adjacencyMatrix=NULL, 
                                   thresholdType="hard", 
                                   thresholdValue=0.8, 
                                   useAbs=TRUE, 
                                   useWeighted=FALSE,
                                   verbose=FALSE) {
  if(is.null(adjacencyMatrix)) {
    stop("prepareAdjacencyMatrix: adjacencyMatrix is a required parameter")
  }
  if(verbose) cat("Begin prepareAdjacencyMatrix\n")
  if(useAbs) {
    if(verbose) cat("Taking absolute value of adjacencyMatrix\n")
    adjacencyMatrix <- abs(adjacencyMatrix)
  }
  # thresholding
  if(thresholdType == "soft") {
    if(verbose) cat("Soft threshold adjacencyMatrix^", thresholdValue, "\n", sep="")
    adjacencyMatrix <- adjacencyMatrix ^ thresholdValue
  } else {
    if(verbose) cat("Hard threshold abs(adjacencyMatrix) >", thresholdValue, "\n")
    passThreshold <- adjacencyMatrix > thresholdValue
    # create adjacency matrix from thresholding
    if(useWeighted) {
      if(verbose) cat("Keeping weighted values that pass threshold [", thresholdValue, "]\n")
      adjacencyMatrix[!passThreshold] <- 0
    } else {
      if(verbose) cat("Keeping values that pass threshold [", thresholdValue, "] as binary 0/1\n")
      adjacencyMatrix[passThreshold] <- 1
      adjacencyMatrix[!passThreshold] <- 0
    }
  }
  if(verbose) cat("End prepareAdjacencyMatrix\n")
  adjacencyMatrix
}

# ----------------------------------------------------------------------------
#' Print stats about the IGraph g.
#' 
#' \code{printIGraphStats}
#' 
#' @family network functions
#' @param g IGraph graph object.
#' @export
printIGraphStats <- function(g) {
  cat("Number of nodes:", igraph::vcount(g), "\n")
  cat("Maximum degree: ", max(igraph::degree(g)), "\n")
  cat("Minimum degree: ", min(igraph::degree(g)), "\n")
  cat("Number of edges:", igraph::ecount(g), "\n")
  cat("Radius:         ", igraph::radius(g), "\n")
  cat("Density:        ", igraph::graph.density(g), "\n")
  cat("Transitivity:   ", igraph::transitivity(g), "\n")
  cat("Avg path length:", igraph::average.path.length(g), "\n")
  cat("Assortativity:  ", igraph::assortativity.degree(g), "\n")
  
  return(TRUE)
}

# ----------------------------------------------------------------------------
#' Random network simulation.
#'
#' Generates a Erdos-Renyi random network as an adjacency matrix.
#' 
#' \code{randomNetworkSim}
#' 
#' @family network functions
#' @family simulation functions
#' @param n \code{numeric} number of nodes.
#' @param p \code{numeric} probability of attachment.
#' @param doFitPlot \code{logical} plot sclae-free fit.
#' @param doNetworkPlot \code{logical} plot the network.
#' @param doHistPlot \code{logical} plot histogram of node degree.
#' @param useIgraph \code{logical} use the IGraph library to generate the graph.
#' @param numBins \code{numeric} number of bins to use in the degree histogram.
#' @param filePrefix \code{string} file prefix for output files.
#' @param verbose \code{logical} to send runtime messages to stdout.
#' @return network adjacency \code{matrix}.
#' @examples
#' net <- randomNetworkSim(n=1000)
#' @export
randomNetworkSim <- function(n=100, p=0.1, doFitPlot=F, doNetworkPlot=F, doHistPlot=F, 
                             useIgraph=FALSE, numBins=10, filePrefix="random_network_sim",
                             verbose=FALSE) {
  # Erdos-Renyi
  # usage: A <- random_network_sim(100,.1,1)
  if(useIgraph) {
    g <- igraph::erdos.renyi.game(n, p)
    printIGraphStats(g)
    A <- igraph::get.adjacency(g)
  } else {
    A <- matrix(ncol=n, nrow=n, data=runif(n * n))
    # undirected no self-connections, no weights
    A <- as.matrix(Matrix::triu(A < p, 1)) # take upper triangle (k=1) of (random nxn < p)
    #A = sparse(A)
    # now make A symmetric, undirected
    A <- A + t(A)
  }
  degrees <- rowSums(A)
  
  # calculates the scale free parameter and shows network
  if(useIgraph) {
    deg_counts <- floor(igraph::degree.distribution(g)*100)
  } else {
    k_rows <- rowSums(A)
    bins <- numBins
    histObj <- hist(k_rows, bins, plot=doHistPlot)
    deg_counts <- histObj$counts
  }
  #print(deg_counts)
  # find non zero entries for log plot  
  nz_idx <- which(deg_counts != 0, arr.ind=T)
  nz_deg_counts <- deg_counts[nz_idx]
  
  # Igraph's discrete power law fit
  powerLawFit <- igraph::power.law.fit(nz_deg_counts)
  if(verbose) cat("Erdos-Renyi Graph Estimate of x_min:", powerLawFit$xmin, "\n")
  if(verbose) cat("Erdos-Renyi Graph Estimate of x_min KS value:", powerLawFit$KS.stat, "\n")
  if(verbose) cat("Erdos-Renyi Graph Estimate of scaling parameter alpha:", powerLawFit$alpha, "\n")
  if(doFitPlot) {
    png(paste(filePrefix, "_er_powerlaw_fit.png", sep=""), width=1024, height=768)
    # data
    plot(powerLawFit, xlab="Degree", ylab="CDF", 
         main="Erdos-Renyi Degree Distribution with Power Law Fit")
    # fit
    lines(powerLawFit, col=2)
    dev.off()
  }
  
  # viz with lower triangle removed.
  if(doNetworkPlot) {
    g <- igraph::graph.adjacency(A)
    igraph::V(g)$size <- scaleAB(degrees, 10, 20)
    png(paste(filePrefix, "_er_network.png", sep=""), width=1024, height=768)
    plot(g, layout=igraph::layout.fruchterman.reingold, edge.arrow.mode=0)
    dev.off()
  }
  
  return(A)
}

# ----------------------------------------------------------------------------
#' Scale free network simulation.
#' 
#' Generates a scale free network as an adjacency matrix.
#' 
#' \code{scaleFreeSim}
#' 
#' @family network functions
#' @family simulation functions
#' @param n \code{numeric} number of nodes.
#' @param doFitPlot \code{logical} plot sclae-free fit.
#' @param doNetworkPlot \code{logical} plot the network.
#' @param useIgraph \code{logical} use the IGraph library to generate the graph.
#' @param numBins \code{numeric} number of bins to use in the degree histogram.
#' @param filePrefix \code{string} file prefix for output files.
#' @param verbose \code{logical} to send runtime messages to stdout.
#' @return network adjacency \code{matrix}.
#' @examples
#' net <- scaleFreeSim(n=1000)
#' @export
scaleFreeSim <- function(n=100, doFitPlot=F, doNetworkPlot=F, useIgraph=F, 
                         numBins=10, filePrefix="scale_free_sim", verbose=FALSE) {
  
  # Baraba'si-Albert (BA) model generation of scale-free network
  if(useIgraph) {
    g <- igraph::barabasi.game(n, directed=F)
    foo <- printIGraphStats(g)
    A <- igraph::get.adjacency(g)
  } else {
    A <- matrix(ncol=n, nrow=n, data=c(0))
    deg <- rep(0, n)
    A[1,2] <- 1      # Connect nodes 1 and 2 to seed the network.
    A[2,1] <- 1      # Make symmetric
    deg[1] <- 1      
    deg[2] <- 1      
    sumdeg <- 2
    # Other nodes will be added one at a time and given a chance to connect to
    # existing nodes. Existing nodes will get more chances to make a connection.
    for(j in 3:n) {        # Start adding the remaining nodes j
      for(i in 1:(j-1)) {  # Give j a chance to connect to existing nodes 1:j-1
        p_i <- deg[i] / sumdeg
        if(runif(1) < p_i) {
          A[i, j] <- 1
          A[j, i] <- 1
          deg[i] <- deg[i] + 1
          deg[j] <- deg[j] + 1
          sumdeg <- sumdeg + 2
        }
      }
    }
  }
  degrees <- rowSums(A)
  
  # calculates the scale free parameter and shows network
  if(useIgraph) {
    deg_counts <- floor(igraph::degree.distribution(g)*100)
  } else {
    k_rows <- rowSums(A)
    bins <- numBins
    histObj <- hist(k_rows, bins, plot=F)
    deg_counts <- histObj$counts
  }

  # find non zero entries for log plot  
  nz_idx <- which(deg_counts != 0, arr.ind=T)
  nz_deg_counts <- deg_counts[nz_idx]
  
  # Igraph's discrete power law fit
  powerLawFit <- igraph::power.law.fit(nz_deg_counts)
  if(verbose) cat("Scale Free Graph Estimate of x_min:", powerLawFit$xmin, "\n")
  if(verbose) cat("Scale Free Graph Estimate of x_min KS value:", powerLawFit$KS.stat, "\n")
  if(verbose) cat("Scale Free Graph Estimate of scaling parameter alpha:", powerLawFit$alpha, "\n")
  if(doFitPlot) {
    png(paste(filePrefix, "_ba_powerlaw_fit.png", sep=""), width=1024, height=768)
    # data
    plot(powerLawFit, xlab="Degree", ylab="CDF", 
         main="Scale Free Graph Degree Distribution with Power Law Fit")
    # fit
    lines(powerLawFit, col=2)
    dev.off()
  }
  
  # viz with lower triangle removed
  #A[upper.tri(A)] <- 0
  if(doNetworkPlot) {
    g <- igraph::graph.adjacency(A)
    igraph::V(g)$size <- scaleAB(degrees, 10, 20)
    png(paste(filePrefix, "_ba_network.png", sep=""), width=1024, height=768)
    plot(g, layout=igraph::layout.fruchterman.reingold, edge.arrow.mode=0)
    dev.off()
  }
  
  return(A)
}
