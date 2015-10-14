# ----------------------------------------------------------------------------
# inbixNetwork.R - Bill White - 10/10/15
#
# Rinbix package network/graph functions.

# ----------------------------------------------------------------------------
# 11/30/13 - get stats about the IGraph g
getIGraphStats <- function(g) {
  return(
    data.frame(
      nnodes=vcount(g), 
      maxDegree=max(degree(g)),
      minDegree=min(degree(g)),
      nedges=ecount(g), 
      radius=radius(g), 
      diameter=diameter(g), 
      density=graph.density(g), 
      transitivity=transitivity(g), 
      avgPathLen=average.path.length(g), 
      assortativeDegree=assortativity.degree(g)))
}

# ----------------------------------------------------------------------------
# 11/24/13 - print stats about the IGraph g
printIGraphStats <- function(g) {
  cat("Number of nodes:", vcount(g), "\n")
  cat("Maximum degree: ", max(degree(g)), "\n")
  cat("Minimum degree: ", min(degree(g)), "\n")
  cat("Number of edges:", ecount(g), "\n")
  cat("Radius:         ", radius(g), "\n")
  cat("Diameter:       ", diameter(g), "\n")
  cat("Density:        ", graph.density(g), "\n")
  cat("Transitivity:   ", transitivity(g), "\n")
  cat("Avg path length:", average.path.length(g), "\n")
  cat("Assortativity:  ", assortativity.degree(g), "\n")
  
  return(TRUE)
}

# ----------------------------------------------------------------------------
# 11/19/13 - random network simulation 
# from Brett's Matlab function of the same name
randomNetworkSim <- function(n=100, p=0.1, doFitPlot=F, doNetworkPlot=F, doHistPlot=F, 
                               useIgraph=F, numBins=10, filePrefix="random_network_sim") {
  # Erdos-Renyi
  # usage: A <- random_network_sim(100,.1,1)
  if(useIgraph) {
    g <- erdos.renyi.game(n, p)
    printIGraphStats(g)
    A <- get.adjacency(g)
  } else {
    A <- matrix(ncol=n, nrow=n, data=runif(n * n))
    # undirected no self-connections, no weights
    A <- triu(A < p, 1) # take upper triangle (k=1) of (random nxn < p)
    #A = sparse(A)
    # now make A symmetric, undirected
    A = A + t(A)
  }
  degrees <- rowSums(A)
  
  # calculates the scale free parameter and shows network
  if(useIgraph) {
    deg_counts <- floor(degree.distribution(g)*100)
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
  powerLawFit = displ$new(nz_deg_counts)
  xminEst = estimate_xmin(powerLawFit)
  cat("Erdos-Renyi Graph Estimate of x_min:", xminEst$xmin, "\n")
  cat("Erdos-Renyi Graph Estimate of x_min KS value:", xminEst$KS, "\n")
  cat("Erdos-Renyi Graph Estimate of scaling parameter alpha:", xminEst$pars, "\n")
  powerLawFit$setXmin(xminEst)
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
    g <- graph.adjacency(A)
    V(g)$size <- scale_a_b(degrees, 10, 20)
    png(paste(filePrefix, "_er_network.png", sep=""), width=1024, height=768)
    plot(g, layout=layout.fruchterman.reingold, edge.arrow.mode=0)
    dev.off()
  }
  
  return(A)
}

# ----------------------------------------------------------------------------
# 11/13/13 - scale free network simulation 
# from Brett's Matlab function of the same name
scaleFreeSim <- function(n=100, doFitPlot=F, doNetworkPlot=F, useIgraph=F, 
                           numBins=10, filePrefix="scale_free_sim") {
  
  # Baraba'si-Albert (BA) model generation of scale-free network
  if(useIgraph) {
    g <- barabasi.game(n, directed=F)
    foo <- printIGraphStats(g)
    A <- get.adjacency(g)
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
    deg_counts <- floor(degree.distribution(g)*100)
  } else {
    k_rows <- rowSums(A)
    bins <- numBins
    histObj <- hist(k_rows, bins, plot=F)
    deg_counts <- histObj$counts
  }
  #print(deg_counts)
  # find non zero entries for log plot  
  nz_idx <- which(deg_counts != 0, arr.ind=T)
  nz_deg_counts <- deg_counts[nz_idx]
  
  # Igraph's discrete power law fit
  powerLawFit = displ$new(nz_deg_counts)
  xminEst = estimate_xmin(powerLawFit)
  cat("Baraba'si-Albert Graph Estimate of x_min:", xminEst$xmin, "\n")
  cat("Baraba'si-Albert Graph Estimate of x_min KS value:", xminEst$KS, "\n")
  cat("Baraba'si-Albert Graph Estimate of scaling parameter alpha:", xminEst$pars, "\n")
  powerLawFit$setXmin(xminEst)
  if(doFitPlot) {
    png(paste(filePrefix, "_ba_powerlaw_fit.png", sep=""), width=1024, height=768)
    # data
    plot(powerLawFit, xlab="Degree", ylab="CDF", 
         main="Baraba'si-Albert Graph Degree Distribution with Power Law Fit")
    # fit
    lines(powerLawFit, col=2)
    dev.off()
  }
  
  # viz with lower triangle removed
  #A[upper.tri(A)] <- 0
  if(doNetworkPlot) {
    g <- graph.adjacency(A)
    V(g)$size <- scale_a_b(degrees, 10, 20)
    png(paste(filePrefix, "_ba_network.png", sep=""), width=1024, height=768)
    plot(g, layout=layout.fruchterman.reingold, edge.arrow.mode=0)
    dev.off()
  }
  
  return(A)
}

