# ----------------------------------------------------------------------------
# inbixViz.R - Bill White - 10/31/15
#
# Rinbix package for visualization.

# ----------------------------------------------------------------------------
#' Create an Igraph network from an adjacency matrix.
#' 
#' \code{adjacencyToNetList} 
#' 
#' @param Aadj Matrix correlation matrix.
#' @param thresholdType String threshold type: "hard" or "soft".
#' @param thresholdValue Numeric hard threshold correlation value or soft threshold power.
#' @param useAbs Flag take absolute value of the correlation matrix.
#' @param useWeighted Flag use weighted adjacency matrix versus binary.
#' @param verbose Flag send verbose messages to standard out .
#' @param plotFilename String image file for network plot.
#' @param groups Vector group assignments for network nodes.
#' @return list with igraph network, edge list and node list.
#' @examples
#' data("testdata10")
#' predictors <- testdata10[, -ncol(testdata10)]
#' Acorr <- cor(predictors)
#' netlist <- adjacencyToNetList(Acorr, 
#'                               thresholdType="hard", 
#'                               thresholdValue=0.2, 
#'                               useAbs=TRUE, 
#'                               useWeighted=TRUE, 
#'                               verbose=TRUE, 
#'                               groups=NULL)
#' @export
adjacencyToNetList <- function(Aadj=NULL, thresholdType="hard", thresholdValue=0.8, 
                         useAbs=TRUE, useWeighted=FALSE, verbose=FALSE,
                         groups=NULL) {
  if(is.null(Aadj)) {
    stop("adjacencyToNetList: Aadj is a required parameter")
  }
  if(verbose) cat("Begin adjacencyToNetList\n")
  if(thresholdType == "soft") {
    if(verbose) cat("Soft Threshold Aadj ^", thresholdValue, "\n")
    Aadj <- Aadj ^ thresholdValue
  } else {
    if(useAbs) {
      if(verbose) cat("Threshold abs(Aadj) >", thresholdValue, "\n")
      passThreshIdx <- abs(Aadj) > thresholdValue
    } else {
      if(verbose) cat("Threshold Aadj >", thresholdValue, "\n")
      passThreshIdx <- Aadj > thresholdValue
    }
    # create binary adjacency matrix from thresholding
    if(useWeighted) {
      if(verbose) cat("Making filtered weighted\n")
      Aadj[!passThreshIdx] <- 0
    } else {
      if(verbose) cat("Making filtered binary\n")
      Aadj[passThreshIdx] <- 1
      Aadj[!passThreshIdx] <- 0
    }
  }
  nodeLabels <- colnames(Aadj)
  numNodes <- length(nodeLabels)
  nodeGroups <- factor(seq(1, numNodes))
  if((!is.null(groups)) && (length(groups) != numNodes)) {
    nodeGroups <- factor(groups)
  }
  numGroups <- nlevels(nodeGroups)
  if(verbose) cat("Number of group levels", numGroups, "\n")

  # -----------------------------------------------------------------------------
  # create an igraph object from the edge weights matrix
  if(verbose) cat("Creating igraph object from weighted adjacency matrix\n")
  adjacencyNet <- igraph::graph_from_adjacency_matrix(Aadj, mode="upper", 
                                              weighted=TRUE, diag=FALSE)
  # degree is sum of connections
  nodeDegrees <- igraph::degree(adjacencyNet)
  # strength is weighted degree
  nodeStrengths <- igraph::strength(adjacencyNet)
  if(verbose) {
    print("Degrees")
    print(nodeDegrees)
    print("Weighted Degrees (Strengths)")
    print(nodeStrengths)
  }
  
  # node labels
  igraph::V(adjacencyNet)$label <- nodeLabels
  igraph::V(adjacencyNet)$shape <- "circle"

  # node colors from RColorBrewer
  nodeColors <- RColorBrewer::brewer.pal(numGroups, "Set3") 
  igraph::V(adjacencyNet)$color <- nodeColors[groups]
    
  igraph::V(adjacencyNet)$size <- nodeDegrees * 5
  igraph::V(adjacencyNet)$size2 <- V(adjacencyNet)$size

  # edge weights computed above in pairwise sums of connections
  igraph::E(adjacencyNet)$width <- igraph::E(adjacencyNet)$weight
  #E(adjacencyNet)$width <- (E(adjacencyNet)$weight / max(E(adjacencyNet)$weight)) * 7
  #E(adjacencyNet)$width <- 1

  # Offset vertex labels for smaller points (default=0).
  igraph::V(adjacencyNet)$label.dist <- ifelse(igraph::V(adjacencyNet)$size >= 5, 0, 0.2)

  # filterPercentile <- 0.95
  # filteredNetwork <- igraph::delete.edges(adjacencyNet, 
  #                                  igraph::E(adjacencyNet)[ weight < quantile(weight, filterPercentile) ])

  # -----------------------------------------------------------------------------
  # build graph data structures (node and edge lists) for D3 or other network viz
  if(verbose) cat("Creating node and edge lists\n")
  # edge list
  graphLinks <- NULL
  for(i in 1:numNodes) {
    nodeStart <- nodeLabels[i]
    for(j in 1:numNodes) {
      if(j <= i) { next }
      if(Aadj[i, j] > 0) {
        nodeEnd <- nodeLabels[j]
        edgeWeight <- Aadj[i, j]
        graphLinks <- rbind(graphLinks, data.frame(src=i, target=j, value=edgeWeight))
      }
    }
  }
  rownames(graphLinks) <- seq(1, nrow(graphLinks))
  uniqueNodes <- unique(c(unique(graphLinks$src), unique(graphLinks$target)))
  numUniqueNodes <- length(uniqueNodes)  
  # node list
  graphNodes <- NULL
  for(uniqIdx in 1:numUniqueNodes) {
    graphNodeIdx <- uniqueNodes[uniqIdx]
    graphNodes <- rbind(graphNodes, data.frame(idn=graphNodeIdx,
                                               name=nodeLabels[graphNodeIdx], 
                                               group=nodeGroups[graphNodeIdx], 
                                               size=nodeDegrees[graphNodeIdx],
                                               stringsAsFactors=F))
  }
  rownames(graphNodes) <- uniqueNodes
  
  # return data frames suitable for D3 or other viz
  edgeList <- data.frame(from=as.numeric(factor(graphLinks$src)) - 1, 
                         to=as.numeric(factor(graphLinks$target)) - 1 )
  nodeList <- cbind(idn=factor(graphNodes$idn, 
                    levels=graphNodes$idn), graphNodes)

  # return the Igraph object, edge list and node list
  if(verbose) cat("End adjacencyToNetList\n")
  list(net=adjacencyNet, links=edgeList, nodes=nodeList, groups=groups)
}
