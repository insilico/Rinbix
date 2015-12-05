## ---- echo=FALSE, message=FALSE------------------------------------------
library(Rinbix)
library(igraph)
library(networkD3)
library(RColorBrewer)

## ------------------------------------------------------------------------
# get a small data set and get the correlation between predictors
data("testdata10")
predictors <- testdata10[, -ncol(testdata10)]
Acorr <- cor(predictors)

# create a network from the correlation/adjacency matrix 
netlist <- adjacencyToNetList(Acorr, 
                              thresholdType="hard", 
                              thresholdValue=0.2, 
                              useAbs=TRUE, 
                              useWeighted=TRUE, 
                              verbose=TRUE)
printIGraphStats(netlist$net)

# Igraph circle layout
filteredNetwork <- netlist$net
netLayout <- igraph::layout_in_circle(filteredNetwork, order=igraph::V(filteredNetwork))
plot(filteredNetwork, layout=netLayout, edge.curved=TRUE,
     vertex.label.color="blue", vertex.label.cex=1.0)

# D3 layout - simple
networkD3::simpleNetwork(netlist$links, Source="Source", Target="Target", 
                         linkColour="black", 
                         charge=-2000)

# D3 layout - force 
networkD3::forceNetwork(Links=netlist$links, 
                        Nodes=netlist$nodes, 
                        Source="Source", Target="Target",
                        NodeID="NodeID", Group="Group", 
                        Nodesize="Size", 
                        linkColour="#afafaf", 
                        fontSize=10, zoom=T, legend=T,
                        opacity=0.9, charge=-1000, 
                        width=600, height=600)

