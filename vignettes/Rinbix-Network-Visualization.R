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
netLayout <- layout_in_circle(filteredNetwork, order=V(filteredNetwork))
plot(filteredNetwork, layout=netLayout, edge.curved=TRUE,
     vertex.label.color="blue", vertex.label.cex=1.0)

# D3 layout
simpleNetwork(netlist$links, Source="from", Target="to", 
              linkColour="black", 
              charge=-2000)
forceNetwork(Links=netlist$links, 
             Nodes=netlist$nodes, 
             Source="from", Target="to",
             NodeID="name", Group="group", 
             Nodesize="size", 
             linkColour="#afafaf", 
             fontSize=10, zoom=T, legend=T,
             opacity=0.9, charge=-1000, 
             width=600, height=600)

