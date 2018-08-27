## ------------------------------------------------------------------------
library(Rinbix)
data("testdata100ME4")

## ------------------------------------------------------------------------
rankedVarsRelief <- rankRelieff(testdata100ME4)
knitr::kable(head(rankedVarsRelief), row.names=F)

## ------------------------------------------------------------------------
rankedVarsIterRelief <- rankIterativeRelieff(testdata100ME4, 
                                             percentRemovePerIteration=10,
                                             targetNumAttributes=10, 
                                             verbose=FALSE)
knitr::kable(head(rankedVarsIterRelief$all.scores), row.names=F)

