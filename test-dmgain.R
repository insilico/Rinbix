#!/usr/bin/Rscript
#
# test-dmgain.R - Bill White - 7/29/14

rm(list=ls())
options(echo=F)

# insilico bioinformatics toolbox
library(Rinbix)

dsRows <- 100
colRange <- c(100)
for(i in 1:length(colRange)) {
  thisCols <- colRange[i]
  ds <- createRandomRegressionDataset(dsRows, thisCols)
  writeRegressionDataAsInbixNumeric(ds, "random")
  ptm <- proc.time()
  dcResults <- dcgain(ds)
  dmResults <- dmgain(ds)
  rgResults <- regainParallel(ds)
  cat("Rows:", dsRows, ", cols:", thisCols, ", elapsed time:", proc.time() - ptm, "\n")
}

print(dcResults$scores[1:5,1:5])
print(dmResults$scores[1:5,1:5])
print(rgResults[1:5,1:5])
