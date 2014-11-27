#!/usr/bin/Rscript
#
# test-dmgain.R - Bill White - 7/29/14

rm(list=ls())
options(echo=F)

# insilico bioinformatics toolbox
library(Rinbix)

dsRows <- 10
thisCols <- 10

ds <- createRandomRegressionDataset(dsRows, thisCols)
writeRegressionDataAsInbixNumeric(ds, "random")

ptm <- proc.time()
rgResults <- regainParallel(ds, stdBetas = TRUE, absBetas = TRUE)
cat("Rows:", dsRows, ", cols:", thisCols, ", elapsed time:", proc.time() - ptm, "\n")
print(rgResults[1:5,1:5])

ptm <- proc.time()
rgResultsCpp <- regainInbix(ds, stdBetas = TRUE, absBetas = TRUE)
cat("Rows:", dsRows, ", cols:", thisCols, ", elapsed time:", proc.time() - ptm, "\n")
print(rgResultsCpp$reGAIN[1:5,1:5])
