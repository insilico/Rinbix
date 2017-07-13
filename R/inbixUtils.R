# inbixUtils.R - Bill White - 10/10/15
#
# Rinbix package utility/miscellaneous functions.

#' Fisher correlation transformation function R to Z distribution.
#' 
#' \code{fisherRtoZ}
#' 
#' @keywords math
#' @family utility functions
#' @param x \code{numeric} Correlation value -1 to 1.
#' @return \code{numeric} Transformed correlation value.
#' @examples
#' fisherRtoZ(seq(from = -1, to = 1, by = 0.25))
#' @export
fisherRtoZ <- function(x) { 
  # catch +/-Inf in transform
  x[x > 0.9999999999] <- 0.9999999999
  x[x < -0.9999999999] <- -0.9999999999
  atanh(x)
}

#' Remove genes with low coefficient of variation.
#' 
#' \code{geneLowCoefOfVarFilter} removes genes with coefficient of variation less than value.
#' 
#' @family utility functions
#' @family filter functions
#' @family microarray functions
#' @param dataMatrix \code{data.frame} with genes in rows and samples in columns.
#' @param coefOfVarThreshold \code{numeric} threshold below which genes will be removed.
#' @return \code{list} with the mask used and filtered data frame.
#' @examples
#' data(testdata100ME4)
#' lowCVFiltered <- geneLowCoefOfVarFilter(t(testdata100ME4[, -ncol(testdata100ME4)]))
#' @export
geneLowCoefOfVarFilter <- function(dataMatrix, coefOfVarThreshold = 0.1) {
  coefOfVars <- apply(dataMatrix, 1, function(x) { sd(x) / abs(mean(x)) })
  # the smaller the threshold, the higher the experimental effect relative 
  # to the measurement precision
  # filter the data matrix
  fdata <- dataMatrix[coefOfVars < coefOfVarThreshold, ]
  # return the filtered data
  list(fdata = fdata)
}

#' Remove genes with low absolute values.
#' 
#' \code{geneLowValueFilter} removes genes with values in the lowest percentile.
#' 
#' @family utility functions
#' @family filter functions
#' @family microarray functions
#' @family inbix synonym functions
#' @param dataMatrix \code{data.frame} with genes in rows and samples in columns.
#' @param percentile \code{numeric} percentile threshold below which genes will be removed.
#' @return \code{list} with the mask used and filtered data frame.
#' @examples
#' data(testdata100ME4)
#' lowValFiltered <- geneLowValueFilter(t(testdata100ME4[, -ncol(testdata100ME4)]))
#' @export
geneLowValueFilter <- function(dataMatrix, percentile = 0.1) {
  # Remove gene profiles with low absolute values in dataMatrix. Returns:
  # 1) a logical vector mask identifying gene expression profiles in dataMatrix
  #    that have absolute expression levels in the lowest 10% of the data set.
  # 2) a data matrix containing filtered expression profiles.
  threshold <- quantile(as.matrix(dataMatrix), c(percentile))
  mask <- apply(dataMatrix, 1, function(x) all(x < threshold))
  fdata <- dataMatrix[!mask, ]
  
  # return the row mask and filtered data
  list(mask = !mask, fdata = fdata)
}

#' Remove genes with low expression variance.
#' 
#' \code{geneLowVarianceFilter} removes genes with variance below a 
#' variance percentile threshold.
#' 
#' MATLAB: ... calculates the variance for each gene expression profile, 
#' which identifies the gene expression profiles with a variance less than the 10th 
#' percentile. Mask is a logical vector with one element for each row in Data. 
#' The elements of Mask corresponding to rows with a variance greater than the threshold 
#' have a value of 1, and those with a variance less than the threshold are 0.
#'
#' @family utility functions
#' @family filter functions
#' @family microarray functions
#' @family inbix synonym functions
#' @param dataMatrix \code{data.frame} with genes in rows and samples in columns.
#' @param percentile \code{numeric} variance percentile threshold below which genes will be removed.
#' @return \code{list} with the mask used and filtered data frame
#' @examples
#' data(testdata100ME4)
#' lowVarFiltered <- geneLowVarianceFilter(t(testdata100ME4[, -ncol(testdata100ME4)]))
#' @export
geneLowVarianceFilter <- function(dataMatrix, percentile = 0.1) {
  variances <- apply(as.matrix(dataMatrix), 1, var)
  threshold <- quantile(variances, c(percentile))
  mask <- apply(dataMatrix, 1, function(x) var(x) > threshold)
  fdata <- dataMatrix[mask, ]
  
  # return the row mask and filtered data
  list(mask = mask, fdata = fdata)
}

#' Logarithmic spiral coordinate generator for igraph network node layout.
#'
#' \code{logSpiral}
#' 
#' Ported from Javascript found here: 
#' \url{http://www.pixelwit.com/blog/2008/05/how-to-draw-logarithmic-spiral/}
#'
#' @keywords math
#' @family utility functions
#' @param centerX \code{numeric} X origin of the spiral.
#' @param centerY \code{numeric} Y origin of the spiral.
#' @param radius \code{numeric} distance from origin to outer arm.
#' @param sides \code{numeric} points or sides along the spiral's arm.
#' @param coils \code{numeric} coils or full rotations. (Positive numbers spin clockwise, negative numbers spin counter-clockwise)
#' @param rotation \code{numeric} overall rotation of the spiral. ('0'=no rotation, '1'=360 degrees, '180/360'=180 degrees)
#' @return \code{matrix} of x-y coordinates 'sides' rows and two columns.
#' @export
#' @keywords internal
logSpiral <- function(centerX, centerY, radius, sides, coils, rotation) {
  # Start at the center.
  #moveTo(centerX, centerY);
  # How far to rotate around center for each side.
  aroundStep <- coils / sides # 0 to 1 based.
  # Convert aroundStep to radians.
  aroundRadians <- aroundStep * 2 * pi
  # Convert rotation to radians.
  rotation <- rotation * 2 * pi
  results <- NULL
  # For every side, step around and away from center.
  for (i in 1:sides) {
    # How far away from center
    away <- radius^(i / sides)
    # How far around the center.
    around <- i * aroundRadians + rotation;
    # Convert 'around' and 'away' to X and Y.
    x <- centerX + cos(around) * away;
    y <- centerY + sin(around) * away;
    #cat(x, "\t", y, "\n", sep = "")
    results <- rbind(results, data.frame(x = x, y = y))
  }
  as.matrix(results)
}

#' Look up gene names for a list of Ensembl IDs.
#' 
#' \code{lookupGenesFromEnsemblIds}
#' 
#' @family utility functions
#' @param ensemblIds \code{character} vector of Ensembl IDs.
#' @return \code{character} vector of gene names/symbols.
#' @examples
#' ensemblIds <- c("ENSG00000000003", "ENSG00000000005")
#' geneNames <- lookupGenesFromEnsemblIds(ensemblIds)
#' @export
lookupGenesFromEnsemblIds <- function(ensemblIds) {
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  geneNames <- biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), 
                              filters = 'ensembl_gene_id', 
                              values = ensemblIds, 
                              mart = ensembl)
  geneNames
}

#' Scales a numeric vector to a range (a, b).
#' 
#' \code{scaleAB}
#' 
#' @keywords math array
#' @family utility functions
#' @param v \code{numeric} vector.
#' @param a \code{numeric} minimum value in range.
#' @param b \code{numeric} maximum value in range.
#' @return \code{numeric} v scaled to range (a, b).
#' @examples
#' scaleAB(1:10, 0, 1)
#' @export
scaleAB <- function(v, a, b) {
  v <- v - min(v)
  v <- v / max(v)
  v <- v * (b - a)
  v + a
}

#' Compute the power series of a matrix.
#' 
#' A <- A^1 + A^2 + . . . + A^n
#' 
#' \code{sumOfPowers} 
#' 
#' @keywords math array
#' @family utility functions
#' @param A \code{matrix} adjacency.
#' @param n \code{numeric} maximum power of series adjacencyMatrix^n.
#' @param verbose \code{logical} to send messages to stdout.
#' @return \code{matrix} A^n.
#' @examples
#' A <- matrix(1:9, 3)
#' sumOfPowers(A, 3)
#' @export
sumOfPowers <- function(A, n, verbose = FALSE) {
  # g = A + A^2 + A^3 + ... A^n
  # Creates weighted graph g that includes direct connections (A),
  # all paths with 1 intermediate node (A^2)     -- first order indirect
  # all paths with 2 intermediate nodes (A^3)    -- second order indirect
  # ...                                          ...
  # up to n-1 intermediate connections (A^(n-1)) -- n-1 order indirect
  createStr <- "Creating g <- A"
  if (verbose) cat(createStr, "\n")
  g <- A
  currPowerOfA <- A
  if (n >= 2) { # this should be a pretty efficient way to create g
    for (i in 2:n) {
      if (verbose) cat(paste(createStr," + A^", i, "\n", sep = ""))
      currPowerOfA <- currPowerOfA %*% A
      g <- g + currPowerOfA
    }
  }
  g
}
