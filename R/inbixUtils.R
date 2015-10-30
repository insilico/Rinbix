# -----------------------------------------------------------------------------
# inbixUtils.R - Bill White - 10/10/15
#
# Rinbix package utility/miscellaneous functions.

# ----------------------------------------------------------------------------
#' Fisher correlation transformation function R to Z distribution.
#' 
#' \code{fisherRtoZ}
#' 
#' @param x Correlation value -1 to 1.
#' @return Transformed correlation value.
#' @examples
#' fisherRtoZ(seq(from=-1, to=1, by=0.25))
#' @export
fisherRtoZ <- function(x) { 
  atanh(x)
}

# -----------------------------------------------------------------------------
#' Remove genes with low absolute values.
#' 
#' \code{geneLowValueFilter} removes genes with values in the lowest percentile.
#' 
#' @param dataMatrix Data frame with genes in rows and samples in columns.
#' @param percentile Numeric percentile threshold below which genes will be removed.
#' @return List with the mask used and filtered data frame.
#' @examples
#' data(testdata100ME4)
#' lowValFiltered <- geneLowValueFilter(testdata100ME4[, -ncol(testdata100ME4)])
#' @export
geneLowValueFilter <- function(dataMatrix, percentile=0.1) {
  # Remove gene profiles with low absolute values in dataMatrix. Returns:
  # 1) a logical vector mask identifying gene expression profiles in dataMatrix
  #    that have absolute expression levels in the lowest 10% of the data set.
  # 2) a data matrix containing filtered expression profiles.
  threshold <- quantile(as.matrix(dataMatrix), c(percentile))
  mask <- apply(dataMatrix, 1, function(x) all(x < threshold))
  fdata <- dataMatrix[!mask, ]
  
  # return the row mask and filtered data
  list(mask=!mask, fdata=fdata)
}

# -----------------------------------------------------------------------------
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
#' @param dataMatrix Data frame with genes in rows and samples in columns.
#' @param percentile Numeric variance percentile threshold below which genes will be removed.
#' @return List with the mask used and filtered data frame
#' @examples
#' data(testdata100ME4)
#' lowVarFiltered <- geneLowVarianceFilter(testdata100ME4[, -ncol(testdata100ME4)])
#' @export
geneLowVarianceFilter <- function(dataMatrix, percentile=0.1) {
  variances <- apply(as.matrix(dataMatrix), 1, var)
  threshold <- quantile(variances, c(percentile))
  mask <- apply(dataMatrix, 1, function(x) var(x) > threshold)
  fdata <- dataMatrix[mask, ]
  
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

# -----------------------------------------------------------------------------
#' Logarithmic spiral coordinate generator for igraph network node layout.
#'
#' \code{logSpiral}
#' 
#' Ported from Javascript found here: 
#' \url{http://www.pixelwit.com/blog/2008/05/how-to-draw-logarithmic-spiral/}
#'
#' @param centerX Numeric X origin of the spiral.
#' @param centerY Numeric Y origin of the spiral.
#' @param radius Numeric distance from origin to outer arm.
#' @param sides Numeric points or sides along the spiral's arm.
#' @param coils Numeric coils or full rotations. (Positive numbers spin clockwise, negative numbers spin counter-clockwise)
#' @param rotation Numeric overall rotation of the spiral. ('0'=no rotation, '1'=360 degrees, '180/360'=180 degrees)
#' @return matrix of x-y coordinates 'sides' rows and two columns.
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
  for(i in 1:sides) {
    # How far away from center
    away <- radius^(i / sides)
    # How far around the center.
    around <- i * aroundRadians + rotation;
    # Convert 'around' and 'away' to X and Y.
    x <- centerX + cos(around) * away;
    y <- centerY + sin(around) * away;
    #cat(x, "\t", y, "\n", sep="")
    results <- rbind(results, data.frame(x=x, y=y))
  }
  as.matrix(results)
}

# ----------------------------------------------------------------------------
#' Scales a numeric vector to a range (a, b).
#' 
#' \code{scaleAB}
#' 
#' @param v Numeric vector.
#' @param a Numeric minimum value in range.
#' @param b Numeric maximum value in range.
#' @return v scaled to range (a, b).
#' @examples
#' scaleAB(1:10, 0, 1)
#' @export
scaleAB <- function(v, a, b) {
  v <- v - min(v)
  v <- v / max(v)
  v <- v * (b - a)
  v + a
}

# ----------------------------------------------------------------------------
#' Compute the power series of a matrix.
#' 
#' A <- A^1 + A^2 + . . . + A^n
#' 
#' \code{sumOfPowers} 
#' 
#' @param A Matrix adjacency.
#' @param n Numeric maximum power of series adjacencyMatrix^n.
#' @param verbose Flag to send messages to stdout.
#' @return A^n.
#' @examples
#' A <- matrix(1:9, 3)
#' sumOfPowers(A, 3)
#' @export
sumOfPowers <- function(A, n, verbose=FALSE) {
  # g = A + A^2 + A^3 + ... A^n
  # Creates weighted graph g that includes direct connections (A),
  # all paths with 1 intermediate node (A^2)     -- first order indirect
  # all paths with 2 intermediate nodes (A^3)    -- second order indirect
  # ...                                          ...
  # up to n-1 intermediate connections (A^(n-1)) -- n-1 order indirect
  createStr <- "Creating g <- A"
  if(verbose) cat(createStr, "\n")
  g <- A
  currPowerOfA <- A
  if(n >= 2) { # this should be a pretty efficient way to create g
    for(i in 2:n){
      if(verbose) cat(paste(createStr," + A^", i, "\n", sep=""))
      currPowerOfA <- currPowerOfA %*% A
      g <- g + currPowerOfA
    }
  }
  g
}
