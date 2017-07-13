# ----------------------------------------------------------------------------
# inbixSNP.R - Bill White - 12/4/15
#
# Rinbix package functions for SNPs.

#' Read a PLINK binary file of SNP genotypes and meta information.
#' 
#' \code{readPlinkBinary} 
#' 
#' @family SNP functions
#' @param plinkBasename \code{character} PLINK base filename
#' @return \code{list} of genotypes matrix and meta information data frames.
#' @examples
#' \dontrun{ 
#' plinkObj <- readPlinkBinary("wgas1") 
#' }
#' @export
readPlinkBinary <- function(plinkBasename) {
  snpStats::read.plink(bed = paste(plinkBasename, ".bed", sep = ""))
}
