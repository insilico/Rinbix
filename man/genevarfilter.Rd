\name{genevarfilter}
\alias{genevarfilter}
\title{Remove gene profiles with low expression variance}
\usage{
genevarfilter(dataMatrix, percentile)
}
\arguments{
  \item{dataMatrix}{a data frame with genes in rows and
  samples in columns}

  \item{percentile}{variance threshold below which genes
  will be removed}
}
\value{
a list with the mask used and filtered data frame
}
\description{
\code{genelowvarfilter} removes genes with variance below a
variance percentile threshold
}

