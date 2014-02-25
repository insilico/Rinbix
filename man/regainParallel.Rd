\name{regainParallel}
\alias{regainParallel}
\title{Parallel execution of the regression genetic association network algorithm}
\usage{
regainParallel(regressionData)
}
\arguments{
  \item{regressionData}{a data frame with gene in columns
  and samples in rows; the last column should be labeled
  'Class' and be 0 or 1 values}
}
\value{
regainMatrix a matrix of gene by gene regression
coefficients
}
\description{
\code{regainParallel}
}

