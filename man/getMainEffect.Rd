\name{getMainEffect}
\alias{getMainEffect}
\title{Get the main effect of a variable using generalized linear regression - glm}
\usage{
getMainEffect(data, geneName, depVarName, regressionFamily)
}
\arguments{
  \item{data}{is the data frame with genes in columns and
  samples in rows}

  \item{geneName}{is the column name of the gene used in
  the regression formula}

  \item{depVarName}{is the phenotype column name used in
  the regression formula}

  \item{regressionFamily}{is the regression family name for
  the glm}
}
\value{
the regression beta coefficient
}
\description{
\code{runMainEffectsTest}
}

