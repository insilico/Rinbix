\name{runMainEffectsTest}
\alias{runMainEffectsTest}
\title{Get the main effect of a variable using generalized linear regression - glm}
\usage{
runMainEffectsTest(data, variableName, depVarName, regressionFamily,
  numCovariates)

runMainEffectsTest(data, variableName, depVarName, regressionFamily,
  numCovariates)
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

  \item{data}{data frame with genes in columns and samples
  in rows}

  \item{variableName}{name of the variable to consider}

  \item{depVarName}{name of the phenotype variable}

  \item{regressionFamily}{glm regression family name}

  \item{numCovariates}{number of included covariates}
}
\value{
the regression beta coefficient

data frame with gene, convergence status, beta coefficient,
p-value, standard error and standardized beta columns
}
\description{
\code{runMainEffectsTest}

\code{runMainEffectsTest}
}

