\name{runMainEffectsTest}
\alias{runMainEffectsTest}
\title{Get the main effect of a variable using generalized linear regression - glm}
\usage{
runMainEffectsTest(data, variableName, depVarName, regressionFamily,
  numCovariates)
}
\arguments{
  \item{data}{data frame with genes in columns and samples
  in rows}

  \item{variableName}{name of the variable to consider}

  \item{depVarName}{name of the phenotype variable}

  \item{regressionFamily}{glm regression family name}

  \item{numCovariates}{number of included covariates}
}
\value{
data frame with gene, convergence status, beta coefficient,
p-value, standard error and standardized beta columns
}
\description{
\code{runMainEffectsTest}
}

