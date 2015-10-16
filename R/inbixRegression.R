# -----------------------------------------------------------------------------
# inbixRegression.R - Bill White - 10/10/15
#
# Rinbix package regression functions.

# -----------------------------------------------------------------------------
#' Get the main effect coefficent of a variable using generalized linear regression - glm
#' 
#' \code{getMainEffect} 
#' 
#' @param data Data frame with genes in columns and samples in rows.
#' @param geneName String column name of the gene used in the regression formula.
#' @param depVarName String phenotype column name used in the regression formula.
#' @param regressionFamily Stringregression family name for the glm.
#' @return regression beta coefficient.
#' @export
getMainEffect <- function(data, geneName, depVarName, regressionFamily) {
  regressionFormula <- as.formula(paste(depVarName, "~", paste("`", geneName, "`", sep=""), sep=" "))
  mainEffectModel <- glm(regressionFormula, family=regressionFamily, data=data)
  
  as.numeric(mainEffectModel$coef[2])
}

# -----------------------------------------------------------------------------
#' Run main effect and interaction models on two variables versus a phenotype.
#' 
#' \code{runGlms} 
#' 
#' @param geneA Vector of gene expression levels.
#' @param geneB Vector of gene expression levels.
#' @param pheno Vector of casee-control phenotypes 0/1.
#' @return Data frame of main effect and interaction coefficients.
#' @export
runGlms <- function(geneA, geneB, pheno) {
  # linear fits: binomial logit link function by default
  phenoGeneAFit <- glm(pheno ~ geneA, family=binomial)
  maineffectA_beta <- summary(phenoGeneAFit)$coefficients[2, "Estimate"]
  maineffectA_pval <- summary(phenoGeneAFit)$coefficients[2, "Pr(>|z|)"]
  
  phenoGeneBFit <- glm(pheno ~ geneB, family=binomial)
  maineffectB_beta <- summary(phenoGeneBFit)$coefficients[2, "Estimate"]
  maineffectB_pval <- summary(phenoGeneBFit)$coefficients[2, "Pr(>|z|)"]
  
  fullmodel <- glm(pheno ~ geneA*geneB, family=binomial)
  fullmodel_betaA <- summary(fullmodel)$coefficients[2, "Estimate"]
  fullmodel_pvalA <- summary(fullmodel)$coefficients[2, "Pr(>|z|)"]
  
  fullmodel_betaB <- summary(fullmodel)$coefficients[3, "Estimate"]
  fullmodel_pvalB <- summary(fullmodel)$coefficients[3, "Pr(>|z|)"]
  
  fullmodel_betaAB <- summary(fullmodel)$coefficients[4, "Estimate"]
  fullmodel_pvalAB <- summary(fullmodel)$coefficients[4, "Pr(>|z|)"]
  
  puremodel <- glm(pheno ~ geneA:geneB, family=binomial)
  puremodel_beta <- summary(puremodel)$coefficients[2, "Estimate"]
  puremodel_pval <- summary(puremodel)$coefficients[2, "Pr(>|z|)"]
  
  resultsP <- cbind(maineffectA_pval, maineffectB_pval, 
                    fullmodel_pvalA, fullmodel_pvalB, fullmodel_pvalAB, 
                    puremodel_pval)
  resultsB <- cbind(maineffectA_beta, maineffectB_beta, 
                    fullmodel_betaA, fullmodel_betaB, fullmodel_betaAB, 
                    puremodel_beta)
  result <- rbind(resultsB, resultsP)
  colnames(result) <- c("geneA_main", "geneB_main", "geneA_full", 
                        "geneB_full", "A*B_full", "A*B_pure")
  rownames(result) <- c("Betas", "p-values")
  
  return(result)
}
