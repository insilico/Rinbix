# ----------------------------------------------------------------------------
# inbixGSEA.R - Bill White - 10/31/15
#
# Rinbix package for gene set enrichment analysis (GSEA).

# ----------------------------------------------------------------------------
#' Get the Reactome pathway descriptions for a vector of gene symbols.
#' 
#' \code{getReactomePathways} 
#' 
#' @param geneSymbols \code{vector} gene symbols.
#' @return \code{\link[DOSE]{enrichResult-class}} enriched pathways with FDR control.
#' @family GSEA functions
#' @examples
#' data(geneListSymbols)
#' reactomePathDesc <- getReactomePathways(geneListSymbols)
#' @export
getReactomePathways <- function(geneSymbols) {
  geneIdsEntrez <- AnnotationDbi::mget(geneSymbols, 
                                       AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  pathwayEnrich <- ReactomePA::enrichPathway(gene=geneIdsEntrez, 
                                             pvalueCutoff=1, 
                                             qvalueCutoff=0.9, 
                                             minGSSize=2)
  pathwayEnrich
}

# ----------------------------------------------------------------------------
#' Get the KEGG enrichment analysis for a vector of gene symbols.
#' 
#' \code{getKEGGAnalysis} 
#' 
#' @param geneSymbols \code{vector} gene symbols.
#' @return \code{\link[DOSE]{enrichResult-class}} enriched pathways with FDR control.
#' @family GSEA functions
#' @examples
#' data(geneListSymbols)
#' keggEnrichment <- getKEGGAnalysis(geneListSymbols)
#' @export
getKEGGAnalysis <- function(geneSymbols) {
  geneIdsEntrez <- AnnotationDbi::mget(geneSymbols, 
  	                                   AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  geneIdsEntrez[is.na(geneIdsEntrez)] = " "
  keggEnrichment <- clusterProfiler::enrichKEGG(geneIdsEntrez, "human", 
                                                pvalueCutoff=0.05,
                                                use_internal_data=TRUE)
  keggEnrichment
}

# ----------------------------------------------------------------------------
#' Get the GO enrichment analysis for a vector of gene symbols.
#' 
#' \code{getGOAnalysis} 
#' 
#' @param geneSymbols \code{vector} gene symbols.
#' @return \code{\link[DOSE]{enrichResult-class}} enriched pathways with FDR control.
#' @family GSEA functions
#' @examples
#' data(geneListSymbols)
#' goEnrichment <- getGOAnalysis(geneListSymbols)
#' @export
getGOAnalysis <- function(geneSymbols) {
  geneIdsEntrez <- AnnotationDbi::mget(geneSymbols, 
  	                                   AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  goEnrichment <- clusterProfiler::enrichGO(geneIdsEntrez, "human", "MF")
  goEnrichment
}
