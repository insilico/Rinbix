# ----------------------------------------------------------------------------
# inbixPathway.R - Bill White - 10/31/15
#
# Rinbix package for pathway analysis.

# ----------------------------------------------------------------------------
#' Get the Reactome pathway descriptions for a vector of gene symbols.
#' 
#' \code{getReactomePathways} 
#' 
#' @param geneSymbols vector gene symbols.
#' @return enrichResult enriched pathways with FDR control.
#' @export
getReactomePathways <- function(geneSymbols) {
  geneIdsEntrez <- AnnotationDbi::mget(geneSymbols, 
                                       AnnotationDbi::revmap(org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  pathwayEnrich <- ReactomePA::enrichPathway(gene=geneIdsEntrez, 
                                             pvalueCutoff=1, 
                                             qvalueCutoff=0.9, 
                                             minGSSize=2)
}

# ----------------------------------------------------------------------------
#' Get the KEGG enrichment analysis for a vector of gene symbols.
#' 
#' \code{getKEGGAnalysis} 
#' 
#' @param geneSymbols vector gene symbols.
#' @return enrichResult enrichment KEGG categories with FDR control.
#' @export
getKEGGAnalysis <- function(geneSymbols) {
  geneIdsEntrez <- AnnotationDbi::mget(geneSymbols, AnnotationDbi::revmap(org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  geneIdsEntrez[is.na(geneIdsEntrez)] = " "
  keggEnrichment <- clusterProfiler::enrichKEGG(geneIdsEntrez, "human", 
                                                pvalueCutoff=0.05,
                                                use_internal_data=TRUE)
}

# ----------------------------------------------------------------------------
#' Get the GO enrichment analysis for a vector of gene symbols.
#' 
#' \code{getGOAnalysis} 
#' 
#' @param geneSymbols vector gene symbols.
#' @return enrichResult enrichment GO categories with FDR control.
#' @export
getGOAnalysis <- function(geneSymbols) {
  geneIdsEntrez <- AnnotationDbi::mget(geneSymbols, AnnotationDbi::revmap(org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  goEnrichment <- clusterProfiler::enrichGO(geneIdsEntrez, "human", "MF")
}
