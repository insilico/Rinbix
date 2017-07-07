# ----------------------------------------------------------------------------
# inbixGSEA.R - Bill White - 10/31/15
#
# Rinbix package for gene set enrichment analysis (GSEA).

# ----------------------------------------------------------------------------
#' Lookup a gene description from a gene symbol.
#' 
#' \code{lookupGeneDescBiomart} 
#' 
#' @family SNP functions
#' @param genes.list \code{character} vector of gene symbols
#' @return A list containing:
#' \describe{
#'   \item{run.results}{data frame of run results of each sim type}
#'   \item{Gene}{gene symbol}
#'   \item{Description}{gene description}
#' }
#' @examples
#' \dontrun{ 
#' gene.info <- lookupGeneDescBiomart(c("CEPN", "ATMIN"))
#' }
#' @export
lookupGeneDescBiomart <- function(genes.list=NULL) {
  if(is.null(genes.list)) {
    stop("No gene symbols in gene.list")
  }
  # lookup the ENSEMBL ID and get the gene symbol
  # load the gene info database
  ensembl2gene.biomart <- biomaRt::useEnsembl(biomart="ensembl", 
                                              dataset="hsapiens_gene_ensembl")
  gene.info <- biomaRt::getBM(attributes=c('ensembl_gene_id', 
                                           'hgnc_symbol', 
                                           'description'), 
                              filters='hgnc_symbol', 
                              values=genes.list,
                              uniqueRows=TRUE,
                              mart=ensembl2gene.biomart)
  colnames(gene.info) <- c("Gene.ID", "Gene.Symbol", "Description")
  # sort data frame returned by snp then gene ID
  gene.info[order(gene.info$Gene.Symbol), ]
}

# ----------------------------------------------------------------------------
#' Map a list of SNP rs# to ensembl genes with description.
#' 
#' \code{mapSNPsToGenesBiomart} 
#' 
#' @family SNP functions
#' @param snps.list \code{character} vector of rs# SNP IDs
#' @return A list containing:
#' \describe{
#'   \item{run.results}{data frame of run results of each sim type}
#'   \item{RefSNP}{RefSNP rs number}
#'   \item{Ensembl}{ensembl gene ID}
#'   \item{Gene}{gene symbol}
#'   \item{Description}{gene descriptions}
#' }
#' @examples
#' snp.info <- mapSNPsToGenesBiomart(c("rs1048194", "rs1553460"))
#' @export
mapSNPsToGenesBiomart <- function(snps.list=NULL) {
  if(is.null(snps.list)) {
    stop("No SNPs")
  }
  # load the SNP info database
  snp2ensembl.biomart <- biomaRt::useMart("ENSEMBL_MART_SNP", 
                                          dataset="hsapiens_snp")
  snp.id.info <- biomaRt::getBM(c("refsnp_id", "ensembl_gene_stable_id", 
                                  "ensembl_transcript_stable_id",
                                  "chr_name", "chrom_start", "chrom_end"),
                                filters="snp_filter",
                                values=snps.list,
                                uniqueRows=FALSE,
                                mart=snp2ensembl.biomart)
  map.ensembl.clean <- snp.id.info[snp.id.info$ensembl_gene_stable_id != "", ]
  # lookup the ENSEMBL ID and get the gene symbol
  # load the gene info database
  ensembl.id <- map.ensembl.clean$ensembl_transcript_stable_id
  ensembl2gene.biomart <- biomaRt::useEnsembl(biomart="ensembl", 
                                              dataset="hsapiens_gene_ensembl")
  snp.gene.info <- biomaRt::getBM(attributes=c('ensembl_gene_id', 
                                               'hgnc_symbol', 
                                               'description'), 
                                  filters='ensembl_transcript_id', 
                                  values=ensembl.id,
                                  uniqueRows=TRUE,
                                  mart=ensembl2gene.biomart)
  snp.genes.annot <- cbind(map.ensembl.clean[match(snp.gene.info$ensembl_gene_id, 
                                                   ensembl.id), 1], 
                           snp.gene.info)
  colnames(snp.genes.annot) <- c("RefSNP", "Ensembl", "Gene", "Description")
  # remove snps with no associated gene symbol
  snp.genes.annot <- snp.genes.annot[snp.genes.annot$Gene != "", ]
  # sort data frame returned by snp then gene ID
  snp.genes.annot[order(snp.genes.annot$RefSNP, snp.genes.annot$Ensembl), ]
}

# ----------------------------------------------------------------------------
#' Get the Entrez IDs for a vector of gene symbols.
#' 
#' \code{geneSymbolToEntrezID} 
#' 
#' @param gene.symbols \code{vector} gene symbols.
#' @return  \code{vector} of Entrez IDs.
#' @family GSEA functions
#' @examples
#' data(geneListSymbols)
#' symbolsEntrez <- geneSymbolToEntrezID(geneListSymbols)
#' @export
geneSymbolToEntrezID <- function(gene.symbols) {
  geneIdsEntrez <- AnnotationDbi::mget(gene.symbols, 
                                       AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL), 
                                       ifnotfound=NA)
  geneIdsEntrez
}

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
  geneIdsEntrez <- geneSymbolToEntrezID(geneSymbols)
  pathwayEnrich <- ReactomePA::enrichPathway(gene=geneIdsEntrez, 
                                             pvalueCutoff=1, 
                                             qvalueCutoff=0.9, 
                                             minGSSize=2,
                                             readable=TRUE)
  pathwayEnrich
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
  geneIdsEntrez <- geneSymbolToEntrezID(geneSymbols)
  goEnrichment <- clusterProfiler::enrichGO(gene=geneIdsEntrez, 
                                            OrgDb="org.Hs.eg.db", 
                                            ont="MF")
  goEnrichment
}
