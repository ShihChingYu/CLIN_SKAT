#' @title relate significant variants to pathways in MSigDB database
#' @description relate a list of genes, convert ensemble ID to gene ID and pathway enrichment analysis
#' @name relate2GeneDisease
#'
#' @param SNPdata significant variants of interest from statistical analysis
#' @param species model organisms such as hamun and mouse
#' @param category gene sets including "H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"
#' @import biomaRt
#' @import msigdbr
#' @import clusterProfiler
#' @return a table including pathways, gene ratio, background ratio, and p-value
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata", "asso001.csv", package = "skatvcf"),
#'              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' gene_disease<-relate2GeneDisease(SNPdata=dat, species="Homo sapiens", category="H")
#'

relate2GeneDisease<-function(SNPdata, species=c("Homo sapiens", "Mus musculus"),
                             category=c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")){
  regions <- paste(SNPdata$CHR, SNPdata$BP, SNPdata$BP, sep=":")
  mart <- biomaRt::useEnsembl("ensembl","hsapiens_gene_ensembl")
  results <- biomaRt::getBM(attributes = c("entrezgene_id", "chromosome_name",
                                           "start_position","end_position"),
                            filters = c("chromosomal_region"),
                            values=regions,
                            mart=mart)

  m_t2g <- msigdbr::msigdbr(species = species, category = category) %>%
    dplyr::select(gs_name, entrez_gene)

  gene<-results$entrezgene_id
  em <- clusterProfiler::enricher(gene, TERM2GENE=m_t2g)
  em_df<-as.data.frame(em)
  return(em_df)
}



