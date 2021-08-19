#' @title relate a list of genes to disease table based on Clinvar
#' @description relate a list of genes to disease table based on Clinvar
#' @name relatetoDisease
#'
#' @param genelist genes of interest corresponding to disease table
#' @import dplyr
#' @return a table including genes and corresponding diseases
#' @export
#' @examples
#' genes<-c("BCOR",  "NAA10")
#' gene_disease<-relatetoDisease(genelist=genes)
#'

relatetoDisease<-function(genelist){
  df<-data.frame(paste(genelist, collapse=", "))
  colnames(df)<-"gene"
  diseases<-read.csv(system.file("extdata", "merge_duDisease.csv", package = "skatvcf"),
                     stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
  g_d<-dplyr::left_join(df, diseases, by=c("gene" = "AssociatedGenes"))
  return(g_d)
}

