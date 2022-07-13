#' @title LD plot
#' @description users can set threshold to plot the figures
#' @name clin_ld
#'
#' @param gwas dataset for making the figure after analysis
#' @import BEDMatrix
#' @import gaston
#' @return a Q-Q Plot
#' @export
#' @examples
#' \dontrun{
#' gwas<-read.csv(system.file("extdata", "all_assoc_logistic.csv", package = "CLINSKAT"))
#' clin_ld(gwas, p_val=0.001, geno=Brs_sample.bed)
#' }
#'

clin_ld <- function(gwas, p_val=0.001, geno=Brs_sample.bed){
  gwas <- gwas[which(gwas$P < p_val),]
  gwas$SNP_list<-paste0(gwas$SNP,'_',gwas$A1)
  SNP_list <- gwas$SNP_list

  g <- BEDMatrix::BEDMatrix('Brs_sample.bed')
  x_gwas <- g[,SNP_list]
  x <- gaston::as.bed.matrix(x_gwas)

  ld.x <- gaston::LD(x, c(1,ncol(x)))
  jpeg('Brs_1e4_LD.jpg',width=1800,height=1200)
  LD.plot(ld.x)
  dev.off()
}


