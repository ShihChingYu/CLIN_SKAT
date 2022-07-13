#' @title barplot of significant genes after skat analysis
#' @description users can use the output from skat_asso() to draw the plot
#' @name skat_gene_bar
#'
#' @param asso_result dataset for making the figure after analysis
#' @param pval the threshold of p value
#' @param width bar width
#' @param space bar space
#' @return a bar Plot
#' @export
#' @examples
#' \dontrun{
#' gwas<-read.csv(system.file("extdata", "plink_clinical_82samples_620snp_skatD.csv", package = "CLINSKAT"))
#' skat_gene_bar(asso_result)
#' }
#'

skat_gene_bar <- function(asso_result, pval=0.0001, width=0.5, space=1){
  asso_result1<-asso_result[asso_result$P.value<pval, ]
  asso_result1$logP<-(-asso_result1$P.value)

  barplot(asso_result1$logP, names.arg = asso_result1$SetID, col=c("royalblue", "seagreen3"),
          width = width, space = space, ylab = "-log10(P)", main="CLIN_SKAT significant gene")
}



