#' @title Manhattan
#' @description users can set threshold to plot the figures
#' @name clin_manhattan
#'
#' @param gwas dataset for making the figure after analysis
#' @param sug_val the value to draw suggestive line
#' @param geno_val the value to draw genome-wide significant line
#' @import qqman
#' @return a Manhattan Plot
#' @export
#' @examples
#' \dontrun{
#' gwas<-read.csv(system.file("extdata", "all_assoc_logistic.csv", package = "CLINSKAT"))
#' clin_manhattan(gwas)
#' }
#'

clin_manhattan <- function(gwas, sug_val=3, geno_val=5){
  qqman::manhattan(gwas, main = "Manhattan Plot", annotateTop = TRUE, ylim = c(0, 9), cex = 0.8, cex.axis = 0.9,
            col = c("blue4", "orange3"), suggestiveline = sug_val, genomewideline = geno_val, chrlabs = c(1:22,'P','Q'))
}


