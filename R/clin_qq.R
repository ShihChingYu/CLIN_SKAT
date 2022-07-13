#' @title Q-Q plot
#' @description users can set threshold to plot the figures
#' @name clin_qq
#'
#' @param gwas dataset for making the figure after analysis
#' @import qqman
#' @return a Q-Q Plot
#' @export
#' @examples
#' \dontrun{
#' gwas<-read.csv(system.file("extdata", "all_assoc_logistic.csv", package = "CLINSKAT"))
#' gg(gwas$P)
#' }
#'

clin_qq <- function(gwas){
  qqman::qq(gwas)
}


