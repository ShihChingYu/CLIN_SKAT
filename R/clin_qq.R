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
  chisq = qchisq(gwas$P,1,lower.tail=FALSE)
  lambda<-median(chisq,na.rm=TRUE)/qchisq(0.5,1)

  qqman::qq(gwas, main=paste("lambda = " , lambda))
}


