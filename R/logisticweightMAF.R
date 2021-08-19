#' @title weighting based on minor allele frequency
#' @description read vcf.gz and transform into Plink to do SAKT analysis
#' @name Get_Logistic_Weights_MAF
#'
#' @param dat input data with column chr, start, end, ref, alt to retrieve MAF in specified population
#' @param pop population frequency in TransAT package. Default is "db_gnomAD_exome_freq"
#' @param par1 default is 0.07
#' @param par2 default is 150
#' @import TransAT
#' @return a table of weighted value of MAF
#' @export
#' @examples
#' \dontrun{
#' dat<-read.csv(system.file("extdata", "anno_freq_data.csv", package = "skatvcf"),
#'              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' weights_MAF<-Get_Logistic_Weights_MAF(dat, pop="db_gnomAD_exome_freq", par1= 0.07, par2=150)
#' }
#'
Get_Logistic_Weights_MAF<-function(dat, pop="db_gnomAD_exome_freq", par1= 0.07, par2=150){

  dat1<-TransAT::pop_freq(dat, pop)
  match_col<-grep("freq$", colnames(dat1), ignore.case = T)
  dat2<-dat1[,match_col]

  variants_maf<-list()
  maf<-list()
  for (i in 1:nrow(dat2)){
    variants<-dat2[i, ]
    for (j in 1:ncol(dat2)){
      if(j %% 2 == 1){
        maf[j]<-min(variants[,j], variants[,j+1])
      }
      maf_table<-do.call("cbind", maf)
    }
    variants_maf[[i]]<-maf_table
  }

  var_maf_table<-do.call("rbind", variants_maf)
  colnames(var_maf_table)<-unlist(unique(strsplit(colnames(dat2), "_([^_]*_[^_]*)$")))

  n<-length(var_maf_table)
  weights<-rep(0,n)
  IDX<-which(var_maf_table > 0)
  if(length(IDX) == 0){
    stop("No polymorphic SNPs")
  } else {
    x1<-(var_maf_table[IDX] - par1) * par2
    weights[IDX]<-exp(-x1)/(1+exp(-x1))
  }

  weights_df<-data.frame(matrix(weights, nrow=3))
  colnames(weights_df)<-colnames(var_maf_table)
  return(weights_df)
}
