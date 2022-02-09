#' @title Convert vcf data into plink format
#' @description read vcf.gz and transform into Plink to do SKAT analysis
#' @name convert_vcf
#'
#' @param vcf A vcf.gz format to read into the function
#' @param outputfile file name when write out
#' @param formula an object of class formula
#' @param out_type an indicator of the outcome type. "C" for the continuous outcome and "D" for
#' the dichotomous outcome
#' @param SNP_weight output object from Read_SNP_WeightFile (default=NULL)
#' @import vcfR
#' @import SKAT
#' @import data.table
#' @return a plink format with .map .bim .fam .bed .SetID for Generate_SSD_SetID
#' @export
#' @examples
#' \dontrun{
#' dat<-system.file("extdata", "sub10k.vcf.gz", package = "CLINSKAT")
#' formula=as.formula(y~1)
#' plink_data<-convert_vcf(vcf=dat, outputfile="subbackToPLINK", formula=formula, out_type="C")
#' }
#'

convert_vcf<-function(vcf, outputfile="subbackToPLINK", formula, out_type="C", SNP_weight=NULL){

  data_command<-paste0("C:/Users/user2/Desktop/survival_skat/plink --vcf ", vcf, " --allow-no-sex --make-bed --recode --out ", outputfile)
  system(data_command)

  #read backtoPLINK files
  dat_map<-paste0(getwd(), "/", outputfile, ".map")
  dat_bim<-paste0(getwd(), "/", outputfile, ".bim")
  dat_fam<-paste0(getwd(), "/", outputfile, ".fam")
  dat_bed<-paste0(getwd(), "/", outputfile, ".bed")

  bim<-data.table::fread(dat_bim, sep = "\t")
  bim[,2]<-data.frame(paste0("SNP", 1:nrow(bim)))
  data.table::fwrite(bim, dat_bim, col.names = F, row.names = F, sep = "\t")

  fam<-data.table::fread(dat_fam)
  fam$V6<-sample(1:2, 53, replace = T) #create artificial numbers in terms of gender
  fam$V7<-sample(30:80, 53, replace = T) #create artificial numbers in terms of age
  data.table::fwrite(fam[,2:7], dat_fam, col.names = F, row.names = F, sep = "\t")

  dat_vcf<-vcfR::read.vcfR(vcf, "hg19")
  vcf_format<-vcfR::extract_info_tidy(dat_vcf)
  vcf_format_gene<-data.frame(vcf_format$Gene.refGene)
  SetID<-cbind(vcf_format_gene$vcf_format.Gene.refGene, bim$V2)
  write.table(SetID, paste0(getwd(), "/", outputfile, ".SetID"), col.names = F, row.names = F, sep = "\t", quote = F)

  dat_SetID<-paste0(getwd(), "/", outputfile, ".SetID")
  dat_SSD<-paste0(getwd(), "/", outputfile, ".SSD")
  dat_Info<-paste0(getwd(), "/", outputfile, ".SSD.info")

  SKAT::Generate_SSD_SetID(dat_bed, dat_bim, dat_fam, dat_SetID, dat_SSD, dat_Info)

  SSD.INFO<-SKAT::Open_SSD(dat_SSD, dat_Info)
  FAM<-SKAT::Read_Plink_FAM(dat_fam, Is.binary=FALSE)
  y<-FAM$Phenotype
  obj<-SKAT::SKAT_Null_Model(formula, out_type=out_type)
  out<-SKAT::SKAT.SSD.All(SSD.INFO, obj, obj.SNPWeight=SNP_weight)
  return(out)
  SKAT::Close_SSD()
}




