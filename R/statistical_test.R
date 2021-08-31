#' @title statistical analysis with SNP information in addition to meta data
#' @description users can choose different statistical methods to test significant variants
#' @name statistical_test
#'
#' @param data a data with SNP and meta data
#' @param method statistical methods including "fisher", "chisq", "lm", "survfit", "coxph"
#' @param formula a symbolic descripton of the model to be fitted
#' @import survival
#' @import survminer
#' @return a table containing varaibles and p-value
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata", "merged_clinical_variant_data1_md.csv", package = "skatvcf"),
#'              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' sur_model<-as.formula("survival::Surv(Follow_up_months, Last_Status) ~ X3.142172064C.G")
#' KM_result<-stat_test(data=dat, method="survfit", formula = sur_model)
#'
#' dat$Gender<-as.factor(dat$Gender)
#' dat$X1.161298236A.G<-as.factor(dat$X1.161298236A.G)
#' mt<-table(dat$Gender, dat$X1.161298236A.G)
#' fisher_test<-stat_test(data=mt, method="fisher")
statistical_test<-function(data, method=c("fisher", "chisq", "lm", "survfit", "coxph"), formula){
  if (method == "fisher"){
    res<-fisher.test(data)$p.value
  }else if (method == "chisq"){
    res<-chisq.test(data)$p.value
  }else if (method == "lm"){
    res<-summary(lm(formula, data = data))$coefficients[,4]
  }else if (method == "survfit"){
    res<-survminer::surv_pvalue(survminer::surv_fit(formula, data=data))
  }else if (method == "coxph"){
    res<-summary(survival::coxph(formula, data=data))$coefficients[,5]
  }else stop("invalid type")
  return(res)
}




