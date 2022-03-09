#' @title statistical analysis with SNP information in addition to meta data
#' @description users can choose different statistical methods to test significant variants
#' @name functional_analysis
#'
#' @param data a data with SNP and meta data
#' @param method statistical methods including "fisher", "chisq", "glm", "survfit", "coxph"
#' @param formula a symbolic descripton of the model to be fitted
#' @import survival
#' @import survminer
#' @import DescTools
#' @importFrom dplyr mutate arrange left_join group_by select summarise summarize
#' @importFrom magrittr `%>%`
#' @importFrom ggplot2 aes theme_bw theme element_blank
#' @importFrom grDevices pdf dev.off
#' @return a table containing varaibles and p-value
#' @export
#' @examples
#' dat<-read.csv(system.file("extdata", "merged_clinical_variant_data1_md.csv", package = "CLINSKAT"),
#'              stringsAsFactors = FALSE, encoding = "UTF-8", row.names = NULL, sep = ",")
#' sur_model<-as.formula("survival::Surv(Follow_up_months, Last_Status) ~ X3.142172064C.G")
#' KM_result<-functional_analysis(data=dat, method="survfit", formula = sur_model)
#'
#' dat$Gender<-as.factor(dat$Gender)
#' dat$X1.161298236A.G<-as.factor(dat$X1.161298236A.G)
#' mt<-table(dat$Gender, dat$X1.161298236A.G)
#' fisher_test<-functional_analysis(data=mt, method="fisher")
#' CAtest<-functional_analysis(data=mt, method="trend_test")
functional_analysis<-function(data, method=c("fisher", "chisq", "glm", "survfit", "coxph", "t_test", "trend_test"), formula){
  if (method == "fisher"){
    res<-fisher.test(data)$p.value
  }else if (method == "chisq"){
    res<-chisq.test(data)$p.value
  }else if (method == "glm"){
    res<-coef(summary(glm(formula, data = data)))[,4]
  }else if (method == "survfit"){
    res<-survminer::surv_pvalue(survminer::surv_fit(formula, data=data))
  }else if (method == "coxph"){
    res<-summary(survival::coxph(formula, data=data))$coefficients[,5]
  }else if (method == "t_test"){
    res<-t.test(dat$Age, dat$X1.161298236A.G)$p.value
  }else if (method == "trend_test"){
    res<-DescTools::CochranArmitageTest(data)$p.value
  }else stop("invalid type")

  #creat table for Manhattan plot
  if (any(colnames(dat)=="CHR")){
    don <- res %>%

      # Compute chromosome size
      group_by(CHR) %>%
      summarise(chr_len=max(BP)) %>%

      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%

      # Add this info to the initial dataset
      left_join(res, ., by=c("CHR"="CHR")) %>%

      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate( BPcum=BP+tot)

    axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    #make Manhattan plot
    output.pdf<-function(don){
      pdf(file = paste(getwd(), "Manhattanplot.pdf", sep="/"), onefile = TRUE)
      Manhattan_plot<-ggplot2::ggplot(don, aes(x=BPcum, y=-log10(P))) +

        # Show all points
        ggplot2::geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        ggplot2::scale_color_manual(values = rep(c("hotpink1", "skyblue", "tan1", "slateblue", "yellow2"), 5 )) +

        # custom X axis:
        ggplot2::scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, name = "Chromosome" ) +
        ggplot2::scale_y_continuous(limits = c(0, 6) ) +     # remove space between plot area and x axis

        #horizontal line
        ggplot2::geom_hline(yintercept=2, linetype="dashed", color = "dark grey") +

        # Custom the theme:
        theme_bw() +
        theme(
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
      output.pdf(pop_result)
      grDevices::dev.off()
      }
  }
  return(res)
}





