#' A causal mediation method with a single CpG site as the mediator
#'
#' @name mediation_single
#' @aliases mediation_single
#' @param pheno A vector of continuous or binary phenotypes (class: numeric).
#' @param predictor A vector of values for the exposure variable (class: numeric).
#' @param cpg A vector of a CpG (class: numeric).
#' @param covariate A matrix of covariates. Each column is a covariate (class: data.frame).
#' @param family "gaussian" for continuous outcome or "binomial" for binary outcome.
#' @importFrom MASS ginv
#' @importFrom stats glm lm pf anova var
#' @return 1. pval$TE:   total effect (TE) p-value \cr
#'         2. pval$DE:   direct effect (DE) p-value \cr
#'         3. pval$IE:   indirect effect (IE) p-value \cr
#'         4. pval_MX:   p-value for the association between methylation and exposure \cr
#' @examples
#' ################
#' ### Examples ###
#' ################
#' data("example_data")
#' predictor = data$exposure
#' cpg = data[,9] #any number in c(7:dim(data)[2])
#' covariates = subset(data, select=c("age","gender"))
#' # binary outcome
#' pheno_bin = data$pheno_bin
#' mediation_single(pheno_bin, predictor, cpg, covariate=covariates, family="binomial")
#' # continuous outcome
#' pheno_con = data$pheno_con
#' mediation_single(pheno_con, predictor, cpg, covariate=covariates, family="gaussian")
#' @export

mediation_single <- function(pheno, predictor, cpg, covariate, family="gaussian")
  {
   covariate = as.matrix(covariate)
   predictor = as.vector(as.matrix(predictor)[,1])
   interaction = predictor*cpg

   ###########################
   if (family=="gaussian") {
       # cpg ~ covariate + predictor
       coef2_all = summary(lm(cpg ~ covariate + predictor))$coefficients 
       alpha_cov = coef2_all[c(1:(dim(covariate)[2]+1)),1]
       coef2_all_nocov = coef2_all[-c(1:(dim(covariate)[2]+1)), ]
       alpha_predictor = coef2_all_nocov[1]
       pval2 = coef2_all_nocov[4]
 
       pval <- list()
       # TE null hypothesis: b_predictor = b_cpg = b_interaction = 0
       test_part = cbind(predictor, cpg, interaction)
       untest_part = covariate
       fit = glm(pheno ~ untest_part + test_part, family="gaussian")
       pval$TE = anova(fit, test = "F")[3,6]
    
       # DE null hypothesis: b_predictor = b_interaction = 0
       test_part = cbind(predictor, interaction)
       untest_part = cbind(covariate, cpg)
       fit = glm(pheno ~ untest_part + test_part, family="gaussian")
       pval$DE = anova(fit, test = "F")[3,6]
    
       # IE null hypothesis: b_cpg = b_interaction = 0
       test_part = cbind(cpg, interaction)
       untest_part = cbind(covariate, predictor)
       fit = glm(pheno ~ untest_part + test_part, family="gaussian")
       pval$IE = anova(fit, test = "F")[3,6]
   }

   if (family=="binomial") {
       # cpg ~ covariate + predictor using controls
       cpg_c = cpg[which(pheno==0)]
       covariate_c = covariate[which(pheno==0),]
       covariate_cdqr = qr(covariate_c)
       covariate_cindex = covariate_cdqr$pivot[1:covariate_cdqr$rank]
       covariate_c = covariate_c[, covariate_cindex]
       diff = colnames(covariate)[!(colnames(covariate)%in%colnames(covariate_c))]
       predictor_c = as.vector(predictor)[which(pheno==0)]
       coef2_all = summary(lm(cpg_c ~ covariate_c + predictor_c))$coefficients 
       alpha_cov = coef2_all[c(1:(dim(covariate_c)[2]+1)),1]
       coef2_all_nocov = coef2_all[-c(1:(dim(covariate_c)[2]+1)), ]
       alpha_predictor = coef2_all_nocov[1]
       pval2 = coef2_all_nocov[4]
 
       pval <- list()
       # TE null hypothesis: b_predictor = b_cpg = b_interaction = 0
       test_part = cbind(predictor, cpg, interaction)
       untest_part = covariate
       fit = glm(pheno ~ untest_part + test_part, family="binomial")
       pval$TE = anova(fit, test = "Rao")[3,6]

       # DE null hypothesis: b_predictor = b_interaction = 0
       test_part = cbind(predictor, interaction)
       untest_part = cbind(covariate, cpg)
       fit = glm(pheno ~ untest_part + test_part, family="binomial")
       pval$DE = anova(fit, test = "Rao")[3,6]

       # IE null hypothesis: b_cpg = b_interaction = 0
       test_part = cbind(cpg, interaction)
       untest_part = cbind(covariate, predictor)
       fit = glm(pheno ~ untest_part + test_part, family="binomial")
       pval$IE = anova(fit, test = "Rao")[3,6]
 }
   
   return(list(pval=pval, pval_MX=pval2))
 }
