#' A causal mediation method with methylated region as the mediator
#' 
#' @name mediation
#' @aliases mediation
#' @param pheno A vector of continuous or binary phenotypes (class: numeric).
#' @param predictor A vector of values for the exposure variable (class: numeric).
#' @param region A matrix of CpGs in a region. Each column is a CpG (class: data.frame).
#' @param pos A vector of CpG locations from the defined region and they are from the same chromosome (class: integer).
#' @param order A value for the order of bspline basis. 1: constant, 2: linear, 3: quadratic and 4: cubic.
#' @param gbasis A value for the number of basis being used for functional transformation on CpGs.
#' @param covariate A matrix of covariates. Each column is a covariate (class: data.frame).
#' @param base "bspline" for B-spline basis or "fspline" for Fourier basis.
#' @param family "gaussian" for continuous outcome or "binomial" for binary outcome.
#' @import fda
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
#' region = data[,7:dim(data)[2]]
#' covariates = subset(data, select=c("age","gender"))
#' # binary outcome
#' pheno_bin = data$pheno_bin
#' mediation(pheno_bin, predictor, region, pos, covariate=covariates, order=4, 
#' gbasis=4, base="bspline", family="binomial")
#' # continuous outcome 
#' pheno_con = data$pheno_con
#' mediation(pheno_con, predictor, region, pos, covariate=covariates, order=4, 
#' gbasis=4, base="bspline", family="gaussian")
#' @export

#library(fda)
#library(MASS) #ginv

mediation <- function(pheno, predictor, region, pos, order, gbasis, covariate, base="bspline", family="gaussian")
  {
   predictor[is.na(predictor)] = 0
   region[is.na(region)] = 0
   covariate[is.na(covariate)] = 0
   idx     = is.na(pheno)
   pheno   = pheno[!idx]
   region    = region[!idx,]
   covariate = covariate[!idx,]
   dqr     = qr(region)
   index   = dqr$pivot[1:dqr$rank]
   region    = region[, index]
   pos     = pos[index]
   nsample = nrow(region)
   nsnp    = ncol(region)

   if(max(pos) > 1) pos = (pos - min(pos)) / (max(pos) - min(pos))

   betabasis  = create.bspline.basis(norder = 1, nbasis = 1)
   if (base == "bspline"){
      regionbasis  = create.bspline.basis(norder = order, nbasis = gbasis)
      } else if (base == "fspline"){
      regionbasis  = create.fourier.basis(c(0,1), nbasis = gbasis)
      }else { }

   region = as.matrix(region)
   B = eval.basis(pos, regionbasis)
   to_mul = ginv(t(B) %*% B) %*% t(B)
   U      = region %*% t( to_mul )
   J      = inprod(regionbasis, betabasis)
   UJ     = matrix( U %*% J, ncol = ncol(J) )
   region_transformed = as.vector(UJ)
   covariate = as.matrix(covariate)
   predictor = as.matrix(predictor)
   interaction = as.vector(predictor)*region_transformed

   ###########################
   if (family=="gaussian") {
       # region_transformed ~ covariate + predictor
       coef2_all = summary(lm(region_transformed ~ covariate + predictor))$coefficients 
       alpha_cov = coef2_all[c(1:(dim(covariate)[2]+1)),1]
       coef2_all_nocov = coef2_all[-c(1:(dim(covariate)[2]+1)), ]
       alpha_predictor = coef2_all_nocov[1]
       pval2 = coef2_all_nocov[4]
 
       pval <- list()
       # TE null hypothesis: b_predictor = b_region_transformed = b_interaction = 0
       test_part = cbind(predictor, region_transformed, interaction)
       untest_part = covariate
       fit = glm(pheno ~ untest_part + test_part, family="gaussian")
       pval$TE = anova(fit, test = "F")[3,6]
    
       # DE null hypothesis: b_predictor = b_interaction = 0
       test_part = cbind(predictor, interaction)
       untest_part = cbind(covariate, region_transformed)
       fit = glm(pheno ~ untest_part + test_part, family="gaussian")
       pval$DE = anova(fit, test = "F")[3,6]
    
       # IE null hypothesis: b_region_transformed = b_interaction = 0
       test_part = cbind(region_transformed, interaction)
       untest_part = cbind(covariate, predictor)
       fit = glm(pheno ~ untest_part + test_part, family="gaussian")
       pval$IE = anova(fit, test = "F")[3,6]
   }

   if (family=="binomial") {
       # region_transformed ~ covariate + predictor using controls
       region_transformed_c = region_transformed[which(pheno==0)]
       covariate_c = covariate[which(pheno==0),]
       covariate_cdqr = qr(covariate_c)
       covariate_cindex = covariate_cdqr$pivot[1:covariate_cdqr$rank]
       covariate_c = covariate_c[, covariate_cindex]
       diff = colnames(covariate)[!(colnames(covariate)%in%colnames(covariate_c))]
       predictor_c = as.vector(predictor)[which(pheno==0)]
       coef2_all = summary(lm(region_transformed_c ~ covariate_c + predictor_c))$coefficients 
       alpha_cov = coef2_all[c(1:(dim(covariate_c)[2]+1)),1]
       coef2_all_nocov = coef2_all[-c(1:(dim(covariate_c)[2]+1)), ]
       alpha_predictor = coef2_all_nocov[1]
       pval2 = coef2_all_nocov[4]
 
       pval <- list()
       # TE null hypothesis: b_predictor = b_region_transformed = b_interaction = 0
       test_part = cbind(predictor, region_transformed, interaction)
       untest_part = covariate
       fit = glm(pheno ~ untest_part + test_part, family="binomial")
       pval$TE = anova(fit, test = "Rao")[3,6]

       # DE null hypothesis: b_predictor = b_interaction = 0
       test_part = cbind(predictor, interaction)
       untest_part = cbind(covariate, region_transformed)
       fit = glm(pheno ~ untest_part + test_part, family="binomial")
       pval$DE = anova(fit, test = "Rao")[3,6]

       # IE null hypothesis: b_region_transformed = b_interaction = 0
       test_part = cbind(region_transformed, interaction)
       untest_part = cbind(covariate, predictor)
       fit = glm(pheno ~ untest_part + test_part, family="binomial")
       pval$IE = anova(fit, test = "Rao")[3,6]
 }
  return(list(pval=pval, pval_MX=pval2))
 }
