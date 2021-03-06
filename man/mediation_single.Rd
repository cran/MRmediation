% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediation_single.R
\name{mediation_single}
\alias{mediation_single}
\title{A causal mediation method with a single CpG site as the mediator}
\usage{
mediation_single(pheno, predictor, cpg, covariate, family = "gaussian")
}
\arguments{
\item{pheno}{A vector of continuous or binary phenotypes (class: numeric).}

\item{predictor}{A vector of values for the exposure variable (class: numeric).}

\item{cpg}{A vector of a CpG (class: numeric).}

\item{covariate}{A matrix of covariates. Each column is a covariate (class: data.frame).}

\item{family}{"gaussian" for continuous outcome or "binomial" for binary outcome.}
}
\value{
1. pval$TE:   total effect (TE) p-value \cr
        2. pval$DE:   direct effect (DE) p-value \cr
        3. pval$IE:   indirect effect (IE) p-value \cr
        4. pval_MX:   p-value for the association between methylation and exposure \cr
}
\description{
A causal mediation method with a single CpG site as the mediator
}
\examples{
################
### Examples ###
################
data("example_data")
predictor = data$exposure
cpg = data[,9] #any number in c(7:dim(data)[2])
covariates = subset(data, select=c("age","gender"))
# binary outcome
pheno_bin = data$pheno_bin
mediation_single(pheno_bin, predictor, cpg, covariate=covariates, family="binomial")
# continuous outcome
pheno_con = data$pheno_con
mediation_single(pheno_con, predictor, cpg, covariate=covariates, family="gaussian")
}
