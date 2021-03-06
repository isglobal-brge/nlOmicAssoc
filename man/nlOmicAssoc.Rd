% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlOmicAssoc.R
\name{nlOmicAssoc}
\alias{nlOmicAssoc}
\title{Fits the specified multivariable model to a data.frame or an expressionSet}
\usage{
nlOmicAssoc(object, covars = NULL, model = c("partdsa", "exwassp", "mfp",
  "gam", "gamboost", "rforest", "nnetwork"), imp_method = "pmm",
  perc_imp = 70, cores = 1L, verbose = FALSE, ...)
}
\arguments{
\item{object}{numeric matrix, with samples as columns and observations as rows, or \code{ExpressionSet} or \code{SummarizedExperiment}}

\item{covars}{data.frame with samples as rows and variables as columns}

\item{imp_method}{imputation method for missing values in variables. Default is pmm}

\item{perc_imp}{numerical in [0,100], the maximum \% of a variable allowed to have NAs. Default is 70}

\item{cores}{cores in case of parallelization (no windows)}

\item{verbose}{logical to verbose (comment) the steps of the function. Default is FALSE}
}
\value{
nlAssoc

a file results with model, variables
}
\description{
Prepares data and fits specified model using specific call for each case
Returns an object on nlAssoc class
}
\examples{
to be built
}
