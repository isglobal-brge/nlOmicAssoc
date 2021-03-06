% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.ewas.R
\name{fit.ewas}
\alias{fit.ewas}
\title{Function to fit ewas to a data.frame or an expressionSet}
\usage{
fit.ewas(data = data_m, vars_df = vars_df, df = 3L, fdr = 0.05,
  cores = 1L, verbose = FALSE, ...)
}
\arguments{
\item{data}{\code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)}

\item{vars_df}{data.frame with sampes as rows and variables as columns}

\item{df}{degrees of freedom to apply to model (number of knots in the natural cubic splines)}

\item{fdr}{threshold to select linear significative vars in terms of FDR adjusted p.value}

\item{cores}{cores in case of parallelization (no windows)}

\item{verbose}{logical to verbose (comment) the steps of the function, default(FALSE)}
}
\description{
For each variable (for each row in data) a natural cubic spline with 3 knots is fitted using splines package.
Only significative variables obtained after multiple comparison are then included in the
multivariate model.
Models are fitted using lm and simple least squares.
Model is compared to linear model using AIC
}
