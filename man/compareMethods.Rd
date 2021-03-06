% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareMethods.R
\name{compareMethods}
\alias{compareMethods}
\title{Compare nlAssoc objects obtained by nlOmicAssoc()}
\usage{
compareMethods(list.meth, filter = TRUE, param = "adj.p", comp = "below",
  threshold = 0.05)
}
\arguments{
\item{list.meth}{a NAMED list of nlAssoc objects as obtained by nlOmicAssoc(), usually applying different methods on the same data}

\item{filter}{if nlAssoc object has not been previously filtered it can be done here. Default is TRUE}

\item{param}{parameter of the table to filter with (vars_n, aic, Cor2, p or adj.p). Default is "adj.p"}

\item{comp}{comparison param comp threshold. Default is "below"}

\item{threshold}{numerical, the threshold to select associated variables related to the specified param. Default is 0.05}
}
\description{
Compare nlAssoc objects obtained by nlOmicAssoc()
}
