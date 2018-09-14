#' Function to calculate the adj_p based on the number of estimated effective tests
#'
#' Function to calculate the adj_p based on the number of estimated effective tests. Hidden to user
#'
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param select threshold for selection

adj.p <- function(vars_df, select = 0.05){


  cormat <- cor(vars_df, use = "pairwise.complete.obs")
  #if there is an na in cormat stop
  if (anyNA(cormat)) stop("Some NAs in correlation matrix, try to impute some of the missing values")
  M <- ncol(cormat)
  lambdas <- eigen(cormat,only.values = T)$values
  Vobs <- sum(((lambdas - 1)^2)) / (M - 1)
  Meff <- M * (1 - (M - 1) * Vobs / M^2) # number of effective tests
  adj_p <- select <- 1 - (1 - 0.05)^(1 / Meff) # corrected significance level (for global 0.05)
  rm(cormat, M, lambdas, Vobs, Meff)
  return(adj_p)
  gc(reset = TRUE)
}
