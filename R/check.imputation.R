#' Function to check whether data needs imputation: NO CAL
#'
#' @param vars_df data.frame to be checked for completeness (if it needs imputation)
#' @param perc number in (0,100) that is required for each variable to be complete
#' @import mice
#'
check.imputation <- function(vars_df,perc_imp=70,imp_method){
  # Removing samples and variables (in this order) that have >= 70% NAs
  pNA <- function(x){sum(is.na(x))/length(x)}

  #stop if too many missing values
  if (any(apply(vars_df,1,pNA)>(perc_imp/100))) stop("There are patients with more than 70% missing values, please remove or impute those cases and rerun function")
  if (any(apply(vars_df,2,pNA)>(perc_imp/100))) stop("There are variables with more than 70% missing values, please remove or impute those cases and rerun function")

  if (anyNA(vars_df)){
    #Imputation of missing values
    cat("Missing values will be imputed using function specified in imp.method of the mice package\n")
    cat("Imputing missing variables\n")
    vars_df_mice <- mice(vars_df, m=1, maxit = 10, method = 'fastpmm', seed = 500)
    vars_df <- mice::complete(vars_df_mice) #same object as original
    return(vars_df)
  }
}
