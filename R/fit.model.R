#' Fits the specified multivariable model to a data.frame or an expressionSet
#'
#' Prepares data and fits specified model using specific call for each case
#'
#' @param data_m \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!). rownames!
#' @param vars_df data.frame with samples as rows and variables as columns
#' @param df degrees of freedom for mfp, default 4, 1 would be linear
#' @param select alpha: sets the FP selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param select_adj select: sets the variable selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param cores cores in case of parallelization (no windows)
#' @param imp_method imputation method for missing values in variables
#' @param perc_imp
#'
#' @import stats
#' @import mice
#' @import stringr
#' @import mfp
#' @import Biobase
#' @import parallel
#' @import stats
#'
#' @return a file results with model, variables
#' @examples to be built

#' @export fit.model
#'

fit.model <- function(data_m,
                      vars_df,
                      model = c("dsa","ewas","mfp","gam","gamboost","rforest","nnetwork"),
                      imp_method = 'fastpmm',
                      perc_imp = 70,
                      cores = 1L,
                      ...)
{
  #incloure df com a param gral? es fa servir a mfp ewas,

  #Argument Checkings
  # # Corrected select by correlations structure for a global 0.05(default): ->fit.mfp
  # if( !is.null(select_adj) ){
  #   if(select_adj > 1 | select_adj < 0) {
  #     warning("Invalid select.adjust! It will be adjusted for a global 0.05 select value by the correlation matrix of vars_df")
  #     select_adj <- NULL
  #   }
  #
  # }

  #match.call match.arg match.fun CHECK
  call <- match.call()

  #ARGCHECK: data_m data.frame
  #ARGCHECK: numeric data_m
  #ARGCHECK: numeric vars_df
  #ARGCHECK: model in list
  #ARGCHECK: cores numeric
  #ARGCHECK: imp_method='fastpmm', other?
  #ARGCHECK: perc_imp numeric in [0,100]
  #ARGCHECK: vars_df data.frame with common patients of data_m
  if (all(!rownames(vars_df) %in% colnames(data_m)) )
    stop("vars_df rownames should match to data_m columns")

  data_m <- data_m[, rownames(vars_df)] #check si ha de ser aix? per tots el m?todes

  # Removing samples and variables (in this order) that have >= 70% NAs
  #Function check.imputation? no cal pq aix? ho hem de fer sempre
  #ci<-check.imputation(vars_df=vars_df,perc_imp=perc_imp,imp_method=imp_method)

  pNA <- function(x){sum(is.na(x))/length(x)}

  #we first remove patients having too many missing values too
  if (any(apply(vars_df,1,pNA)>(perc_imp/100))) stop("There are patients with more than 70% missing values, please remove or impute those cases and rerun function")
  if (any(apply(vars_df,2,pNA)>(perc_imp/100))) stop("There are variables with more than 70% missing values, please remove or impute those cases and rerun function")

  if (anyNA(vars_df)){
    #Imputation of missing values
    cat("Missing values will be imputed using function specified in imp.method of the mice package\n")
    cat("Imputing missing variables\n")
    vars_df_mice <- mice::mice(vars_df, m=1, maxit = 10, method = 'pmm', seed = 500) #a 3.4.2 ja no existeix fastpmm???
    vars_df <- mice::complete(vars_df_mice) #same object as original
  }

  if( is.null(select_adj) ){ #AKI no varia molt si calculem directament Bonferroni, el nombre de tests efectius redueix poc!
    #espec?fic de mfp, per tant es pot posar directament a la fn fit.mfp
    # cormat <- cor(vars_df, use = "pairwise.complete.obs")
    # #if there is an na in cormat stop
    # if (anyNA(cormat)) stop("Some NAs in correlation matrix, try to impute some of the missing values")
    # M <- ncol(cormat)
    # lambdas <- eigen(cormat,only.values = T)$values
    # Vobs <- sum(((lambdas - 1)^2)) / (M - 1)
    # Meff <- M * (1 - (M - 1) * Vobs / M^2) # number of effective tests
    # select_adj <- select <- 1 - (1 - 0.05)^(1 / Meff) # corrected significance level (for global 0.05)
    # rm(cormat, M, lambdas, Vobs, Meff)
    #
    # gc(reset = TRUE)
    select_adj<-adj.p(vars_df)
    print(select_adj)
  }


  #AKI la crida a cada model
  if (model=="dsa"){

    res <- fit.partDSA(data_m, vars_df, df = 3, cores = 1L,...)

  } else if (model=="ewas"){

    res <- fit.ewas(data_m, vars_df, df = 3, cores = 1L,...)

  } else if (model=="mfp"){

    res<-fit.mfp(data_m, vars_df, df = 4, select = 0.05, select_adj = NULL, cores = 1L,...)
      #filter results in terms of a param, cor?
  } else if (model=="gam"){

    res <- fit.gam(data_m, vars_df, df = 3, cores = 1L,...)

  } else if (model=="gamboost"){

    res <- fit.gamboost(data_m, vars_df, df = 3, cores = 1L,...)

  } else if (model=="rforest"){

    res <- fit.rforest(data_m, vars_df, seed = NULL, maxnodes = 15, ntree = 100, cores = 1L,...)

  } else if (model=="nnetwork"){

    res <- fit.nnetwork(data_m, vars_df, size = 2, cores = 1L,...)

  }

  #modificar resultats? SI: cal seleccionar-los
  res

}


# ===============================================================================

