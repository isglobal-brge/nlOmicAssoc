#' Fits the specified multivariable model to a data.frame or an expressionSet
#'
#' Prepares data and fits specified model using specific call for each case
#' Returns an object on nlAssoc class
#'
#' @param set numeric matrix, or \code{ExpressionSet} or \code{SummarizedExperiment} with samples as columns and observations as rows could be an ExpressionSet(modif!). rownames!
#' @param vars_df data.frame with samples as rows and variables as columns
#' @param df degrees of freedom for mfp, default 4, 1 would be linear
#' @param select alpha: sets the FP selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param select_adj select: sets the variable selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param cores cores in case of parallelization (no windows)
#' @param imp_method imputation method for missing values in variables
#' @param perc_imp numerical in [0,100], the maximum % of a variable allowed to have Nas
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)
#'
#' @return nlAssoc
#'
#' @import stats
#' @import Biobase
#' @import parallel
#' @import mice
#' @import stringr
#' @import caret
#' @import splines
#' @import mfp
#' @import partDSA
#' @import gam
#' @import mboost
#' @import randomForest
#' @import nnet
#' @import NeuralNetTools
#' @import SummarizedExperiment
#'
#' @return a file results with model, variables
#' @examples to be built

#' @export fit.model


fit.model <- function(set,
                      vars_df=NULL,
                      model = c("dsa","ewas","mfp","gam","gamboost","rforest","nnetwork"),
                      imp_method = 'pmm',
                      perc_imp = 70,
                      cores = 1L,
                      verbose = FALSE,
                      ...)
{

#########################
##Argument checking #####
#########################
  call<-match.call() #to be returned at the end
  #globalVariables("data_m")

  if (is(set, "matrix")){

    #ARGCHECK: numeric data_m
    if (!is.numeric(set) ) stop("set should be a numeric matrix, an ExpressionSet or a SummarizedExperiment")
    data_m <- set

    #ARGCHECK: vars_df as data.frame
    if (!is.data.frame(vars_df) ) stop("vars_df should be a data.frame")

    #ARGCHECK: vars_df common patients of data_m
    if (all(!rownames(vars_df) %in% colnames(data_m)))  stop("vars_df rownames should match to set columns")

  } else if (is(set, "ExpressionSet")){

    data_m <- Biobase::exprs(set)
    vars_df <- Biobase::pData(set)

  } else if (is(set, "SummarizedExperiment")){

    data_m <- SummarizedExperiment::assay(set)
    vars_df <- SummarizedExperiment::colData(set)

  } else {

    stop("set must be a numeric matrix, an ExpressionSet or a SummarizedExperiment")
  }

  vars_class<-sapply(vars_df, function(x) class(x))
  vars_class_group<-split(names(vars_class),vars_class)

  #ARGCHECK: numeric variables
  if(is.null(vars_class_group$numeric))
    stop("variables must contain at least one continuous variable")

  #ARGCHECK: model in list
  if (!model %in% c("dsa","ewas","mfp","gam","gamboost","rforest","nnetwork"))
    stop("model has to be one of the following: dsa, ewas, mfp, gam, gamboost, rforest, nnetwork")

  #ARGCHECK: cores numeric
  if (!is.numeric(cores)) stop("cores should be numeric")


  #ARGCHECK: perc_imp numeric in [0,100]
  if (!is.numeric(perc_imp) | perc_imp<0 | perc_imp>100)
    stop("perc_imp should be numeric and between 0 and 100")

#########################
##NAs and imputation ####
#########################

  # Removing samples and variables (in this order) that have >=perc_imp NAs

  pNA <- function(x){sum(is.na(x))/length(x)}

  #if there are rows or cols with more than perc_imp NAs notify and stop, put perc
  if (any(apply(vars_df,2,pNA)>(perc_imp/100)))

    stop(paste0("There are variables with more than ",perc_imp,"% missing values, please remove or impute those cases and rerun function"))

  if (any(apply(vars_df,1,pNA)>(perc_imp/100)))

    stop(paste0("There are patients with more than ",perc_imp,"% missing values, please remove or impute those cases and rerun function"))

  #Imputation of missing values
  if (anyNA(vars_df)){

    cat("Missing values will be imputed using function specified in imp.method of the mice package\n")
    cat("Imputing missing variables through package mice\n")
    vars_df_mice <- mice::mice(vars_df, m=1, maxit = 10, method = imp_method, seed = 500) #a 3.4.2 ja no existeix fastpmm???
    vars_df_mice <- mice::complete(vars_df_mice) #same object as original excepf for rownames
    rownames(vars_df_mice)<-rownames(vars_df)
    vars_df <- vars_df_mice
    rm(vars_df_mice)

  }

  data_m <- data_m[, rownames(vars_df)]

#########################
##### model call ########
#########################

  if (model=="dsa"){

    res <- fit.partDSA(data_m, vars_df, df = 3L, cores = cores, verbose=verbose, ...) #not checked

  } else if (model=="ewas"){

    res <- fit.ewas(data_m, vars_df, df = 3L, cores = cores, verbose=verbose, ...)

  } else if (model=="mfp"){

    res<-fit.mfp(data_m, vars_df, df = 4L, cores = cores, verbose=verbose, ...)

  } else if (model=="gam"){

    res <- fit.gam(data_m, vars_df, cores = cores,df=4L, verbose=verbose, ...)

  } else if (model=="gamboost"){

    res <- fit.gamboost(data_m, vars_df, df = 3L, cores = cores, verbose=verbose, ...) #si!

  } else if (model=="rforest"){

    res <- fit.rforest(data_m, vars_df, seed = NULL, maxnodes = 15L, ntree = 100L, cores = cores, verbose=verbose, ...)

  } else if (model=="nnetwork"){

    res <- fit.nnetwork(data_m, vars_df, size = 2L,cores = cores, verbose=verbose, ...)

  }

  faux <- function(x){

    taula <- x[[1]]
    if(class(taula) == "try-error") taula <- data.frame(cbind(probe = taula$probe, vars_n = NA, Cor2 = NA, p = NA, aic = NA))
    taula

  }

  taula <- as.data.frame(do.call(rbind, lapply(res, faux)))
  rownames(taula) <- NULL
  taula$vars_n <- as.numeric(taula$vars_n)
  taula$Cor2<- as.numeric(taula$Cor2)
  taula$p <- as.numeric(taula$p)
  taula$adj.p <- p.adjust(taula$p)
  taula$aic<- as.numeric(taula$aic)

  selected_vars <-  lapply(res, function(x) x[2])
  final_formula <-  lapply(res, function(x) x[3])

  results <- list(table = taula, selected_vars = selected_vars, analysis_vars= colnames(vars_df), final_formula = final_formula)

  results$call <- call
  results$set.def <- data_m
  results$vars_df.def <- vars_df
  class(results)<-"nlAssoc"

  return(results)

}


# ===============================================================================



