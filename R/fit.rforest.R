#' Function to fit rforest to a data.frame or an expressionSet
#'
#' Random forest models are fitted for each row in data using package randomForest
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param seed
#' @param maxnodes
#' @param ntree
#' @param cores cores in case of parallelization (no windows)
#' @param df degrees of freedom to apply to model
#'
#' @import stats
#' @import splines
#' @import parallel
#' @import randomForest
#' @import caret
#'
#' @export fit.rforest
#  no retorna formula!

fit.rforest <- function(data=data_m,
                    vars_df = vars_df,
                    seed = NULL,
                    maxnodes = 15,
                    ntree = 100,
                    cores = 1L,
                    ...)
{

  fmodel <- function(y, vars_df, probe) #canviar el nom per a que no sigui com mfp?
  {
    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))
    xnames <- colnames(vars_df)

    selrf <- selRF(xdata = data[,-1], Y = data[,1], ntree = ntree)
    formrforest<-as.formula(paste0("y ~ ", selrf$selected.model))

    mod <- try(randomForest(formrforest, data = data, ntree = ntree),silent=T) #posar ...?check param que passa do.models

    # Get model performance
    if(class(mod)[1] == "try-error" ){
      vars <- NA
      p <- NA
      cor2 <- NA
      AIC <- NA
    } else {
      vars = names(mod$importance[order(mod$importance[,1], decreasing = TRUE),1])
      p = try(as.integer(sum(!is.na(vars))), TRUE)
      cor2 <- ifelse(p == 0, 0, try(cor(mod$predicted, data$y, use = "pairwise.complete.obs")^2, TRUE))
      AIC <- "NA"
      gc(reset = TRUE)
    }


    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, p = p, Cor2 = cor2), stringsAsFactors = FALSE), TRUE)

    result <- try(list(table = taula, selected_vars = vars, final_formula = NA, aic = aic), TRUE)
    result
  }

  if(parallel_ind){
    results <- try(mcmapply(function(y, ny) fmodel(y = y, vars_df = vars_df, probe = ny),
                            y = apply(data_m, 1, list),
                            ny = rownames(data_m),
                            SIMPLIFY = FALSE, mc.cores = cores), TRUE)
  } else {
    results <- try(mapply(function(y, ny) fmodel(y = y, vars_df = vars_df, probe = ny),
                          y = apply(data_m, 1, list),
                          ny = rownames(data_m),
                          SIMPLIFY = FALSE), TRUE)
  }

  faux <- function(x){
    taula <- x[[1]]
    if(class(taula) == "try-error") taula <- data.frame(cbind(probe = taula$probe, p = NA, Cor2 = NA))
    taula
  }

  taula <- as.data.frame(do.call(rbind, lapply(results, faux)))
  rownames(taula) <- NULL
  #taula should be a table and not a list
  selected_vars <-  lapply(results, function(x) x[2])
  final_formula <-  lapply(results, function(x) x[3])

  return(list(table = taula, selected_vars = selected_vars, analysis_vars= colnames(vars_df), final_formula = final_formula, aic = aic))
}
