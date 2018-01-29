#' Function to fit gam with boosting to a data.frame or an expressionSet
#'
#' GAM model is fitted for each row in data using mboost package
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param cores cores in case of parallelization (no windows)
#' @param df degrees of freedom to apply to model
#'
#' @import stats
#' @import splines
#' @import parallel
#' @import mboost #hi ha un altre paquet GAMboost CHECK
#'
#' @export fit.gamboost
#'
fit.gamboost<- function(data = data_m,
                    vars_df = vars_df,
                    df=3,
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

    formulaboost <- as.formula(paste0("y ~ ", paste0(colnames(data)[-1], collapse =" + ")))
    mod <- try(gamboost(formulaboost, data = data, dfbase = 3),silent=T)

    # Get model performance
    if(class(mod)[1] == "try-error" ){
      vars <- NA
      p <- NA
      cor2 <- NA
      p_LRT <- NA
    } else {
      vars = variables = as.character(extract(mod, what = "variable.names"))
      p = try(as.integer(sum(!is.na(vars))), TRUE)
      cor2 <- ifelse(p == 0, 0, try(cor(predict(mod), data$y, use = "pairwise.complete.obs")^2, TRUE))

      formula.lin = as.formula(paste0("y"," ~ ", paste0(vars,collapse = " + ")))
      mod.lin <- try(lm(formula.lin, data = data),silent=TRUE)
      #compare both models using AIC
      aic.mod<-AIC(mod) #specific from mboost
      aic.lin<-AIC(mod.lin)
      aic<-matrix(c(aic.lin,aic.mod),2,1,
                  dimnames=list(c("mod.lin","mod"),"AIC"))

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
   selected_vars <-  lapply(results, function(x) x[2])
  final_formula <-  lapply(results, function(x) x[3])

  return(list(table = taula, selected_vars = selected_vars, analysis_vars= colnames(vars_df), final_formula = final_formula, aic = aic))
}
