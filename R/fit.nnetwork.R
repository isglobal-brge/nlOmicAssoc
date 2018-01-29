#' Function to fit neural network to a data.frame or an expressionSet
#'
#' Neural networks are adjusted for each row in data using package nnet
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param size nnet param: number of units in the hidden layer. Can be zero if there are skip-layer units
#' @param cores cores in case of parallelization (no windows)
#' @import stats
#' @import splines
#' @import parallel
#' @import nnet
#' @export fit.nnetwork
#  no retorna formula!
#REVISAR

fit.nnetwork <- function(data = data_m,
                    vars_df = vars_df,
                    size = 2,
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

    formulannet <- as.formula(paste0("y ~ ", paste0(colnames(data)[-1], collapse =" + ")))
    #check fn selnet de les simulacions, per a escollir el millor model
    mod <- try(nnet(formulannet, data = data,size=size),silent=T)

    # Get model performance
    if(class(mod)[1] == "try-error" ){
      vars <- NA
      p <- NA
      cor2 <- NA
      AIC <- NA
    } else {
      vari<-varImp(mod)
      vars = names(vari[order(vari[,1], decreasing = TRUE),])
      p = try(as.integer(sum(!is.na(vars))), TRUE)
      cor2 <- ifelse(p == 0, 0, try(cor(predict(mod,data[,1],type="raw"), data$y, use = "pairwise.complete.obs")^2, TRUE)) #no xuta
      AIC <- NA
      gc(reset = TRUE)
    }


    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, p = p, Cor2 = cor2), stringsAsFactors = FALSE), TRUE)

    result <- try(list(table = taula, selected_vars = vars, final_formula = NA), TRUE)
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
