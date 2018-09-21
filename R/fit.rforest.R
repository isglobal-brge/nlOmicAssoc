#' Function to fit rforest to a data.frame or an expressionSet
#'
#' Random forest models are fitted for each row in data using package randomForest
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param cores cores in case of parallelization (no windows)
#' @param seed set a see, for reproducibility
#' @param maxnodes maximum number of terminal nodes trees in the forest can have
#' @param ntree number of trees to grow
#' @param df degrees of freedom to apply to model
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)
#'
#' @export fit.rforest


fit.rforest <- function(data=data_m,
                    vars_df = vars_df,
                    seed = NULL,
                    maxnodes = 15L,
                    ntree = 2000,
                    cores = 1L,
                    verbose = FALSE,
                    ...)
{

  fmodel <- function(y, vars_df, probe)
  {

    if(verbose) print(paste("Analyzing probe var ", probe))

    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))
    xnames <- colnames(vars_df)

    # Best model
    selrf <- bestRF(xdata = data[,-1], Y = data[,1], ntree = ntree ,vars.drop.frac = 0.2,
                    verbose = verbose)

    formula<-as.formula(paste0("y ~ ", selrf))
    mod <- try(randomForest(formula, data = data, ntree = ntree),silent=T)

    # Get model performance
    if(class(mod)[1] == "try-error" ){

      vars <- NA
      vars_n <- NA
      cor2 <- NA
      p <- NA
      aic <- NA

    } else {

      vars = names(mod$importance[order(mod$importance[,1], decreasing = TRUE),1])
      vars_n = try(as.integer(sum(!is.na(vars))), TRUE)
      y_hat <- mod$predicted
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic <- NA #potser podr?em retornar oob.err
      gc(reset = TRUE)

    }

    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, vars_n = vars_n, aic =aic, Cor2 = cor2, p = p), stringsAsFactors = FALSE), TRUE)

    result <- try(list(table = taula, selected_vars = vars, final_formula = NA), TRUE)
    result

  }

  if(cores>1){

    results <- try(mcmapply(function(y, ny) fmodel(y = y, vars_df = vars_df, probe = ny),
                            y = apply(data, 1, list),
                            ny = rownames(data),
                            SIMPLIFY = FALSE, mc.cores = cores), TRUE)

  } else {

    results <- try(mapply(function(y, ny) fmodel(y = y, vars_df = vars_df, probe = ny),
                          y = apply(data, 1, list),
                          ny = rownames(data),
                          SIMPLIFY = FALSE), TRUE)

  }

  return(results)

}
