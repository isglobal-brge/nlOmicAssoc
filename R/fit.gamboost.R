#' Function to fit gam with boosting to a data.frame or an expressionSet
#'
#' GAM model is fitted for each row in data using mboost package
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param cores cores in case of parallelization (no windows)
#' @param df degrees of freedom to apply to model
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)
#'
#' @export fit.gamboost
#'
fit.gamboost<- function(data = data_m,
                    vars_df = vars_df,
                    df=3,
                    cores = 1L,
                    verbose = FALSE,
                    ...)
{

  # Formula structure for the models: only for numerical vars
  vars_class<-sapply(vars_df, function(x) class(x))
  vars_class_group<-split(names(vars_class),vars_class)
  formula <- as.formula(paste0("y ~ ", paste0(vars_class_group$numeric, collapse = " + ")))

  fmodel <- function(y, vars_df, probe)
  {

     if(verbose) print(paste("Analyzing probe var ", probe))

    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))
    xnames <- colnames(vars_df)

    mod <- try(gamboost(formula, data = data, dfbase = df),silent=T)

    # Get model performance
    if(class(mod)[1] == "try-error" ){

      vars <- NA
      vars_n <- NA
      cor2 <- NA
      p <- NA
      aic <- NA

    } else {

      vars = as.character(extract(mod, what = "variable.names"))
      vars_n = try(as.integer(sum(!is.na(vars))), TRUE)
      y_hat <- predict(mod)
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic<-round(AIC(mod),4) #specific from mboost
      gc(reset = TRUE)

    }

    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, vars_n = vars_n, aic = aic, Cor2 = cor2, p = p), stringsAsFactors = FALSE), TRUE)

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
