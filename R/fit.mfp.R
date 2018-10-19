#' Function to fit a multivariable fractional polynomial model to a data.frame or an expressionSet
#'
#' Fits and mfp model for each row in data using mfp package
#'
#' @param data_m \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param df degrees of freedom for mfp, default 4, 1 would be linear
#' @param alpha sets the FP selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param select select: sets the variable selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param cores cores in case of parallelization (no windows)
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)
#'
#' @examples to be built

fit.mfp <- function(data = data_m,
                    vars_df = vars_df,
                    df = 4L,
                    alpha = 0.05,
                    select = 0.05,
                    cores = 1L,
                    verbose = FALSE,
                    ...)

{

  # Formula structure for the models: only for numerical vars
  vars_numeric <- names(vars_df)
  formula <- as.formula(paste0("y ~ ", paste0("fp(", vars_numeric, ", df = ", df, " )", collapse = " + ")))

  fmodel <- function(y, vars_df, probe)
  {

    if(verbose) print(paste("Analyzing probe var ", probe))

    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))

    # Fit the model:
    mod <- try(mfp(formula, alpha=alpha, select = select, rescale = FALSE, verbose = FALSE, data = data), TRUE)

    # Get model performances:
    if(class(mod)[1] == "try-error"){

      vars <- NA
      vars_n <- NA
      form <- NA
      cor2 <- NA
      p <- NA
      aic <- NA

    } else {

      vars <- try(unique(as.character(str_sub(word(names(mod$coefficients)[-1],1), start = 1, end = -3))), TRUE) #aix? nom?s extrau totes les vars continues
      vars_n <- try(as.integer(sum(!is.na(vars))), TRUE)
      form <- try(mod$formula)
      y_hat <- predict(mod)
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic <- round(mod$aic,4)
      gc(reset = TRUE)

    }

    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, vars_n = vars_n, aic = aic, Cor2 = cor2, p = p), stringsAsFactors = FALSE), TRUE)

    result <- try(list(table = taula, selected_vars = vars, final_formula = form), TRUE)
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


# ===============================================================================

