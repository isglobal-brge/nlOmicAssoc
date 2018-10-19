#' Function to fit gam to a data.frame or an expressionSet
#'
#' GAM models are fitted for each row in data using lm and backfitting using bruto function in mda package.
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param cores cores in case of parallelization (no windows)
#' @param df degrees of freedom to apply to model
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)


fit.gam <- function(data = data_m,
                    vars_df = vars_df,
                    #vars_class = vars_class,
                    df=4L,
                    cores = 1L,
                    verbose = FALSE,
                    ...)
{

  # Formula structure for the models: only for numerical vars
  # vars_class<-sapply(vars_df, function(x) class(x))
  # vars_class_group <- split(names(vars_class),vars_class)
  # vars_numeric <- vars_class_group$numeric
  vars_numeric <- names(vars_df)
  formula <- as.formula(paste0("y ~ ", paste0("s(", vars_numeric, ", df = ", df, " )", collapse = " + ")))

  fmodel <- function(y, vars_df, probe)
  {

    if(verbose) print(paste("Analyzing probe var ", probe))

    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y",vars_numeric)
    res.df = T

    if(nrow(data)<length(vars_numeric)*df) {

      print("num of observations is higher than df")
      res.df = FALSE

    }
    mod <- try(gam(formula,data=data ),silent=T)

    # Get model performance
    if(class(mod)[1] == "try-error" ){

      vars <- NA
      vars_n <- NA
      form = NA
      cor2 <- NA
      p <- NA
      aic <- NA

    } else if (!(res.df)){

      vars <- NA
      vars_n <- NA
      form = NA
      y_hat <- predict(mod)
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic <- round(mod$aic,4)

    } else {

      vars = try(vars_numeric[anova(mod)[,3][-1] < 0.05],TRUE)
      vars_n = try(as.integer(sum(!is.na(vars))), TRUE)
      form <- as.formula(paste0("y ~ ", paste0("s(", vars, ", df = ", df, " )", collapse = " + ")))
      y_hat <- predict(mod)
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic <- round(mod$aic,4)

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
