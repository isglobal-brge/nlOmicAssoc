#' Function to fit ewas to a data.frame or an expressionSet
#'
#' For each variable (for each row in data) a natural cubic spline with 3 knots is fitted using splines package.
#' Only significative variables obtained after multiple comparison are then included in the
#' multivariate model.
#' Models are fitted using lm and simple least squares.
#' Model is compared to linear model using AIC
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param cores cores in case of parallelization (no windows)
#' @param df degrees of freedom to apply to model (number of knots in the natural cubic splines)
#' @param fdr threshold to select linear significative vars in terms of FDR adjusted p.value
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)
#'
#' @export fit.ewas


fit.ewas <- function(data = data_m,
                     vars_df = vars_df,
                     df = 3L,
                     fdr = .05,
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


    formulesuni <- lapply(xnames, function(x) as.formula(paste0("y", " ~ ", paste0("ns(", x, ", df = ",df,")"))))

    modelsuni <- lapply(formulesuni, function(f) try(lm(formula = f, data = data),silent=T))

    #extract p-values
    p <- unlist(lapply(modelsuni, function(x) try(anova(object = x, test = "LRT")$"Pr(>F)"[1],silent=T)))

    padj <- p.adjust(p, method = "fdr")
    names(padj) <- names(p) <- names(formulesuni) <-  names(modelsuni) <- xnames

    # Significative variables: sigincomes <- names(p)[p <= .05]
    sigincomes <- names(padj)[!is.na(padj) & padj <= fdr]

    if(length(sigincomes) == 0) {

      formula <- as.formula(paste("y"," ~ 1"))

    } else {

      formula <- as.formula(paste0("y"," ~ ", paste0("ns(", sigincomes, ", df = ",df,")", collapse = " + ")))

    }

    # Multivariate model
    mod <- try(lm(formula, data = data),silent=TRUE)

  # Get model performance
    if(class(mod)[1] == "try-error" | length(sigincomes) == 0){

      vars <- NA
      vars_n <- NA
      form <- NA
      cor2 <- NA
      p <- NA
      aic <- NA

    } else {

      vars <- sigincomes
      vars_n <- try(as.integer(sum(!is.na(vars))), TRUE)
      form <- try(formula(mod))
      y_hat <- predict(mod)
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic<-round(AIC(mod),4) #specific from mboost
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




