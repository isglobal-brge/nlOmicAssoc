#' Function to fit a multivariable fractional polynomial model to a data.frame or an expressionSet
#'
#' Fits and mfp model for each row in data using mfp package
#'
#' @param data_m \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param df degrees of freedom for mfp, default 4, 1 would be linear
#' @param select alpha: sets the FP selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param select_adj select: sets the variable selection level for all predictors. Values for individual predictors may be changed via the fp function in the formula.
#' @param cores cores in case of parallelization (no windows)
#' @param imp_method imputation method for missing values in variables

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
#'
#' @export fit.mfp

fit.mfp <- function(data_m = data_m,
                    vars_df = vars_df,
                    df = 4,
                    select = 0.05,
                    select_adj = NULL,
                    cores = 1L,
                    ...)
{

  # Corrected select by correlations structure for a global 0.05(default):
  if( !is.null(select_adj) ){
    if(select_adj > 1 | select_adj < 0) {
      warning("Invalid select.adjust! It will be adjusted for a global 0.05 select value by the correlation matrix of vars_df")
      select_adj <- NULL
    }
  } else { #AKI no varia molt si calculem directament Bonferroni, el nombre de tests efectius redueix poc!
    cormat <- cor(vars_df, use = "pairwise.complete.obs")
    #if there is an na in cormat stop
    if (anyNA(cormat)) stop("Some NAs in correlation matrix, try to impute some of the missing values")
    M <- ncol(cormat)
    lambdas <- eigen(cormat,only.values = T)$values
    Vobs <- sum(((lambdas - 1)^2)) / (M - 1)
    Meff <- M * (1 - (M - 1) * Vobs / M^2) # number of effective tests
    select_adj <- select <- 1 - (1 - 0.05)^(1 / Meff) # corrected significance level (for global 0.05)
    rm(cormat, M, lambdas, Vobs, Meff)

    gc(reset = TRUE)
    select_adj<-adj.p(vars_df)
    print(select_adj)
  }

  # if (all(!rownames(vars_df) %in% colnames(data_m)) )
  # stop("vars_df rownames should match to data_m columns")

  data_m <- data_m[, rownames(vars_df)]

  # Formula structure for the models:
  formula <- as.formula(paste0("y ~ ", paste0("fp(", colnames(vars_df), ", df = ", df, " )", collapse = " + ")))

  #needed model
  fmodel <- function(y, vars_df, probe)
  {
    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))

    # Fit the model:
    mod <- try(mfp(formula, alpha=select, select = select_adj, rescale = FALSE, verbose = FALSE, data = data,...), TRUE)


    # Get model performances:
    if(class(mod)[1] == "try-error"){
      vars <- NA
      p <- NA
      cor2 <- NA
      aic <- NA
    } else {
      vars <- try(unique(as.character(str_sub(word(names(mod$coefficients)[-1],1), start = 1, end = -3))), TRUE)
      form<-try(mod$formula)
      p <- try(as.integer(sum(!is.na(vars))), TRUE)
      cor2 <- ifelse(p == 0, 0, try(cor(predict(mod), data$y, use = "pairwise.complete.obs")^2, TRUE))

      # Likelihood ratio test (LRT) :
      # f_null <- glm(data$y ~ 1)
      # f_full <- mod$fit
      # df.diff <- f_null$df.residual - f_full$df.residual
      # vals <- (sum(residuals(f_null)^2) - sum(residuals(f_full)^2))/sum(residuals(f_full)^2) * f_full$df.residual
      # p_LRT <- pchisq(vals, df.diff, lower.tail = FALSE)
      # #es pot fer directament:
      #rms:lrtest(f_null,f_full) #pero dona molt diferent, de fet la p mes alta
      # rm(f_null, f_full, df.diff, vals)

      formula.lin = as.formula(paste0("y"," ~ ", paste0(vars,collapse = " + ")))
      mod.lin <- try(lm(formula.lin, data = data),silent=TRUE)
      aic.mod<-AIC(mod$fit) #specific from mboost
      aic.lin<-AIC(mod.lin)
      aic<-matrix(c(aic.lin,aic.mod),2,1,
                  dimnames=list(c("mod.lin","mod"),"AIC"))
      gc(reset = TRUE)
    }


    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, p = p, Cor2 = cor2), stringsAsFactors = FALSE), TRUE)

    result <- try(list(table = taula, selected_vars = vars, final_formula = form, aic = aic), TRUE)
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


  p_adj <- p.adjust(as.double(as.character(taula$LRT_pval)), method = "BH")
  taula$adj_pval <- format.pval(p_adj)
  taula$LRT_pval <- format.pval(as.double(as.character(taula$LRT_pval)))
  taula <- taula[ aux <- base::order(as.double(taula$LRT_pval), decreasing = FALSE),]
  rownames(taula) <- NULL
  #taula should be a table and not a list
  selected_vars <-  lapply(results, function(x) x[2])[aux]
  final_formula <-  lapply(results, function(x) x[3])[aux]

  # results table, should be a special class and then filter, mirar MEAL
  return(list(table = taula, selected_vars = selected_vars, analysis_vars= colnames(vars_df), final_formula = final_formula, alpha_select = select_adj, aic = aic))
}


# ===============================================================================

