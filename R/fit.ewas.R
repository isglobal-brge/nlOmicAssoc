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
#'
#' @import stats
#' @import splines
#' @import parallel
#'
#' @export fit.ewas


fit.ewas <- function(data = data_m,
                     vars_df = vars_df,
                     df = 3,
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

    #check si es ny i si cal passar-ho a la formula o es pot deixar y gen?ric
    formulesuni <- lapply(xnames, function(x) as.formula(paste0("y", " ~ ", paste0("ns(", x, ", df = ",df,")"))))
    #les vars que tenen poca variabilitat peten...i=18 per exemple, poso un try!
    modelsuni <- lapply(formulesuni, function(f) try(lm(formula = f, data = data),silent=T))

    #extract
    # p-valors
    p <- unlist(lapply(modelsuni, function(x) try(anova(object = x, test = "LRT")$"Pr(>F)"[1],silent=T)))
    padj <- p.adjust(p, method = "BH") # ajustats
    names(padj) <- names(p) <- names(formulesuni) <-  names(modelsuni) <- xnames

    # Variables significatives: sigincomes <- names(p)[p <= .05]
    sigincomes <- names(padj)[!is.na(padj) & padj <= .05]
    #pot ser que no n'hi hagi cap!
    sigincomes

    #model multivariant!
    # Formula model multivariant: revisar indentat de model.ewas
    #if no linear model
    if(length(sigincomes) == 0) { formula <- as.formula(paste("y"," ~ 1"))
    } else {  formula <- as.formula(paste0("y"," ~ ", paste0("ns(", sigincomes, ", df = ",df,")", collapse = " + "))) }

    # Model multivariant
    mod <- try(lm(formula, data = data),silent=TRUE)

  # Get model performance
    if(class(mod)[1] == "try-error" | length(sigincomes) == 0){
      vars <- NA
      p <- NA
      cor2 <- NA
      aic <- NA

    } else {
      vars <- sigincomes
      form<-try(mod$formula)
      p <- try(as.integer(sum(!is.na(vars))), TRUE)
      cor2 <- ifelse(p == 0, 0, try(cor(predict(mod), data$y, use = "pairwise.complete.obs")^2, TRUE))

      formula.lin = as.formula(paste0("y"," ~ ", paste0(sigincomes,collapse = " + ")))
      mod.lin <- try(lm(formula.lin, data = data),silent=TRUE)
      #compare both models using AIC
      #aic<-AIC(mod.lin,mod) #data.frame
      aic.mod<-AIC(mod) #specific from mboost
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
  rownames(taula) <- NULL
  selected_vars <-  lapply(results, function(x) x[2])
  final_formula <-  lapply(results, function(x) x[3])

  return(list(table = taula, selected_vars = selected_vars, analysis_vars= colnames(vars_df), final_formula = final_formula, aic = aic))
}




