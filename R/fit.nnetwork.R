#' Function to fit neural network to a data.frame or an expressionSet
#'
#' Neural networks are adjusted for each row in data using package nnet
#'
#' @param data \code{ExpressionSet} or data.frame with samples as columns and observations as rows could be an ExpressionSet(modif!)
#' @param vars_df data.frame with sampes as rows and variables as columns
#' @param size nnet param: number of units in the hidden layer. Can be zero if there are skip-layer units
#' @param cores cores in case of parallelization (no windows)
#' @param verbose logical to verbose (comment) the steps of the function, default(FALSE)
#'
#' @export fit.nnetwork

fit.nnetwork <- function(data = data_m,
                    vars_df = vars_df,
                    size = 2L,
                    cores = 1L,
                    verbose = FALSE,
                    ...)

{

  # Formula structure for the models: only for numerical vars
  vars_class<-sapply(vars_df, function(x) class(x))
  vars_class_group<-split(names(vars_class),vars_class)
  formula <- as.formula(paste0("y ~ ", paste0(vars_class_group$numeric, collapse = " + ")))
  vars_df <- vars_df[,vars_class_group$numeric]

  fmodel <- function(y, vars_df, probe,i)
  {

    if(verbose) print(paste("Analyzing probe var ", probe))

    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))

    selnn <- bestNN(xdata = data[,-1], Y =data[,1], size = size, vars.drop.frac = 0.2, verbose = verbose)

    if (selnn == "") {

      mod = ""
      class(mod)="try-error"

    } else {

      mod <- try(nnet(as.formula(paste0( "y ~ ", selnn)),
                      data = data, size = size, decay = 0.001, maxit = 1000, linout = TRUE, trace = FALSE))
    }

    if(class(mod)[1] == "try-error" ){

      vars <- NA
      vars_n <- NA
      cor2 <- NA
      p <- NA
      aic <- NA

    } else {

      vars <- mod$coefnames
      vars_n = try(as.integer(sum(!is.na(vars))), TRUE)
      y_hat <- predict(mod,type="raw")
      ct <- cor.test(y_hat, data$y, use = "pairwise.complete.obs")
      cor2 <- ifelse(vars_n == 0, NA, try(round((ct$estimate)^2,4), TRUE))
      p <- ifelse(vars_n == 0, NA, try((ct$p.value)^2, TRUE))
      aic <- NA
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
