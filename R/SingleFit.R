
# ===========================================================================================
# Funcio per ajustar un model MFP a la expressio indicada amb totes les exposicions
# ===========================================================================================

SingleFit <- function(expset, expo, probe, df = 4, select = 0.05, select_adj = NULL)
{
  # Checkings:
  
  if(is.character(probe)){
    probe_ind <- grep(probe, rownames(exprs(expset)))
    if(length(probe_ind) > 1)
      warning("\nMultiple matching for this probe name.\nFirst element has been selected: ", 
              rownames(exprs(expset))[probe_ind], "\n")
    else if (length(probe_ind) == 0)
      stop("\nNo matchings for this probe!\n")
  } else if(is.integer(probe) | is.double(probe)){
    if(length(probe) > 1)
      warning("\nMultiple probe items entered! Only first probe has been selected!")
    probe_ind <- probe[1]
  } else {
    stop("\nInvalid probe!\n")
    }
  
  
  exprs(expset) <- exprs(expset)[, rownames(expo)]
  
  # Required packages:
  require(mfp)
  require(stringr)
  require(Biobase)
  
  # Corrected select by correlations structure for a global 0.05(default):
  if( !is.null(select_adj) ){
    if(select_adj > 1 | select_adj < 0) {
      warning("Invalid select.adjust! It will be adjusted for a global 0.05 select value by the correlation matrix of expo.") 
      select_adj <- NULL
    }
    
  }
  
  if( is.null(select_adj) ){
    cormat <- cor(expo, use = "pairwise.complete.obs")
    M <- ncol(cormat)
    lambdas <- eigen(cormat)$values
    Vobs <- sum(((lambdas - 1)^2)) / (M - 1)
    Meff <- M * (1 - (M - 1) * Vobs / M^2) # number of effective tests
    select_adj <- select <- 1 - (1 - 0.05)^(1 / Meff) # corrected significance level (for global 0.05)
    rm(cormat, M, lambdas, Vobs, Meff)
  }
  
  gc(reset = TRUE)
  
  # Prepare data whith the corresponding outcome (probe):
  data <- data.frame(cbind(y = exprs(expset)[probe_ind[1],], expo))
  data <- data[complete.cases(data),]
  colnames(data) <- c("y", colnames(expo))
  
  # Fit the model:
  formula <- as.formula(paste0("y ~ ", paste0("fp(", colnames(expo), ", df = ", df, " )", collapse = " + ")))
  mod <- try(mfp(formula, select = select_adj, rescale = FALSE, verbose = FALSE, data = data), TRUE)
  mod$probe_name <- rownames(exprs(expset))[probe_ind[1]]
  mod$vars <- try(unique(as.character(str_sub(word(names(mod$coefficients)[-1],1), start = 1, end = -3))), TRUE)
  mod$p <- try(as.integer(sum(!is.na(mod$vars))), TRUE)
  
  class(mod) <- c("mfp.SingleFit", "mfp", "glm", "lm")
  mod

}


summary.mfp.SingleFit <- function(mod){
  # Get model performances:
  if(class(mod)[1] == "try-error"){
    stop("Error occured in SingleFit()!")
  } else {
    vars <- mod$vars
    p <- mod$p
    cor2 <- ifelse(p == 0, 0, try(cor(predict(mod), mod$y, use = "pairwise.complete.obs")^2, TRUE))
    
    # LRT:
    f_null <- glm(mod$y ~ 1)
    f_full <- mod$fit
    df.diff <- f_null$df.residual - f_full$df.residual
    vals <- (sum(residuals(f_null)^2) - sum(residuals(f_full)^2))/sum(residuals(f_full)^2) * f_full$df.residual 
    p_LRT <- pchisq(vals, df.diff, lower.tail = FALSE)
    rm(f_null, f_full, vals)
  }
  
  print(summary.glm(mod))
  
  cat(paste0("Likelihood Ratio Test p-value: ", format.pval(p_LRT) , "(df.diff=", df.diff ,")"), "\n")
  cat(paste0("Correlation between observed and predicted values: ", round(cor2, 4)), "\n")
  
}


plot.mfp.SingleFit <- function(mod, realpoints = FALSE, xlim, ylim, seed = 1234){
  
  if(class(mod)[1] == "try-error"){
    stop("Error occured in SingleFit()!")
  } else {
    
    if(mod$p == 0) stop("This model hasn't associated variables to be plotted!")

    vars <- mod$vars
    p <- mod$p
    X <- mod$X[,colnames(mod$X) %in% vars, drop = FALSE]
    
    cols <- colors()[-c(grep("gray", colors()), grep("white", colors()))]
    set.seed(seed)
    col <- sample(cols, p, rep = FALSE)
    
    if (missing(xlim))
    xlim <- apply(do.call(rbind, lapply(as.data.frame(X), 
                                        function(x) quantile(x, probs = c(.025, .975)))), 2, function(x) max(abs(x))) * c(-1, 1)

    if (missing(ylim))
     ylim <- quantile(mod$y, probs = c(.025, .975))
    
    
    plot(1, type = "n", ylab = expression(f[x]), xlab = "x", rug = FALSE, xlim = xlim, ylim = ylim,
         main = "Partial effects of x", sub = paste0(mod$probe_name), cex.main = 0.9, cex.sub = 0.75)
 
    faux <- function(i){
      n <- nrow(X)
      newX <- cbind(X[,i], apply(as.data.frame(X[,-i]), 2, function(x) rep(mean(x, na.rm = TRUE), n)))
      colnames(newX) <- c(vars[i], vars[-i])
      y_hat <- predict(mod, newdata = as.data.frame(newX))
      
      if(realpoints)
      {
        points(X[,i], predict(mod), col = adjustcolor(col[i], alpha = 0.35), pch = 20)
      }
      
      lines(X[,i][order(X[,i])], y_hat[order(X[,i])], col = col[i], lwd = 1.5)
    }
    
    invisible(sapply(1:p, faux))
    legend("topleft", vars, fill = col, cex = 1 / sqrt(p + 0.5), bty = "n")
  }
}


