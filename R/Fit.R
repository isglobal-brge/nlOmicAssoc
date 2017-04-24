
# ===========================================================================================
# Funcio per ajustar un model MFP a les N primeres expressions amb totes les exposicions
# ===========================================================================================

Fit <- function(expset, expo, df = 4, N = 1:nrow(exprs(expset)), select = 0.05, 
                select_adj = NULL, cores = 1L)
{
  # Required packages:
  require(mfp)
  require(stringr)
  require(Biobase)
  parallel_ind <- (cores != 1L | cores!= 1)
  if(parallel_ind) require(parallel)
  
  # Formula structure for the models:
  formula <- as.formula(paste0("y ~ ", paste0("fp(", colnames(expo), ", df = ", df, " )", collapse = " + ")))
  
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
  
  # Checkings:
  exprs(expset) <- exprs(expset)[, rownames(expo)]
  
  fmodel <- function(y, expo, probe)
  {
    gc(reset = TRUE)
    
    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, expo))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(expo))
    
    # Fit the model:
    mod <- try(mfp(formula, select = select_adj, rescale = FALSE, verbose = FALSE, data = data), TRUE)

    # Get model performances:
    if(class(mod)[1] == "try-error"){
      vars <- NA
      p <- NA
      cor2 <- NA
      p_LRT <- NA
    } else {
      vars <- try(unique(as.character(str_sub(word(names(mod$coefficients)[-1],1), start = 1, end = -3))), TRUE)
      p <- try(as.integer(sum(!is.na(vars))), TRUE)
      cor2 <- ifelse(p == 0, 0, try(cor(predict(mod), data$y, use = "pairwise.complete.obs")^2, TRUE))
      
      # LRT:
      f_null <- glm(data$y ~ 1)
      f_full <- mod$fit
      df.diff <- f_null$df.residual - f_full$df.residual
      vals <- (sum(residuals(f_null)^2) - sum(residuals(f_full)^2))/sum(residuals(f_full)^2) * f_full$df.residual 
      p_LRT <- pchisq(vals, df.diff, lower.tail = FALSE)

      rm(f_null, f_full, df.diff, vals)
    }
    
    
    gc(reset = TRUE)
    
    # Data frame of the results:
    taula <- try(data.frame(cbind(probe = probe, p = p, Cor2 = cor2, LRT_pval = p_LRT)), TRUE)
    
    result <- try(list(table = taula, selected_vars = vars), TRUE)
    result
  }
  
  if(parallel_ind){
    # Usar cores amb mclapply!!!!!!!!!!!!!!!!!!!!!!
    results <- try(mcmapply(function(y, ny) fmodel(y = y, expo = expo, probe = ny), 
                          y = apply(exprs(expset)[N,], 1, list), 
                          ny = rownames(exprs(expset))[N], 
                          SIMPLIFY = FALSE, mc.cores = cores), TRUE)
  } else {
    results <- try(mapply(function(y, ny) fmodel(y = y, expo = expo, probe = ny), 
                          y = apply(exprs(expset)[N,], 1, list), 
                          ny = rownames(exprs(expset))[N], 
                          SIMPLIFY = FALSE), TRUE)
  }
  
  
  faux <- function(x){
    taula <- x[[1]]
    if(class(taula) == "try-error") taula <- data.frame(cbind(probe = taula$probe, p = NA, Cor2 = NA, LRT_pval = NA))
    taula
  }
  
  taula <- as.data.frame(do.call(rbind, lapply(results, faux)))
  
  p_adj <- p.adjust(as.double(as.character(taula$LRT_pval)), method = "BH")
  taula$adj_pval <- format.pval(p_adj)
  taula$LRT_pval <- format.pval(as.double(as.character(taula$LRT_pval)))
  taula <- taula[ aux <- base::order(as.double(taula$LRT_pval), decreasing = FALSE),]
  rownames(taula) <- NULL
  selected_vars <-  lapply(results, function(x) x[[2]])[aux]
  
  
  # taula
  list(table = taula, slected_vars = selected_vars, alpha_select = select_adj)
}


# ===============================================================================

