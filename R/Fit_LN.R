
# ===========================================================================================
# Funcio per ajustar un model MFP a de tots els valors de la matriu data_m (exprs(ExpressionSet) per exemple) amb totes les variables (exposicions)
# ===========================================================================================
#Elimino la N, si volem fer amb algunes que facin el subsetting

load("./nlOmicAssoc/gexp4b.Rdata")
data_m_all<-exprs(gexp4b)
exps<-read.table("./nlOmicAssoc/exposures.tsv",sep = "\t",dec=".",stringsAsFactors = F,header = T)

sel.ids<-colnames(data_m_all)[1:5]
sel.ids<-colnames(data_m_all)[c(1,2,4,5)]
# 
data_m<-data_m_all[1:100,sel.ids]
vars_df<-exps[exps$SubjectID %in% sel.ids,2:5] #3 is NA for all
rownames(vars_df) <- exps[exps$SubjectID %in% sel.ids,1]
fit.all<-Fit(data_m,vars_df,cores=4)

#provo d'aplicar-ho tot 

data_m<-data_m_all[1:1000,]
vars_df<-exps[,-1] #trec la primera i li poso a l'id
rownames(vars_df) <- exps[,1]
#vars_df<-vars_df[,1:20]
# 
# #hem de filtrar variables amb més del 70% de NAs o imputar-les
# 
# library(R.oo)
# equals(rownames(vars_df),colnames(data_m_all)) #TRUE
# 
# cor(vars_df,use="pairwise.complete.obs") #la variable PFBS és tot NA, l'eliminem, però hem de controlar-ho també a l'entrar
# 
Sys.time()
fit.all<-Fit(data_m,vars_df,cores=4)
Sys.time()

complete.cases(vars_df) #ni un!
cor(vars_df)

Fit <- function(data_m, vars_df, df = 4, select = 0.05, select_adj = NULL, cores = 1L)
{
  # Required packages:
  require(mfp)
  require(stringr)
  require(Biobase)
  parallel_ind <- (cores != 1L | cores!= 1)
  if(parallel_ind) require(parallel)
  
  
  # Corrected select by correlations structure for a global 0.05(default):
  if( !is.null(select_adj) ){
    if(select_adj > 1 | select_adj < 0) {
      warning("Invalid select.adjust! It will be adjusted for a global 0.05 select value by the correlation matrix of vars_df") 
      select_adj <- NULL
    }
    
  }
  
  # Check vars structure
  vars.na <- apply(is.na(vars_df),2,sum)
  # just select those variables with < 70% NA
  vars_df <- vars_df[,(vars.na/ncol(vars_df)) <= 0.6]
  cat(paste(ncol(vars_df),"were selected for the analysis"))
  
  # Formula structure for the models:
  formula <- as.formula(paste0("y ~ ", paste0("fp(", colnames(vars_df), ", df = ", df, " )", collapse = " + ")))
  
  if( is.null(select_adj) ){
    cormat <- cor(vars_df, use = "pairwise.complete.obs")
    #if there is an na in cormat stop
    if (anyNA(cormat)) stop("Some NAs in correlation matrix, try to impute some of the missing values")
    M <- ncol(cormat)
    lambdas <- eigen(cormat,only.values = T)$values
    Vobs <- sum(((lambdas - 1)^2)) / (M - 1)
    Meff <- M * (1 - (M - 1) * Vobs / M^2) # number of effective tests
    select_adj <- select <- 1 - (1 - 0.05)^(1 / Meff) # corrected significance level (for global 0.05)
    rm(cormat, M, lambdas, Vobs, Meff)
  }
  
  gc(reset = TRUE)
  
  # Checkings:
  data_m <- data_m[, rownames(vars_df)]
  
  fmodel <- function(y, vars_df, probe)
  {
    gc(reset = TRUE)
    
    # Prepare data whith the corresponding outcome:
    data <- data.frame(cbind(y = y, vars_df))
    data <- data[complete.cases(data),]
    colnames(data) <- c("y", colnames(vars_df))
    
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
    results <- try(mcmapply(function(y, ny) fmodel(y = y, vars_df = vars_df, probe = ny), 
                          y = apply(data_m, 1, list), 
                          ny = rownames(data_m)[], 
                          SIMPLIFY = FALSE, mc.cores = cores), TRUE)
  } else {
    results <- try(mapply(function(y, ny) fmodel(y = y, vars_df = vars_df, probe = ny), 
                          y = apply(data_m, 1, list), 
                          ny = rownames(data_m), 
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
  list(table = taula, selected_vars = selected_vars, analysis_vars= colnames(vars_df), alpha_select = select_adj)
}


# ===============================================================================

