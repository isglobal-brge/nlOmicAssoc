#' plot for an object of class nlAssoc
#'
#' @param res an nlAssoc object as obtained by nlOmicAssoc
#' @param top.probes numeric, number of probes to plot,
#' @param pdf.name name and path of pdf file to store plots. Default NULL

#' @import ggplot2
#' @import stats
#' @import Biobase
#' @import parallel
#' @import mice
#' @import stringr
#' @import caret
#' @import splines
#' @import mfp
#' @import partDSA
#' @import gam
#' @import mboost
#' @import randomForest
#' @import nnet
#' @import NeuralNetTools
#'
#' @include getSignif.R
#'
#' @method plot nlAssoc
#' @export


plot.nlAssoc <- function(res, top.probes=10, pdf.name=NULL)
{
  #for the top probes (), plot predicted vs original, keep original objects??? noo, keep imputed!

  #par computing
  cat("Warning: Plot will be generated reconstructing the models with the results, this might take a while","\n")

  #if results are not writen in pdf just 12
  if (is.null(pdf.name)){
      if (top.probes>20) stop ("Interactive plot only available maximal for 20 results. Please, specify a pdf.name, filter or subset before plotting")
  }

  model <- res$call$model
  res.f <- getSignif(res, threshold = 0.9999) #very relaxing threshold to get something
  res.f
  sorted.probes<-res.f$table$probe[1:top.probes]
  sorted.vars<-res.f$selected_vars[1:top.probes]

  data_m <- res$object
  vars_df <- res$covars

  if (!is.null(pdf.name)) pdf(file=pdf.name)

  if (model=="ewas"){
    #pl <- list()
    for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      formula <- as.formula(paste0("y"," ~ ", paste0("ns(", vars, ", df = ",3,")", collapse = " + ")))
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- lm(formula, data = data)
      mod.pred <- predict(mod)
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  } else if (model=="dsa"){

     for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- partDSA(x = data[,-1], y = data[,1], control = DSA.control(vfold = 2, cut.off.growth = 15, minsplit = 20, MPD = 0.2, save.input = TRUE))
      mod.pred <- mod$pred.test.set.DSA[,ncol(mod$pred.test.set.DSA)]
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  } else if (model=="mfp"){

     for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      formula <-  formula <- as.formula(paste0("y ~ ", paste0("fp(", vars, ", df = ", 4, " )", collapse = " + ")))
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- lm(formula, data = data)
      mod.pred <- predict(mod)
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  } else if (model=="gam"){

     for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      formula <- as.formula(paste0("y ~ ", paste0("s(", vars, ", df = ", 4, " )", collapse = " + ")))
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- gam(formula, data = data)
      mod.pred <- predict(mod)
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  } else if (model=="gamboost"){

     for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      formula <- as.formula(paste0("y ~ ", paste0(vars, collapse = " + ")))
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- gamboost(formula, data = data)
      mod.pred <- predict(mod)
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  } else if (model=="rforest"){

    for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      formula <- as.formula(paste0("y ~ ", paste0(vars, collapse = " + ")))
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- randomForest(formula, data = data, ntree = 2000)
      mod.pred <- mod$predicted
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  } else if (model=="nnetwork"){

     for (i in 1:top.probes){
      probe <- sorted.probes[i]
      vars <- sorted.vars[probe][[1]]$selected_vars
      formula <- as.formula(paste0("y ~ ", paste0(vars, collapse = " + ")))
      data <- data.frame(cbind(y = data_m[probe,], vars_df))
      data <- data[complete.cases(data),]
      mod <- nnet(formula, data = data, size = 2, decay = 0.001, maxit = 1000, linout = TRUE, trace = FALSE)
      mod.pred <-predict(mod,type="raw")
      p <- qplot(x=data_m[probe,rownames(data)],y=mod.pred, xlab = paste("probe",probe), ylab = "predicted")
      p <- p + geom_smooth(method = "loess") #no span
      p <- p + theme_bw()
      print(p)
    }

  }

  if (!is.null(pdf.name)) dev.off()

}
