#' Function to select the best rf, but this could be in the fit.rforest function
#' adapted from varSelRF (varSelRF Ramon Diaz-Uriarte, only for classification)
#'
#' Random forest is obtained with all vars, oob error is computed and the most important variables
#' an iterated process consists on creating a new rf with the (1-vars.drop.frac)*vars more important obtained in the previous step
#' and generating a new rf model. The process stops when abs(mean(oob error)+sd(oob error)) is >= than previous step or only two vars are left
#'
#' @param xdata data set to be analyzed
#' @param Y var to analyze
#' @param ntree number of trees
#' @param vars.drop.frac fraction of variables removed at each iteration step, it selects the most important variables 1-fraction at each step


bestRF <- function(xdata, Y, ntree = 100, vars.drop.frac = 0.2,  verbose = FALSE)
{
  c.sd=1
  max.num.steps <- dim(xdata)[2]
  num.subjects <- dim(xdata)[1]
  keep.forest=T

  if(is.null(colnames(xdata))) colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")

  ## oversize the vectors; will prune later
  n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)

  oobError <- function(rf) {
   rsq <- rf$rsq
   rsq[rsq < 0] <- 0
   er.rsq<-1 - mean(rsq, na.rm = TRUE)
   return(er.rsq)
  }

    mtry <- floor(sqrt(ncol(xdata)))
    rf <- randomForest(x = xdata, y = Y, ntree = ntree, mtry = mtry, importance = TRUE, keep.forest = T)

  m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
  sd.iterated.ob.error <- sd.initial.ob.error <- sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) * (1/num.subjects))

  if(verbose) {
    print(paste("Initial 1 - pseudoR2: mean = ", round(m.initial.ob.error, 4), "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
  }

  importances <- importance(rf, type = 1, scale = FALSE)
  selected.vars <- order(importances, decreasing = TRUE)
  ordered.importances <- importances[selected.vars]

  initialImportances <- importances
  initialOrderedImportances <- ordered.importances

  j <- 1
  n.vars[j] <- dim(xdata)[2]
  vars[j] <- paste(colnames(xdata), collapse = " + ")
  OOB.rf[j] <- m.iterated.ob.error
  OOB.sd[j] <- sd.iterated.ob.error

  var.simplify <- TRUE

  while(var.simplify) {
    last.rf <- rf
    last.vars <- selected.vars
    previous.m.error <- m.iterated.ob.error
    previous.sd.error <- sd.iterated.ob.error

    if(length(selected.vars) <= 2) {
      var.simplify <- FALSE
      break
    }

    num.vars <- length(selected.vars)
    vars.drop <- round(num.vars * vars.drop.frac)

    if(num.vars >= (vars.drop + 2))
    {
      selected.vars <- selected.vars[1: (num.vars - vars.drop)]
      ordered.importances <- ordered.importances[1: (num.vars - vars.drop)]
    }
    else
    {
      selected.vars <- selected.vars[1:2]
      ordered.importances <- ordered.importances[1:2]
    }

    mtry <- floor(sqrt(length(selected.vars)))
    if(mtry > length(selected.vars)) mtry <- length(selected.vars)

    rf <- randomForest(x = xdata[, selected.vars], y = Y, importance = TRUE,  ntree = ntree, mtry = mtry, keep.forest = TRUE)
    m.iterated.ob.error <- oobError(rf)
    sd.iterated.ob.error <- sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) * (1/num.subjects))

    if(verbose) {
      print(paste("..... iteration ", j, "; 1 - R2: mean = ", round(m.iterated.ob.error, 4),
                  "; sd = ", round(sd.iterated.ob.error, 4), "; num. vars = ", length(selected.vars),
                  sep = ""))
    }
    j <- j + 1
    n.vars[j] <- length(selected.vars)
    vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error

  }

    n.vars <- n.vars[1:j]
    vars <- vars[1:j]
    OOB.rf<- OOB.rf[1:j]
    OOB.sd <- OOB.sd[1:j]
    min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
    best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= min.oob.ci)])]

    selected.vars <- sort(unlist(strsplit(vars[best.pos], " + ", fixed = TRUE)))
    selected.model = paste(selected.vars, collapse = " + ")

    return(selected.model)

}


