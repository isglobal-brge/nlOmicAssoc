#' Function to select the best nnet, but this could be in the fit.nnetwork function
#'
#' nnet is obtained with all vars, so-called oob error
#' 1 - cor(predict(Nnet)[complete.cases(predict(Nnet)),1], Y[complete.cases(cbind(Y, xdata[, Nnet$coefnames]))])^2)
#  is computed and the most important variables using varImp(caret) are taken
#' an iterated process consists on creating a new nnet with the (1-vars.drop.frac)*vars more important obtained in the previous step
#' and generating a new nnet model. The process stops when abs(mean(oob error)+sd(oob error)) is >= than previous step or only two vars are left

#' @param xdata data set to be analyzed
#' @param Y var to analyze
#' @param size nnet param: number of units in the hidden layer. Can be zero if there are skip-layer units
#' @param vars.drop.frac fraction of variables removed at each iteration step, it selects the most important variables 1-fraction at each step


bestNN <- function(xdata, Y, size = 2, vars.drop.frac = 0.2, verbose = FALSE)
{

  c.sd = 1
  max.num.steps <- dim(xdata)[2]
  num.subjects <- dim(xdata)[1]

  if(is.null(colnames(xdata))) colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep ="")

  n.vars <- vars <- v.nnet <- v.sd <- rep(NA, max.num.steps)

 vError <- function(Nnet) {
   v.err=1 - cor(predict(Nnet)[complete.cases(predict(Nnet)),1], Y[complete.cases(cbind(Y, xdata[, Nnet$coefnames]))])^2
   return(v.err)
  }

  data.aux <- cbind(Y, xdata)
  formula.aux <- as.formula(paste0(colnames(data.aux)[1], " ~ ", paste0(colnames(data.aux)[-1], collapse = " + ")))
  Nnet <- nnet(formula.aux, data = data.aux, size = size, linout = TRUE, trace = FALSE)

  m.iterated.ob.error <- m.initial.ob.error <- vError(Nnet)
  sd.iterated.ob.error <- sd.initial.ob.error <- sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) * (1/num.subjects))

  if(verbose) {
     print(paste("Initial SE = 1 - cor2 = ", round(m.initial.ob.error, 4), "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
  }

  importances<-abs(varImp(Nnet))[[1]]

  selected.vars <- order(importances, decreasing = TRUE)
  ordered.importances <- importances[selected.vars]

  initialImportances <- importances
  initialOrderedImportances <- ordered.importances

  j <- 1
  n.vars[j] <- dim(xdata)[2]
  vars[j] <- paste(colnames(xdata), collapse = " + ")
  v.nnet[j] <- m.iterated.ob.error
  v.sd[j] <- sd.iterated.ob.error

  var.simplify <- TRUE

  while(var.simplify) {

    last.nnet <- Nnet
    last.vars <- selected.vars
    previous.m.error <- m.iterated.ob.error
    previous.sd.error <- sd.iterated.ob.error

    if(length(selected.vars) <= 2) {
      var.simplify <- FALSE
      break
    }

    num.vars <- length(selected.vars)
    vars.drop <- round(num.vars * vars.drop.frac)

    if(num.vars > (vars.drop) + 1) {
      selected.vars <- selected.vars[1: (num.vars - vars.drop)]
      ordered.importances <- ordered.importances[1: (num.vars - vars.drop)]
    }
    else {
      selected.vars <- selected.vars[1]
      ordered.importances <- ordered.importances[1]
    }

    data.aux <- cbind(Y, xdata[, selected.vars])
    formula.aux <- as.formula(paste0(colnames(data.aux)[1], " ~ ", paste0(colnames(data.aux)[-1], collapse = " + ")))
    Nnet <- invisible(nnet(formula.aux, data = data.aux, size = size, linout = TRUE, trace = FALSE))

    m.iterated.ob.error <- vError(Nnet)
    sd.iterated.ob.error <- sqrt(m.iterated.ob.error * (1 - m.iterated.ob.error) * (1/num.subjects))

    if(verbose) {
      print(paste("..... iteration ", j, "; SE = 1 - cor2 = ", round(m.iterated.ob.error, 4),
            "; sd = ", round(sd.iterated.ob.error, 4), "; num. vars = ", length(selected.vars), sep = ""))
    }
    j <- j + 1

    n.vars[j] <- length(selected.vars)
    vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
    v.nnet[j] <- m.iterated.ob.error
    v.sd[j] <- sd.iterated.ob.error

  }

    n.vars <- n.vars[1:j]
    vars <- vars[1:j]
    v.nnet <- v.nnet[1:j]
    v.sd <- v.sd[1:j]
    min.oob.ci <- min(v.nnet,na.rm=T) + c.sd * v.sd[which.min(v.nnet)]
    best.pos <- which(v.nnet <= min.oob.ci)[which.min(n.vars[which(v.nnet <= min.oob.ci)])]

    selected.vars <- sort(unlist(strsplit(vars[best.pos], " + ", fixed = TRUE)))
    selected.model = paste(selected.vars, collapse = " + ")
    return(selected.model)

}

