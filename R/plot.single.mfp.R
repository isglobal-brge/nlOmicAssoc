#' plot for an object of class nlAssoc, obtatined with model mfp, for just a feature
#'
#' @param res an nlAssoc object as obtained by nlOmicAssoc
#' @param probes character, id of the probe to plot,
#' @param xlim xlim for plot
#' @param ylim ylim for plot
#' @param seed seed for internal random processes
#' @export plot.single.mfp

plot.single.mfp <- function(res,probe, xlim, ylim, seed = 1234){

  if (!("nlAssoc" %in% class(res) )) stop("res must be an nlAssoc objcect")

  model <- res$call$model

  if (model!="mfp") stop("This plot is only for results obtained with the mfp method")

  sel.vars <-res$selected_vars[[probe]]$selected_vars

  if(is.na(sel.vars[1])) stop("This model hasn't associated variables to be plotted!")
  covars <- res$covars
  formul<-res$final_formula[[probe]]$final_formula
  sel.vars <-res$selected_vars[[probe]]$selected_vars
  X <- covars[,colnames(covars) %in% sel.vars, drop = FALSE]
  X <- scale(X, center=F)

  cols <- colors()[-c(grep("gray", colors()), grep("white", colors()))]
  set.seed(seed)
  col <- sample(cols, length(sel.vars), rep = FALSE)

  faux <- function(i){
    n <- nrow(X)
    newX <- cbind(X[,i], apply(as.data.frame(X[,-i]), 2, function(x) rep(mean(x, na.rm = TRUE), n))) #deixa la variable en estudi i fixa la mitjana
    colnames(newX) <- c(sel.vars[i], sel.vars[-i])
    data <- data.frame(cbind(y = res$object[probe,], covars))
    data <- data[complete.cases(data),]
    mod <- mfp(formul, data = data)
    y_hat <- predict(mod, newdata = as.data.frame(newX))
  }
  Y <- invisible(sapply(1:length(sel.vars), faux))

  if (missing(xlim))
    xlim <- range(as.vector(X))
  if (missing(ylim))
    ylim <- quantile(as.vector(Y), probs = c(.1, .9), na.rm = T)
  par(mar=c(4,3,3,1))
  plot(1, type = "n", ylab = expression(f[x]), xlab = "x", xlim = xlim, ylim = ylim, #review params rug
       main = paste0("Partial effects of ",probe), cex.main = 0.7, cex.axis = 0.7, mgp=c(2,1,0))

  faux2 <- function(i){
    print.x <- X[,i][order(X[,i])]
    print.y <- Y[,i][order(X[,i])]
    cc.log <- complete.cases(print.x, print.y)
    lines(x=print.x[cc.log], y=print.y[cc.log], col = col[i], lwd = 1.5)
  }

  invisible(sapply(1:length(sel.vars), faux2))
  legend("topright", sel.vars, pch=20, col=col, cex = 0.4)

}


