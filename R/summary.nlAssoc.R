#' Summary for an object of class nlAssoc
#'
#' @param res an nlAssoc object as obtained in fit.model()
#' @include getSignif.R
#' @include scoreVars.R

#' @method summary nlAssoc
#' @export

summary.nlAssoc <- function(res)
{

  #if filter show different output
  if (res$filter){
    cat(paste0("Filtered nlAssoc object using ", res$filter.param, " with threshold ",res$filter.threshold,"\n"))
    cat(paste0("Number of filtered probes: ", nrow(res$table),"\n"))
    cat(paste0("Number of assessed probes: ", nrow(res$object),"\n"))
    cat(paste0("Number of assessed variables: ", ncol(res$covars),"\n"))
  } else {
    cat("Object of class nlAssoc \n")
    cat(paste0("Number of assessed probes: ", nrow(res$object),"\n"))
    cat(paste0("Number of assessed variables: ", ncol(res$covars),"\n"))

    res.f<-getSignif(res,param="adj.p",threshold = 0.05)
    cat(paste0("Number of probes at an adjusted.p < 0.05: ", nrow(res.f$table), "\n"))

    res.f1<-getSignif(res,param="p",threshold = 0.01)
    cat(paste0("Number of probes at a p < 0.01: ", nrow(res.f1$table), "\n"))

    res.f2<-getSignif(res,param="p",threshold = 0.05)
    cat(paste0("Number of probes at a p < 0.05: ", nrow(res.f1$table), "\n"))

    res.score<-scoreVars(res)
    res.score.l <- length(res.score)
    if (res.score.l > 5) res.score <- res.score[1:5]
    cat(paste0("Top scoring variables : ", paste(names(res.score),collapse=","), "\n"))
    cat("\n")
    cat("Call: \n")
    print(res$call)

  }

}
