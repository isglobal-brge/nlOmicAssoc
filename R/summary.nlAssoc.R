#' Summary for an object of class nlAssoc
#'
#' @param res an nlAssoc object as obtained in fit.model()
#' @include filter.nlAssoc.R
#' @include score.vars.nlAssoc.R
#' @export summary.nlAssoc

summary.nlAssoc <- function(res)
{
  cat(paste0("Call: ",  res$call), "\n")
  cat(paste0("Number of assessed probes: ", nrow(res$table)),"\n")

  res.f<-filter.nlAssoc(res,param="adj.p",threshold = 0.05)
  cat(paste0("Number of probes at an adjusted.p < 0.05: ", nrow(res.f$table)), "\n")

  res.f1<-filter.nlAssoc(res,param="p",threshold = 0.01)
  cat(paste0("Number of probes at a p < 0.01: ", nrow(res.f1$table)), "\n")

  res.f2<-filter.nlAssoc(res,param="p",threshold = 0.05)
  cat(paste0("Number of probes at a p < 0.05: ", nrow(res.f1$table)), "\n")

  res.score<-score.vars.nlAssoc(res)
  res.score.l <- length(res.score)
  if (res.score.l > 5) res.score <- res.score[1:5]
  cat(paste0("Top scoring variables : ", paste(names(res.score),collapse=",")), "\n")

}
