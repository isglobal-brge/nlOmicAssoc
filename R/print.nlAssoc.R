#' Print for an object of class nlAssoc showing only first registries
#'
#' @param res an nlAssoc object as obtained in fit.model()

#' @method print nlAssoc
#' @export

print.nlAssoc <- function(res)
{

  #if filter different output
  if (res$filter){
    cat(paste0("Filtered object of class nlAssoc using ", res$filter.param, " with threshold ",res$filter.threshold,"\n"))
    cat(paste0("Number of filtered probes: ", nrow(res$table),"\n"))
    cat(paste0("Number of assessed probes: ", nrow(res$object),"\n"))
    cat(paste0("Number of assessed variables: ", ncol(res$covars),"\n"))
  } else {
    cat("Object of class nlAssoc \n")
    cat(paste0("Number of assessed probes: ", nrow(res$object),"\n"))
    cat(paste0("Number of assessed variables: ", ncol(res$covars),"\n"))
    cat("\n")
    cat("Call: \n")
    print(res$call)
  }

}
