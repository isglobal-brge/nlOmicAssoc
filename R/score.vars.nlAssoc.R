#' scoring the variables of selected_vars slot in nlAssoc object
#'
#' @param res an nlAssoc object as obtained in fit.model()
#' @export score.vars.nlAssoc

score.vars.nlAssoc <- function (res)
{
  #we will return a score of the repeated vars/number of tested probes
  #if there are many null results, % will decrease
  if (!(is.null(res))) {
    ll<-length(res$selected_vars)
    vars_all <- sort(table(unlist(res$selected_vars)),decreasing=T)/ll

    return(vars_all)
   }
}
