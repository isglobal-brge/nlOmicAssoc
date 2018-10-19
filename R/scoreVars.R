#' scoring the variables of selected_vars slot in an nlAssoc object
#'
#' @param res an nlAssoc object as obtained by nlOmicAssoc()
#'
#' @export scoreVars

scoreVars <- function (res)
{
  #we will return a score of the repeated vars/number of tested probes
  #if there are many null results, % will decrease
  if (!(is.null(res))) {
    ll<-length(res$selected_vars)
    vars_all <- sort(table(unlist(res$selected_vars)),decreasing=T)/ll

    return(vars_all)
   }
}
