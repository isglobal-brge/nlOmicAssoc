#' filtering for an object of class nlAssoc, to get significative results
#'
#' @param res an nlAssoc object as obtained by nlOmicAssoc()
#' @param param parameter of the table to filter with (vars_n, aic, Cor2, p or adj.p). Default is "adj.p"
#' @param comp comparison param comp threshold. Default is "below"
#' @param threshold numerical, the threshold to select associated variables related to the specified param. Default is 0.05
#' @export getSignif


getSignif <- function (res,param="adj.p",comp="below",threshold=0.05, sort.results = TRUE)
{

  if (!param %in% c("vars_n", "aic", "Cor2", "p", "adj.p"))
    stop("param has to be one of the parameters of the table results: vars_n, aic, Cor2, p, adj.p")

  if (!comp %in% c("below", "above"))
    stop("comp has to be one of the followings: above, below")

    if (!(is.null(threshold))) {

    #in two steps in case there are errors
    #res$table[,param] <- as.numeric(res$table$p)
      if (comp == "below") {

        res.table.f <- res$table[!is.na(res$table[,param]) & res$table[,param] < threshold,]
        res.table.f.s <- res.table.f[order(res.table.f[,param]),]
        sel.bool <- res.table.f.s$probe

      } else {

        res.table.f <- res$table[!is.na(res$table[,param]) & res$table[,param] > threshold,]
        res.table.f.s <- res.table.f[order(res.table.f[,param], decreasing=T),]
        sel.bool <- res.table.f.s$probe

      }

      res$table <- res.table.f.s
      res$selected_vars <- res$selected_vars[sel.bool]
      res$final_formula <- res$final_formula[sel.bool]
      res$filter <- TRUE
      res$filter.param <- param
      res$filter.threshold <- threshold

      return(res)

   }
}
