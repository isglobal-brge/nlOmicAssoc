#' Compare nlAssoc objects obtained by nlOmicAssoc()
#'
#' @param list.meth a NAMED list of nlAssoc objects as obtained by nlOmicAssoc(), usually applying different methods on the same data
#' @param filter if nlAssoc object has not been previously filtered it can be done here. Default is TRUE
#' @param param parameter of the table to filter with (vars_n, aic, Cor2, p or adj.p). Default is "adj.p"
#' @param comp comparison param comp threshold. Default is "below"
#' @param threshold numerical, the threshold to select associated variables related to the specified param. Default is 0.05
#'
#' @include getSignif.R
#' @export compareMethods


compareMethods <- function(list.meth,filter=TRUE,param="adj.p",comp="below",threshold=0.05){
  #names <- deparse(substitute(list.meth))
  # els = substitute(list.meth)
  # names=sapply(els, deparse)
  #how to get object names if names is empty?)
  nam = names(list.meth)
  #puc fer un score de quantes vegades surten les variables sobre els probes seleccionats
  list.meth.l=length(list.meth)
  if (filter){
    for (i in 1:list.meth.l){
      #replace current object in list
      list.meth[[i]] <- getSignif(list.meth[[i]], param = param, comp = comp, threshold = threshold)
    }
  }

  int.m<-matrix(NA,nrow=list.meth.l,ncol=list.meth.l) #diag sup, voldria inferior

  if (is.list(list.meth)){
    for (i in 1:list.meth.l){
      for(j in 1:list.meth.l){
        if (j>=i)  int.m[i,j]<-length(intersect(list.meth[[i]]$table$probe,
                                               list.meth[[j]]$table$probe))
      }
    }

  }

  rownames(int.m) <- nam
  colnames(int.m) <- nam

  return(int.m)
}

