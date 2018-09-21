#' filtering for an object of class nlAssoc
#'
#' @param list.meth a NAMED list of nlAssoc objects as obtained in fit.model(), usually applying different methods on the same data
#' @param filter if nlAssoc object has not been previously filtered it can be done here. Default is TRUE
#' @param param parameter of the table in results to filter with. Default is p
#' @param comp comparison param comp threshold, p below 0.05 default. Default is below
#' @param threshold numerical, the threshold to select associated variables related to the specified param. Default is 0.05
#'
#' @include filter.nlAssoc.R
#' @export compare.methods


compare.methods <- function(list.meth,filter=TRUE,param="p",comp="below",threshold=0.05){
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
      list.meth[[i]] <- filter.nlAssoc(list.meth[[i]], param = param, comp = comp, threshold = threshold)
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

