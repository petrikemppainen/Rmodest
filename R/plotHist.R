#' Plots stacked histograms posterior probabilities (Q-values)
#' 
#' This function uses information produced by \code{\link{getQdist}}, \code{\link{Qdist2KDAG}} and \code{\link{plotKDAG}}) and plots stacked histograms of the posterior probabilities for individuals belonging to a cluster (integrated over replicated STRUCTURE runs with different Ks).
#' 
#' More info to come
#' 
#' @param K K at which histogram should be drawn
#' @param modest File produced by \code{\link{importData}}.
#' @param KDAGplot File produced by \code{\link{plotKDAG}}.
#' @param sort.by If \code{NA} (default) individuals appear in the same order as in the original data. If single integer, gives the Qcolumn by which the stacked histogram is ordered by. If a vector, gives the order of the Q-columns.
#' @return Returns a list with three objects. \code{order}: the order for the Qcolumns which can be used to give all subsequent histograms the same order. \code{Qcol_index}: the index of cloums of file \code{modest$Qcols} which represented each vertex in the KDAG (this information comes from file \code{KDAGplot} produced by function \code{\link{plotKDAG}}). \code{hist.data} gives the mean Q-value for each vertex, where individuals are ordered as defined by \code{sort.by}.
#' @export
#' @keywords Qdist2KDAG, getQdist, importData, STRUCTURE, plotKDAG
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K, merge.by="Qdist")
#' KDAGplot <- plotKDAG(modest_MCMV, Qdist_MCMV, KDAG)
#' out <- plotHist(modest_MCMV, KDAGplot, 2, sort.by=1)
#' names(out)
#' out <- plotHist(modest_MCMV, KDAGplot, 4, sort.by=out$order)

plotHist <- function(modest, KDAGplot, K, sort.by=NA){
  vertices <- KDAGplot$vertices
  K_vertices <- KDAGplot$K_vertices
  color <- KDAGplot$color
  Qcols <- modest$Qcols
  
  temp <- do.call('cbind', lapply(vertices[K_vertices==K], function(x) if(length(x)>1) apply(Qcols[,x], 1, sum) else Qcols[,x]))
  hist.data <- apply(temp, 1, function(x) x/sum(x))
  if(!all(is.na(sort.by))){
    if(length(sort.by)!=1){
      new.order <- sort.by
    }else{new.order <- order(hist.data[sort.by,])}  
  }
  barplot(hist.data[,new.order], col=color[K_vertices==K], border=NA, space=0, axes=FALSE, ylab="Probability")
  hist.data
  out <- list(new.order, vertices[K_vertices==K], hist.data)
  names(out) <- c("order", "Qcol_index", "hist.data")
  out
}


