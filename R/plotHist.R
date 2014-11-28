#' Plots stacked histograms posterior probabilities for beloning to a cluster taking into account multi-modality and integrating over multiple Ks.
#' 
#' This function uses information produced by \code{\link{getQdist}}, \code{\link{Qdist2KDAG}} and \code{\link{plotKDAG}}) and plots stacked histograms of the posterior probabilities for individuals belonging to a cluster (integrated over replicated STRUCTURE runs with different Ks).
#' 
#' More info to come
#' 
#' @param K K at which histogram should be drawn
#' @param modest File produced by \code{\link{importData}}.
#' @param KDAGplot File produced by \code{\link{plotKDAG}}.
#' @param sort.by If single integer, gives the Qcolumn by which the stacked histogram is ordered by. If a vector, gives the order of the Q-columns.
#' @return Returns the order for the Qcolumns such that all other histograms can be given the same order
#' @export
#' @keywords Qdist2KDAG, getQdist, importData, STRUCTURE, plotKDAG
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K, merge.by="Qdist")
#' KDAGplot <- plotKDAG(modest_MCMV, Qdist_MCMV, KDAG)
#' plotHist(modest_MCMV, KDAGplot, 4, sort.by=1)

plotHist <- function(modest, KDAGplot, K, sort.by=1){
  vertices <- KDAGplot$vertices
  K_vertices <- KDAGplot$K_vertices
  color <- KDAGplot$color
  Qcols <- modest$Qcols
  
  temp <- do.call('cbind', lapply(vertices[K_vertices==K], function(x) if(length(x)>1) apply(Qcols[,x], 1, sum) else Qcols[,x]))
  hist.data <- apply(temp, 1, function(x) x/sum(x))
  if(length(sort.by)!=1){
    new.order <- sort.by
  }else{new.order <- order(hist.data[sort.by,])}
  barplot(hist.data[,new.order], col=color[K_vertices==K], border=NA, space=0, axes=FALSE, ylab="Probability")
  new.order
}


