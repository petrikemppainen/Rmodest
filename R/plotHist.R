#' Plots stacked histograms of posterior probabilities (Q-values)
#' 
#' This function uses information produced by \code{\link{getQdist}}, \code{\link{Qdist2KDAG}} and \code{\link{plotKDAG}}) and plots stacked histograms of the posterior probabilities for individuals belonging to a cluster
#' 
#' More info to come
#' 
#' @param K K at which histogram should be drawn
#' @param modest File produced by \code{\link{importData}}.
#' @param KDAGplot File produced by \code{\link{plotKDAG}}.
#' @param sort.by If \code{NA} (default) individuals appear in the same order as in the original data. If single integer, gives the Qcolumn by which the stacked histogram is ordered by. If a vector, gives the order of the Q-columns.
#' 
#' @return Returns a list with three objects. 
#' \itemize{
#'  \item \code{order}: the order for the Qcolumns which can be used to give all subsequent histograms the same order
#'  \item \code{Qcol_index}: the index of cloums of file \code{modest$Qcols} (produced by function \code{\link{importData}}) which represented each vertex in the KDAG (this information comes from file \code{KDAGplot} produced by function \code{\link{plotKDAG}}) 
#'  \item \code{hist.data} gives a matrix with the mean Q-value for each vertex (colum), where individuals (rows) are ordered as defined by \code{sort.by}
#'  }
#' @export
#' @keywords Qdist2KDAG, getQdist, importData, STRUCTURE, plotKDAG
#' @seealso \code{\link{importData}}, \code{\link{Qdist2KDAG}} and \code{\link{plotKDAG}}
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K)
#' KDAGplot <- plotKDAG(modest_MCMV, Qdist_MCMV, KDAG)
#' out <- plotHist(modest_MCMV, KDAGplot, 2, sort.by=1)
#' names(out)
#' out <- plotHist(modest_MCMV, KDAGplot, 4, sort.by=out$order)

plotHist <- function(modest, KDAGplot, K, sort.by=NA){
  
  ### get data ###
  vertices <- KDAGplot$vertices
  K_vertices <- KDAGplot$K_vertices
  color <- KDAGplot$color
  Qcols <- modest$Qcols
  
  ### use info from 'vertices' to calculate means ###
  temp <- do.call('cbind', lapply(vertices[K_vertices==K], function(x) if(length(x)>1) apply(Qcols[,x], 1, sum) else Qcols[,x]))
  
  # normalize
  hist.data <- apply(temp, 1, function(x) x/sum(x))
  
  # order
  if(!all(is.na(sort.by))){
    if(length(sort.by)!=1){
      new.order <- sort.by
    }else{new.order <- order(hist.data[sort.by,])}  
  }else{new.order <- c(1:ncol(hist.data))}
  
  # plot
  barplot(hist.data[,new.order], col=color[K_vertices==K], border=NA, space=0, axes=FALSE, ylab="Probability")
  
  # prepare output
  out <- list(new.order, vertices[K_vertices==K], t(hist.data))
  names(out) <- c("order", "Qcol_index", "hist.data")
  out
}


