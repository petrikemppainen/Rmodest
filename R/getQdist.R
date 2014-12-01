#' Calculates distances between all pairwise Qcolumns
#'
#' Based on information produced by \code{\link{importData}}, \code{\link{getQdist}} calculates distances between all pairwise Q-columns
#' 
#' More info to come.
#'
#' @param modest A file produced by \code{importData} that contains parsed data from replicate STRUCTURE runs
#' @param mc.cores The number of cores to use, i.e. how many processes will be spawned (at most). By default all available cores are used.
#' @keywords STRUCTURE, importData, getQdist, 
#' @seealso \code{\link{importData}} and \code{\link{Qdist2KDAG}}
#' @return An upper diagonal matrix with pairwise distances
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @export
#' @importFrom parallel mclapply
#' @importFrom Deducer likelihood.test
#' @examples
#' \dontrun{
#' data(MODEST)
#' Qdist <- getQdist(modest_MCMV)}

getQdist <- function(modest, mc.cores=getOption("mc.cores")){
  infile <- modest$allele_counts
  pair <- combn(1:ncol(infile), 2, simplify = F) 
  
  G_test <- function(x){
    likelihood.test(infile[, c(x[1],x[2])])$statistic
  }
  out <- mclapply(pair, G_test, mc.cores=mc.cores)  
  corm <- cbind(do.call(rbind, pair), unlist(out))
  g=graph.data.frame(as.data.frame(corm))
  E(g)$weight <- corm[,3]
  Qdist <- get.adjacency(g,sparse=FALSE, attr="weight")
  Qdist <- (Qdist - min(na.omit(Qdist)))/diff(range(na.omit(Qdist)))    
  Qdist[lower.tri(Qdist)] <- NA
  diag(Qdist) <- NA
  Qdist
}
