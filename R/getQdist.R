#' Estimate cluster distance (Qdist)
#'
#' Based on information produced by \code{\link{importData}}, \code{\link{getQdist}} calculates cluster distances between all pairwise Q-columns
#' 
#' We use a likelihood ratio test of independence (G-test: \code{Deducer::likelihood.test}) of allele counts between clusters (contribution from each locus is weighted by the proportion of missing). Based on this distance, \code{Qdist}, a single linkage clustering tree is generated (see \code{\link{hclust}}), method="single" for details) which is then parsed into  a K-directed acyclic graph (KDAG) by the function: \code{\link{Qdist2KDAG}}. Also function \code{\link{plotKDAG}} requires \code{Qdist} as an input.
#'
#' @param modest A file produced by \code{importData} that contains parsed data from STRUCTURE output data
#' @param mc.cores The number of cores to use, i.e. how many processes will be spawned (at most). Default is one (1). Using more than one core in GUI or embedded environmentsis is not recommended.
#' @keywords STRUCTURE, importData, getQdist, 
#' @seealso \code{\link{importData}}, \code{\link{Qdist2KDAG}} and \code{\link{plotKDAG}}
#' @return An upper diagonal matrix with pairwise Q-distances (\code{Qdist})
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @export
#' @importFrom parallel mclapply
#' @importFrom Deducer likelihood.test
#' @examples
#' \dontrun{
#' data(MODEST)
#' Qdist <- getQdist(modest_MCMV)}

getQdist <- function(modest, mc.cores=1){
  ### prepare data ###
  infile <- modest$allele_counts
  pair <- combn(1:ncol(infile), 2, simplify = F) 
  
  ### the test ###
  G_test <- function(x){
    likelihood.test(infile[, c(x[1],x[2])])$statistic
  }
  out <- mclapply(pair, G_test, mc.cores=mc.cores)  
  corm <- cbind(do.call(rbind, pair), unlist(out))
  
  ### parse data to a upper diagonal matrix ###
  g=graph.data.frame(as.data.frame(corm))
  E(g)$weight <- corm[,3]
  Qdist <- get.adjacency(g,sparse=FALSE, attr="weight")
  Qdist <- (Qdist - min(na.omit(Qdist)))/diff(range(na.omit(Qdist)))    
  Qdist[lower.tri(Qdist)] <- NA
  diag(Qdist) <- NA
  Qdist
}
