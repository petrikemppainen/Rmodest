#' Calculates distances between all pairwise Qcolumns
#'
#' Based on information produced by \code{\link{importData}}, \code{\link{getQdist}} calculates distances between all pairwise Q-columns
#' 
#' The method \code{Pearson} is the fastest but is only based on the Q-values instead of the allele counts for all loci in the data set and is thereore likely to be missleading. Use this option only for quick testing. More info to come.
#'
#' @param modest A file produced by \code{importData} that contains parsed data from replicate STRUCTURE runs
#' @param method Distance measure, either based on an explicit mode (\code{explicit_model}), the likelihood ratio test (\code{G_test}) or the Pearson distance {\code{Pearson}}
#' @param scaled By default the distance measure is scaled between 0 and 1, but can be turned off (\code{FALSE})
#' @param mc.cores The number of cores to use, i.e. how many processes will be spawned (at most). By default all available cores are used.
#' @keywords STRUCTURE, importData, getQdist, 
#' @seealso \code{\link{importData}} and \code{\link{Qdist2KDAG}}
#' @return An upper diagonal matrix with pairwise distances
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @export
#' @examples
#' data(MODEST)
#' Qdist <- getQdist(modest_MCMV, "Pearson")

getQdist <- function(modest, method="G_test", scaled=TRUE, mc.cores=getOption("cores")){
  G_test <- function(x){
    likelihood.test(infile[, c(x[1],x[2])])$statistic
  }
  Explicit_model <- function(x){
    sum(na.omit(sapply(as.numeric(unique(rownames(infile))), function(y){
      temp <- infile[rownames(infile)==y, c(x[1],x[2])]+1
      if(!is.vector(temp)){
        pA <- temp[,1]/sum(temp[,1])
        pB <- temp[,2]/sum(temp[,2])
        nA <- temp[,1]
        nB <- temp[,2]
        suppressWarnings(DirichletDistance(pA, pB, nA, nB))
      }else{0}    
    })))
  }
  if(method=="G_test"){
    infile <- modest$allele_counts
    pair <- combn(1:ncol(infile), 2, simplify = F) 
    out <- mclapply(pair, G_test, mc.cores=mc.cores)
  }
  if(method=="Explicit_model"){
    infile <- modest$allele_counts
    pair <- combn(1:ncol(infile), 2, simplify = F) 
    out <- mclapply(pair, Explicit_model, mc.cores=mc.cores)
  }
  if(method=="Pearson"){
    Qdist <- 1-cor(modest$Qcols)
  }else{
    corm <- cbind(do.call(rbind, pair), unlist(out))
    g=graph.data.frame(as.data.frame(corm))
    E(g)$weight <- corm[,3]
    Qdist <- get.adjacency(g,sparse=FALSE, attr="weight")
  }
  if(scaled==TRUE){
    Qdist <- (Qdist - min(na.omit(Qdist)))/diff(range(na.omit(Qdist)))    
  }
  Qdist[lower.tri(Qdist)] <- NA
  diag(Qdist) <- NA
  Qdist
}

DirichletDistance <- function(pA, pB, nA, nB){
  pAB <- (pA*nA + pB*nB)/(nB + nA)
  pAB <- pAB/sum(pAB)
  LLdir <- function(p, n){
    log(ddirichlet(p, n))
  }
  (LLdir(pA, pA*nA) + LLdir(pB, pB*nB)) - (LLdir(pAB, pA*nA) + LLdir(pAB, pB*nB))
}