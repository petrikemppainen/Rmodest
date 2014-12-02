#' Parses STRUCTURE output to a K-directed acyclig graph (KDAG)
#' 
#' Producdes a KDAG from  on single linkage clustering tree based on a pairwise matrix of distances between Q-columns (produced by \code{\link{getQdist}}).
#' 
#' More info to come
#' 
#' @param Qdist A matrix of pairwise distanances between Q-columns produced by \code{\link{getQdist}}
#' @param K A vector of K's for each Qcolumn in \code{Qcols}
#' @param threshold Gives the minimum value for the distance below which Q-columns are merged, expressed as fraction of the maximum distance between any two Q-columns (clusters) in the data set (default=0.05).
#' @return returns an edge list representing the KDAG
#' @export
#' @import ape
#' @keywords Qdist2KDAG, getQdist
#' @seealso \code{\link{getQdist}} and \code{\link{plotKDAG}}
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K)


Qdist2KDAG <- function(Qdist, K, threshold=0.05){
  out <- mergeDAG_Qdist(Qdist2KDAG_inner(Qdist, K), Qdist, K, threshold=threshold)
}  


##############################
### Not exported functions ###
##############################

Qdist2KDAG_inner <- function(Qdist, K){
  ### prepare slink (single linkage clulstering tree) ###
  if(is.na(Qdist[2,1])){
    Qdist <- t(Qdist)
  }
  
  slink <- as.phylo(hclust(as.dist(Qdist), method="single")) # get slink from distance matrix and parse it to a 'phylo' object
  
  ### get necessary data for parsing ###
  Ntips <- length(slink$tip.label) 
  Nclust <- slink$Nnode
  el <- slink$edge
  
  # a list of all possible sub-cluster from slink using function 'extract.clade' from 'ape' package
  Cl <- lapply((Ntips+1):max(el[,1]), function(x) extract.clade(slink, x)) 
  
  Cl <- lapply(Cl, function(x) as.numeric(x$tip.label)) 
  Cl <- c(1:Ntips, Cl) 
  Cl_K <- lapply(Cl, function(x) K[x]) 
  Cl_minK <- sapply(Cl_K, min) 
  Cl_minK_Qcol <- lapply(1:length(Cl), function(x) sort(Cl[[x]][Cl_K[[x]] == Cl_minK[[x]]])) 
  
  ### The parser ###
  DAG <- matrix(NA, 0, 2) 
  Eindex <- 0
  index <- NULL 
  parser <- function(x){ 
    ## prepare info
    children <- el[,2][el[,1]==x] # get children
    children <- children[order(Cl_minK[children])] # order children such that delta K is never negative
    # get Ks for children
    child1_K <- unlist(Cl_minK[children][1])
    child2_K <- unlist(Cl_minK[children][2])
    # Do not create edges when K for children is the same, this is the only exception (!)
    if(child1_K!=child2_K){
      index <<- c(index, x)
      # get grandchildren and vertex names
      grandChildren1 <- unlist(Cl_minK_Qcol[children[1]])
      grandChildren2 <- unlist(Cl_minK_Qcol[children[2]])
      child1_names <- grandChildren1
      child2_names <- grandChildren2
      #two nested functions that create all edges
      parseInner <- function(y){
        ## prepare some info
        # get Ks for children
        createEdges <- function(x) {
          createEdgesInner <- function(y){
            from <- child1_names[y]
            to <- child2_names[x]      
            # if delta K ==1 create edge
            # else create necessary empty vertices
            if(diff(c(child1_K, child2_K))==1){ 
              DAG <- rbind(DAG, c(from, to))
            }else{ 
              Eindex <<- Eindex-1
              DAG <- rbind(DAG, c(Eindex, to))
              if(diff(c(child1_K, child2_K))>2){
                for(i in 1:(diff(c(child1_K, child2_K))-2)){
                  DAG <- rbind(DAG, c(Eindex-1, Eindex))
                  Eindex <<- Eindex-1                
                }
              }
              DAG <- rbind(DAG, c(from, Eindex))
            }
          }
          DAG <- rbind(DAG, do.call('rbind', lapply(1:length(child1_names), createEdgesInner)))  
          return(DAG)
        }
        DAG <- do.call('rbind', lapply(1:length(child2_names), createEdges))
        return(DAG)
      }
      DAG <- parseInner()  
    }
  }
  
  ### now parse and prepare output ###
  DAGlist <- lapply((Ntips+1):(Ntips+Nclust), parser) 
  
  names(DAGlist) <- index 
  DAG <- do.call('rbind', DAGlist) 
  DAG <- DAG[!duplicated(DAG),]
  names <- unlist(Cl[which(K==2)])
  DAG <- rbind(DAG, cbind(0, names))
  return(DAG)
}


mergeDAG_Qdist <- function(DAG, Qdist, K, threshold=0.05){
  ### prepare data ###
  if(is.na(Qdist[1,2])){
    Qdist <- t(Qdist)
  }
  
  ### For each K find set of Q-columns (clusters) that are not connected by any value of Q dist below threshold ###
  G_clusters <- list()
  for(i in unique(K)){
    G_K <- Qdist[K==i, K==i]
    g <- graph.adjacency(G_K, mode="upper", diag=FALSE, weighted=TRUE)          
    g <- delete.edges(g, which(E(g)$weight>threshold))
    G_clusters[[i]] <- lapply(decompose.graph(g), function(x) V(x)$name)  
  }
  
  ### now parse the orignal full KDAG into a merged KDAG based on file G_clusters ###
  new.parents <- list(0)
  names(new.parents[[1]]) <- "root"
  currentK <- 1
  MDAG <- matrix(NA, 0, 2)
  Eindex <- 0
  while((currentK)<max(K)){
    parents <- unique(new.parents)
    new.parents <- list()
    currentK <- currentK+1
    for(i in 1:length(parents)){
      if(any(DAG[,1] %in% unlist(parents[[i]])) && any(unique(DAG[,2][DAG[,1] %in% unlist(parents[[i]])])>0)){
        from <- paste(sort(as.numeric(unlist(parents[[i]]))), collapse=",")
        children <- unique(DAG[,2][DAG[,1] %in% unlist(parents[[i]])])
        vertices <- G_clusters[[currentK]][sapply(G_clusters[[currentK]], function(x) any(x %in% children))]
        for(j in 1:length(vertices)){
          to <- paste(sort(as.numeric(unlist(vertices[j]))), collapse=",")
          MDAG <- rbind(MDAG, c(from, to))
          new.parents <- c(new.parents, vertices[j])
        }
      }
    }
    hanging <- G_clusters[[currentK]][!sapply(G_clusters[[currentK]], function(x) any(x %in% unlist(new.parents)))]
    if(length(hanging)>0){
      ## build in hanging vertices
      for(k in 1:length(hanging)){
        vertices <- hanging[[k]]
        vertices <- unique(DAG[,1][DAG[,2] %in% vertices])
        while(all(vertices<0)){
          vertices <- unique(DAG[,1][DAG[,2] %in% vertices])
        }
        vertices <- as.character(vertices)
        from <- strsplit(unique(as.vector(MDAG)), ",")
        from <- from[sapply(from, function(x) any(x %in% vertices))]
        Kfrom <- unique(K[as.numeric(unlist(from))])
        from <- lapply(from, function(x) paste(sort(as.numeric(unlist(x))), collapse=","))
        Kto <- unique(K[as.numeric(hanging[[k]])])
        to <- paste(hanging[[k]], collapse=",")
        
        DAG.temp <- matrix(NA, 0, 2) 
        for(j in 1:length(from)){
          if(diff(c(Kfrom, Kto))==1){ 
            DAG.temp <- rbind(DAG.temp, c(from[[j]], to))
          }else{ 
            Eindex <- Eindex-1
            DAG.temp <- rbind(DAG.temp, c(Eindex, to))
            if(diff(c(Kfrom, Kto))>2){
              for(i in 1:(diff(c(Kfrom, Kto))-2)){
                DAG.temp <- rbind(c(Eindex-1, Eindex), DAG.temp)
                Eindex <- Eindex-1                
              }
            }
            DAG.temp <- rbind(c(from[[j]], Eindex), DAG.temp)
          }
        }  
        MDAG <- rbind(MDAG, DAG.temp)
        new.parents <- c(new.parents, hanging[k])
      }
    }    
  }
  MDAG
}