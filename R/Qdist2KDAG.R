#' Parses STRUCTURE output to a K-directed acyclig graph (KDAG)
#' 
#' Producdes a KDAG from  on single linkage clustering tree based on a pairwise matrix of distances between Q-columns (produced by \code{\link{getQdist}}).
#' 
#' More info to come
#' 
#' @param Qdist A matrix of pairwise distanances between Q-columns produced by \code{\link{getQdist}}
#' @param K A vector of K's for each Qcolumn in \code{Qcols}
#' @param merge.by Weather to merge by topology or by Qdist
#' @param threshold If `merge.by="Qdist"`, threshold gives the minimum value for the distance below which Q-columns are merged expressed as fraction of the maximum distance in the data set (default=0.05).
#' @return returns an edge list representing the KDAG
#' @export
#' @import igraph
#' @import ape
#' @keywords Qdist2KDAG, getQdist
#' @author Petri Kemppainen, \email{petrikemppainen2@@gmail.com} and Stuart Baird, \email{stuartj.e.baird@@gmail.com}
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K, merge.by="Qdist")


Qdist2KDAG <- function(Qdist, K, merge.by="topology", threshold=0.05){
  if(merge.by=="topology"){
    out <- mergeDAG_top(Qdist2KDAG_inner(Qdist, K), K)
  }
  if(merge.by=="Qdist"){
    out <- mergeDAG_Qdist(Qdist2KDAG_inner(Qdist, K), Qdist, K, threshold=threshold)
  }
  out
} 

Qdist2KDAG_inner <- function(Qdist, K){ # the imputs for this function are a matrix of Qcols and a vector of K for each Qcolumn
  ### prepare slink
  if(is.na(Qdist[2,1])){
    Qdist <- t(Qdist)
  }
  
  slink <- as.phylo(hclust(as.dist(Qdist), method="single")) # get slink from distance matrix and parse it to a 'phylo' object
  
  ### prepare some useful info 
  Ntips <- length(slink$tip.label) # number of tips
  Nclust <- slink$Nnode # number of clusters
  el <- slink$edge # shorter name for the edge list which contains the slink tree structure
  
  ### prepare info to be used by the parser
  Cl <- lapply((Ntips+1):max(el[,1]), function(x) extract.clade(slink, x)) # a list of all possible sub-cluster from slink using function 'extract.clade' from 'ape' package
  Cl <- lapply(Cl, function(x) as.numeric(x$tip.label)) # list with tip labels (Qindex) for each cluster 
  Cl <- c(1:Ntips, Cl) # now add all tips as well 
  Cl_K <- lapply(Cl, function(x) K[x]) # List with K:s for each cluster/tip
  Cl_minK <- sapply(Cl_K, min) # minimum K for each cluster/tip
  Cl_minK_Qcol <- lapply(1:length(Cl), function(x) sort(Cl[[x]][Cl_K[[x]] == Cl_minK[[x]]])) # Q cols that have minK for each cluster; these will be vertices in the DAG
  DAG <- matrix(NA, 0, 2) # file to append to
  Eindex <- 0
  index <- NULL # gives the names for each cluster for which the parser was used
  #x <- 360
  ### The parser
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
            if(diff(c(child1_K, child2_K))==1){ # if delta K ==1
              DAG <- rbind(DAG, c(from, to))
            }else{ # if delta K >1 create necessary empty vertices
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
  ### now parse
  DAGlist <- lapply((Ntips+1):(Ntips+Nclust), parser) # use parser for all clusters
  names(DAGlist) <- index # give each object relevant names (equals to node numbers for which the parser was used), used for debugging
  DAG <- do.call('rbind', DAGlist) # append 
  DAG <- DAG[!duplicated(DAG),]
  ##create root
  names <- unlist(Cl[which(K==2)])
  DAG <- rbind(DAG, cbind(0, names)) # add root
  return(DAG)
}
##################
mergeDAG_top <- function(DAG, K){
  #start with root
  temp <- list()
  temp[[1]] <- list(0)
  names(temp[[1]][[1]]) <- "root"
  currentK <- 1
  MDAG <- matrix(NA, 0, 2)
  while((currentK)<max(K)){
    parents <- temp
    currentK <- currentK+1
    temp <- list()
    x <- 0
    #i <- 11
    #collapse nodes with comparable topologies and build merged DAG
    for(i in 1:length(parents)){
      if(any(DAG[,1] %in% unlist(parents[[i]]))){
        x <- x+1
        from <- paste(sort(unlist(parents[[i]])), collapse=",")
        children <- unique(DAG[,2][DAG[,1] %in% unlist(parents[[i]])])
        subDAGs <- lapply(children, function(x) extractDAG(DAG, x))
        subDAGsTop <- lapply(subDAGs, function(x) collapseSinglesEl(x))
        trees <- subDAGsTop
        which.not.tip <- which(subDAGsTop != "tip")
        trees[which.not.tip] <- lapply(subDAGsTop[which.not.tip], DAG2Phylo)
        d.tree = matrix( nrow=length(trees), ncol=length(trees), dimnames=list(children, children) )
        #decision tree for what two trees have comparable topologies
        for(h in 1:length(trees)){
          for(k in h:length(trees)){
            if(any(c(any(trees[[k]]=="tip"), any(trees[[h]]=="tip")))){
              if(all(trees[[k]]=="tip") && any(trees[[h]]=="tip")){
                d.tree[k,h] <- 0 
              }else{
                d.tree[k,h] <- 1 
              }
            }else{
              if(trees[[k]]$Nnode+trees[[h]]$Nnode==2){
                if(length(trees[[k]]$tip.label)==length(trees[[h]]$tip.label)){
                  d.tree[k,h] <- 0
                }else{
                  #  if(singleNodeComparable==TRUE){
                  d.tree[k,h] <- 1
                  #  }else{d.tree[k,h] <- 1}
                }
              }
              if(trees[[k]]$Nnode+trees[[h]]$Nnode==3) d.tree[k,h] <- 1  
              if(trees[[k]]$Nnode+trees[[h]]$Nnode>3){
                if(any(c(trees[[k]]$Nnode,trees[[h]]$Nnode)==1)){
                  d.tree[k,h] <- 1
                }else{
                  d.tree[k,h] <- dist.topo(trees[[k]], trees[[h]])  
                }
              } 
            }
          }
        }
        #find clusters of the same topology, here I use networks
        d.tree[d.tree!=0] <- -1
        d.tree[d.tree==0] <- 1
        g <- graph.adjacency(d.tree, mode="lower", diag=FALSE, weighted=TRUE)
        g <- delete.edges(g, which(E(g)$weight==-1))
        temp2 <- lapply(decompose.graph(g), function(x) V(x)$name)
        #
        for(j in 1:length(temp2)){
          if(any(unlist(temp2[j]) %in% unlist(temp))){
            shared <- unlist(temp2[j])[unlist(temp2[j]) %in% unlist(temp)]
            which.shared <- which(unlist(lapply(temp, function(x) any(x %in% shared))))
            new.name <- paste(sort(unique(c(unlist(temp[which.shared]), unlist(temp2[j])))), collapse=",")
            for(y in which.shared){
              MDAG[,2][MDAG[,2] == paste(unlist(temp[y]), collapse=",")] <- new.name
              temp[[y]] <- sort(unique(c(unlist(temp[which.shared]), unlist(temp2[j]))))
            }
            MDAG <- rbind(MDAG, c(from, new.name))  
            if(j==length(temp2)) x <- x-1
          }else{
            columns <- sort(unlist(temp2[j]))
            names(columns) <- NULL
            temp[[x]] <- columns
            to <- paste(columns, collapse=",")
            MDAG <- rbind(MDAG, c(from, to))
            if(j<length(temp2)) x <- x+1
          }
        }
        #find which columns should be merged and build merged DAG
      }
    }
    temp
  }
  out <- apply(MDAG, 2, function(y){
    temp <- lapply(strsplit(y, ","), as.numeric)
    temp <- lapply(temp, function(x) if(!all(x<0)) {x[x>=0]}else{-1})
    temp[sapply(temp, function(x) any(x<0))] <- 1:length(which(sapply(temp, function(x) any(x<0))))*-1
    temp <- sapply(temp, function(x) paste(sort(x), collapse=","))  
  })  
}
##################
mergeDAG_Qdist <- function(DAG, Qdist, K, threshold=0.05){
  if(is.na(Qdist[1,2])){
    Qdist <- t(Qdist)
  }
  Qdist[1:5, 1:5]
  G_clusters <- list()
  for(i in unique(K)){
    G_K <- Qdist[K==i, K==i]
    g <- graph.adjacency(G_K, mode="upper", diag=FALSE, weighted=TRUE)          
    g <- delete.edges(g, which(E(g)$weight>threshold))
    G_clusters[[i]] <- lapply(decompose.graph(g), function(x) V(x)$name)  
  }
  
  #start with root
  new.parents <- list(0)
  names(new.parents[[1]]) <- "root"
  currentK <- 1
  MDAG <- matrix(NA, 0, 2)
  Eindex <- 0
  while((currentK)<max(K)){
    parents <- unique(new.parents)
    new.parents <- list()
    currentK <- currentK+1
    #i <- 1
    #collapse nodes with comparable topologies and build merged DAG
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
        
        DAG.temp <- matrix(NA, 0, 2) # file to append to
        for(j in 1:length(from)){
          if(diff(c(Kfrom, Kto))==1){ # if delta K ==1
            DAG.temp <- rbind(DAG.temp, c(from[[j]], to))
          }else{ # if delta K >1 create necessary empty vertices
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
############## Functions
extractDAG <- function(DAG, node){
  el <- DAG
  node.first <- node
  vertices <- NULL
  while(!identical(el[,2][el[,1] %in% node], numeric(0))){
    vertices <- c(vertices, el[,2][el[,1] %in% node])
    node <- el[,2][el[,1] %in% node]
  } 
  vertices <- unique(c(node.first, vertices))
  subDAG <- el[el[,1] %in% vertices & el[,2] %in% vertices,]
  subDAG  
}
###
DAG2Phylo <- function(DAG){
  el <- DAG
  el.new <- el
  tips <- unique(el[,2][!el[,2] %in% el[,1]])
  nodes <- unique(as.vector(el)[!as.vector(el) %in% tips])
  nTips <- length(tips)
  nNodes <- length(nodes)
  el.new[,2] <- c((nTips+1):(nNodes+nTips))[match(el[,2], nodes)]
  el.new[,1] <- c((nTips+1):(nNodes+nTips))[match(el[,1], nodes)]
  el.new[,2][is.na(el.new[,2])] <- na.omit(c(1:nTips)[match(el[,2], tips)])
  tree <- list()
  tree$edge <- apply(el.new, 2, as.numeric)
  tree$tip.label <- tips
  tree$node.label <- nodes
  tree$Nnode <- nNodes
  class(tree) <- "phylo"
  tree
}
#####
collapseSinglesEl <- function(el){
  
  if(is.vector(el)){return("tip")
  }else{
    temp <- table(el[,1])
    singles <- as.numeric(names(temp)[which(temp==1)])
    
    if(length(temp)==length(singles)){return("tip")
    }else{
      for(i in singles){
        el[,2][el[,2] == i] <- el[,2][el[,1] == i]
        el <- el[el[,1] != i,]
      }
      el <- el[!duplicated(el),]
      if(is.vector(el)){
        return("tip")
      }else{
        return(el)  
      } 
    }
  }
}
############

