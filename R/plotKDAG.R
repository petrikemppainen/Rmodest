#' Plots K-directed acyclic graph (KDAG)
#'
#' Plots KDAG based on information produced by \code{\link{importData}}, \code{\link{getQdist}} and \code{\link{Qdist2KDAG}} and exports vertex information for function \code{\link{plotHist}}.
#' 
#' To come
#'
#' @param modest A file produced by \code{importData} that contains parsed data from replicate STRUCTURE runs
#' @param Qdist A file produced by \code{getQdist} that contains all pairwise distances between the Q-columns in file \code{modest}
#' @param KDAG A file produced by \code{Qdist2KDAG} that contains an edge list describing a KDAG.
#' @param min.vertex.size To come
#' @param max.vertex.size To come
#' @param intensities To come
#' @param plot.PCA To come
#' @param hue.saturation To come
#' @param use.rgb Color vertices by the three first axes of RGB space
#' @keywords plotKDAG, plotDAG, STRUCTURE, importData, getQdist
#' @seealso \code{\link{importData}}, \code{\link{getQdist}} and \code{\link{Qdist2KDAG}}
#' @return a list with the following objects. \code{vertices}: a list with each vertex in the DAG as an object containing index of each Q-column that have been merged into that vertex. \code{K_vertices}: the K for each vertex in file \code{vertices}. \code{color}: a vector of color codes that were used for each vertex in the KDAG plot.
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @export
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K)
#' plotKDAG(modest_MCMV, Qdist_MCMV, KDAG) # default, more options available


plotKDAG <- function(modest, Qdist, KDAG, min.vertex.size=4, max.vertex.size=20, intensities=c(1, 0.6, 0.3), plot.PCA=NULL, hue.saturation=c(1,2), use.rgb=FALSE){
  NinR <- modest$NinR
  NinK <- modest$NinK
  K <- modest$K
  rep <- modest$file_nr
  
  el <- modifyKDAG_corr(KDAG, K)
  g <- graph.edgelist(el)
  new_vertices <- unique(V(g)$name)
  new_vertices <- strsplit(new_vertices, ",")
  names(new_vertices) <- c(1:length(new_vertices))
  new_vertices <- suppressWarnings(lapply(new_vertices, as.numeric))
  
  which.to.use <- sapply(new_vertices, function(x) !any(x<0))
  which.to.use[is.na(which.to.use)] <- FALSE
  which.to.use[1] <- FALSE
  
  G_temp <- Qdist
  G_temp[lower.tri(G_temp)] <- t(G_temp)[lower.tri(G_temp)]
  
  pair <- combn(1:length(new_vertices[which.to.use]), 2, simplify = F)
  
  means <- lapply(pair, function(x){
    temp <- new_vertices[which.to.use][c(x[1],x[2])]
    temp2 <- G_temp[temp[[1]], temp[[2]]]
    mean(na.omit(temp2))
  })
  
  g.temp=graph.data.frame(as.data.frame(cbind(do.call(rbind, pair))))
  E(g.temp)$weight <- unlist(means)
  means <- get.adjacency(g.temp,sparse=FALSE, attr="weight")
  dimnames(means) <- list(names(new_vertices[which.to.use]), names(new_vertices[which.to.use]))
  
  PCA <- pcoa(as.dist(t(means)), correction="none", rn=NULL)
  PCA$vectors[,1:3] <- apply(PCA$vectors[,1:3], 2, function(x) (x - min(x))/diff(range(x)))
  vertices <- new_vertices[which.to.use]
  K_vertices <- sapply(vertices, function(x) unique(K[x]))
  
  if(use.rgb==FALSE){
    ########## prep intensity
    ## get modes
    
    modes <- list()
    i <- 6
    for(i in unique(K_vertices)-1){
      vertices.temp <- vertices[K_vertices==i+1]
      vertices.rep <- lapply(vertices.temp, function(x) rep[x])
      vertices.rep.unique <- unique(unlist(vertices.rep))
      modes.temp <- sapply(vertices.rep.unique, function(x) paste(sort(names(vertices.rep)[sapply(vertices.rep, function(y) any(y==x))]), collapse=","))
      mode.table <- table(modes.temp)
      mode.table <- mode.table[rev(order(mode.table))]
      temp <-  strsplit(names(mode.table), ",")
      temp <- lapply(temp, as.numeric)
      names(temp) <- as.vector(mode.table)
      modes[[i]] <- temp
    }
    names(modes) <- unique(K_vertices)
    
    #get intensity 
    intensity <-  rep(NA, length(V(g)))
    pre.int <- intensities
    
    for(i in 1:length(modes)){
      to.color <- as.vector(unlist(modes[[i]]))
      names.temp <- names(modes[[i]])
      counts <- unique(names.temp)  
      for(j in 1:length(names.temp)){
        temp <- intensity[unlist(modes[[i]][names.temp==names.temp[j]])]
        new.vert <- unlist(modes[[i]][names.temp==names.temp[j]])[is.na(temp)]
        if(length(counts)==2){
          intensity[new.vert] <- pre.int[if(j<2) j else 3]
        }else{
          intensity[new.vert] <- pre.int[if(j<3) j else 3]
        }
      }
    }
    intensity <- intensity[which.to.use]
    
    ##### get colors
    #hue.saturation <- c(2,1)
    Colors <- sapply(1:nrow(PCA$vectors), function(x) hsv(0.05+PCA$vectors[x,hue.saturation[1]]*0.95, 0.05+PCA$vectors[x,hue.saturation[2]]*0.95, intensity[x], 1))
    names(Colors) <- V(g)$name[which.to.use]
    
  }else{
    Colors <- rgb(PCA$vectors[,1:3])
    names(Colors) <- V(g)$name[which.to.use]
  }
  
  ############### plot PCA
  if(!is.null(plot.PCA)){
    plot(PCA$vectors[,plot.PCA], bg=Colors, pch=21, cex=2, col="black")
  }
  
  ##########
  type <- V(g)$name
  type[sapply(strsplit(V(g)$name, " "), length)==2] <- "Extra"
  type[which(suppressWarnings(as.numeric(V(g)$name)<0))] <- "Empty"
  type[!type %in% c("Empty", "Extra")] <- "Normal"
  type[V(g)$name=="0"] <- "root"
  ####### prep sizse
  temp <- sapply(new_vertices[type %in% "Normal"], length)
  
  scaleValues <- function(x, a, b, max, min){
    ((b-a)*(x-min))/(max-min)+a
  }
  #max.vertex.size <- 20
  #min.vertex.size <- 4
  V(g)$size[type %in% "Normal"] <- scaleValues(temp, min.vertex.size, max.vertex.size, max(temp), min(temp))
  #V(g)$size[type %in% "Normal"] <- mode.length[which.to.use]
  #########
  V(g)$size[type=="Extra"] <- 0
  V(g)$size[type=="Empty"] <- 3
  V(g)$size[type=="root"] <- 5
  V(g)$shape <- ifelse(type=="Extra", "none", "circle")
  V(g)$shape[type=="Empty"] <- "square"
  V(g)$color[type=="Empty"] <- "white"
  V(g)$shape[type=="root"] <- "square"
  V(g)$color[type=="Extra"] <- NA
  V(g)$color[type=="root"] <- "black"
  E(g)[to(V(g)$name[type!="Extra"])]$color <- "black"
  E(g)$edge.arrow.size <- 0
  #########
  
  V(g)$color[match(names(Colors), V(g)$name)] <- Colors
  lay2 <- layout.sugiyama(g, attributes="all", maxiter = 100)
  plot(g, layout=lay2$layout, edge.arrow.size=0, edge.label.family="Helvetica", vertex.label=NA)
  out <- list(vertices, K_vertices, Colors)
  names(out) <- c("vertices", "K_vertices", "color")
  out
}

##############
modifyKDAG_corr <- function(MDAG, K){  
  DAG.temp <- MDAG
  tips <- DAG.temp[,2][!DAG.temp[,2] %in% DAG.temp[,1]]
  # which tips have K<maxK
  tips_notMaxK <- tips[sapply(strsplit(tips, ","), function(x) unique(K[as.numeric(x)]))!=max(K)]
  # build new edges to extend these to max K
  
  Eindex <- 0
  j <- 1
  if(length(tips_notMaxK)!=0){
    for(j in 1:length(tips_notMaxK)){
      Eindex <- Eindex+1
      tip <- tips_notMaxK[j]
      new_edges <- matrix(NA, 0, 2)
      new_edges <- rbind(new_edges, c(tip, paste("E", Eindex)))
      if(max(K)-unique(K[as.numeric(unlist(strsplit(tip, ",")))])>1){  
        for(i in 1:(max(K)-unique(K[as.numeric(unlist(strsplit(tip, ",")))])-1)){
          new_edges <- rbind(new_edges, c(paste("E", Eindex), paste("E", Eindex+1)))
          Eindex <- Eindex+1
        }
      }
      DAG.temp <- rbind(DAG.temp, new_edges)
    }  
  }
  DAG.temp
}


