#' Plots K-directed acyclic graph (KDAG)
#'
#' Plots KDAG based on information produced by \code{\link{importData}}, \code{\link{getQdist}} and \code{\link{Qdist2KDAG}} and exports vertex information for function \code{\link{plotHist}}.
#' 
#' To come
#'
#' @param modest a file produced by \code{importData} that contains parsed data from replicate STRUCTURE runs
#' @param Qdist a file produced by \code{getQdist} that contains all pairwise distances between the Q-columns in file \code{modest}
#' @param KDAG a file produced by \code{Qdist2KDAG} that contains an edge list describing a KDAG.
#' @param min.vertex.size minimum size to be used for any vertex (default=4)
#' @param max.vertex.size maximum size to be used for any vertex (default=20)
#' @param intensities can be used to indicate which vertices represent Q-colums that only occur in rare modes. Is a vector with three numbers; the first give the intensity for the most common run mode, the second for the second most common run mode and the third for the third (or lower) most common rund mode. Default is c(1, 0.6, 0.3). The parameter \code{use.rgb} needs to be set to \code{FALSE}.
#' @param plot.PCA a vector defining which two axes to plot of the PCA. If \code{NULL} (default) no plot is printed.
#' @param hue.saturation if \code{use.rgb} is \code{FALSE}, this parameter can be used to decide which of the two first axes of the PCA is used for hue and saturation, respectively. Default is c(1,2) meaning that the first axes is used to decide hue for the a vertex and the second is used for saturaton.
#' @param use.rgb if \code{TRUE} color of vertices is defined by their separation along the three first axes of the PCA in RGB space.
#' @keywords plotKDAG, plotDAG, STRUCTURE, importData, getQdist
#' @seealso \code{\link{importData}}, \code{\link{getQdist}} and \code{\link{Qdist2KDAG}}
#' @return a list with the following objects:
#' \itemize{
#'  \item \code{vertices}: a list with each vertex in the DAG as an object containing index of each Q-column that have been merged into that vertex. 
#'  \item \code{K_vertices}: the K for each vertex in file \code{vertices}. 
#'  \item \code{color}: a vector of color codes for each vertex.
#'  }
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @export
#' @examples
#' data(MODEST)
#' KDAG <- Qdist2KDAG(Qdist_MCMV, modest_MCMV$K)
#' plotKDAG(modest_MCMV, Qdist_MCMV, KDAG) # default, more options available


plotKDAG <- function(modest, Qdist, KDAG, min.vertex.size=4, max.vertex.size=20, intensities=c(1, 0.6, 0.3), plot.PCA=NULL, hue.saturation=c(1,2), use.rgb=FALSE){
  ### prepeare data ###
  NinR <- modest$NinR
  NinK <- modest$NinK
  K <- modest$K
  rep <- modest$file_nr
  
  ### modify dag for plottin and get vertex names ###
  el <- modifyKDAG_corr(KDAG, K)
  g <- graph.edgelist(el)
  new_vertices <- unique(V(g)$name)
  new_vertices <- strsplit(new_vertices, ",")
  names(new_vertices) <- c(1:length(new_vertices))
  new_vertices <- suppressWarnings(lapply(new_vertices, as.numeric))
  which.to.use <- sapply(new_vertices, function(x) !any(x<0))
  which.to.use[is.na(which.to.use)] <- FALSE
  which.to.use[1] <- FALSE
  
  ### get means Qdist for each vertex on which the PCA is based on ###
  G_temp <- Qdist
  G_temp[lower.tri(G_temp)] <- t(G_temp)[lower.tri(G_temp)]
  pair <- combn(1:length(new_vertices[which.to.use]), 2, simplify = F)
  means <- lapply(pair, function(x){
    temp <- new_vertices[which.to.use][c(x[1],x[2])]
    temp2 <- G_temp[temp[[1]], temp[[2]]]
    mean(na.omit(temp2))
  })
  
  g.temp <- graph.data.frame(as.data.frame(cbind(do.call(rbind, pair))))
  E(g.temp)$weight <- unlist(means)
  means <- get.adjacency(g.temp,sparse=FALSE, attr="weight")
  dimnames(means) <- list(names(new_vertices[which.to.use]), names(new_vertices[which.to.use]))
  
  
  ### do PCA ###
  PCA <- pcoa(as.dist(t(means)), correction="none", rn=NULL)
  PCA$vectors[,1:3] <- apply(PCA$vectors[,1:3], 2, function(x) (x - min(x))/diff(range(x)))
  vertices <- new_vertices[which.to.use]
  K_vertices <- sapply(vertices, function(x) unique(K[x]))
  
  ### prepare colors ###
  if(!use.rgb){
    # get modes
    modes <- list()
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
    
    # get colors
    Colors <- sapply(1:nrow(PCA$vectors), function(x) hsv(0.05+PCA$vectors[x,hue.saturation[1]]*0.95, 0.05+PCA$vectors[x,hue.saturation[2]]*0.95, intensity[x], 1))
    names(Colors) <- V(g)$name[which.to.use]
    
  }else{
    # if use.rgb==TRUE get color from the first three axes of the PCA
    Colors <- rgb(PCA$vectors[,1:3])
    names(Colors) <- V(g)$name[which.to.use]
  }
  
  ### plot PCA ###
  if(!is.null(plot.PCA)){
    plot(PCA$vectors[,plot.PCA], bg=Colors, pch=21, cex=2, col="black")
  }
  
  ### prepare data for plotting KDAG ###
  type <- V(g)$name
  type[sapply(strsplit(V(g)$name, " "), length)==2] <- "Extra"
  type[which(suppressWarnings(as.numeric(V(g)$name)<0))] <- "Empty"
  type[!type %in% c("Empty", "Extra")] <- "Normal"
  type[V(g)$name=="0"] <- "root"
  temp <- sapply(new_vertices[type %in% "Normal"], length)
  
  # scale size
  scaleValues <- function(x, a, b, max, min){
    ((b-a)*(x-min))/(max-min)+a
  }
  V(g)$size[type %in% "Normal"] <- scaleValues(temp, min.vertex.size, max.vertex.size, max(temp), min(temp))
  
  # prepare info for KDAG
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
  V(g)$color[match(names(Colors), V(g)$name)] <- Colors

  ### plot KDAG ###
  lay2 <- layout.sugiyama(g, attributes="all", maxiter = 100)
  plot(g, layout=lay2$layout, edge.arrow.size=0, edge.label.family="Helvetica", vertex.label=NA)
  
  ### prepare output ###
  out <- list(vertices, K_vertices, Colors)
  names(out) <- c("vertices", "K_vertices", "color")
  out
}

##############################
### Not exported functions ###
##############################

### this function is necessary for the sugyama layout to behave properly ###
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


