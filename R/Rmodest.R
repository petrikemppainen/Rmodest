#' Rmodest
#'
#' Rmodest integrates results across replicated STRUCTURE runs spanning a large range of K values (assumed number of populations) and takes into account multimodality of the data.
#' 
#' \code{Rmodest} impliements modest analyses (MODe ESTimation). The input for modest is a collection of STRUCTURE output files (f-files) deriving from replicate runs over a large range of Ks. Given a path to a folder containing these files, the function \code{\link{importData}} extracts relevant information for subsequent modest analyses. In short, all Q-columns are appended into a single matrix and any other relevant information for modest analyses are stored. Next, a dissimilarity measure (\code{Qdist}) is estimated between all the Q-colums (by function \code{\link{getQdist}}) such that results from cluster runs can be compared across replicates as well as between runs at different Ks. The function \code{\link{Qdist2KDAG}} uses the \code{Qdist} to first construct a single linkage clustering tree and transform this tree into a K-directed acyclic graph. Any Q-columns within each K that are more similar than a user defined threshold (expressed as a proportion of the maximum distance between any two Q-clumns in the data) will be merged and considered to represent the same cluster. The results can be visualised by function  \code{\link{plotKDAG}}. This function also assigns a color value (based the \code{Qdist}, inferred by the two first axes of a PCA) for each vertex (here representing a pool of different Q-clumns) and prepares data that can be used by \code{\link{plotHist}}. \code{\link{plotHist}} plots stacked histogram of probabilities that a given individual belongs to a specific cluster. Because of the potential for multimodality of the data an individual can be assingned to more clusters than is defined by K.
#' 
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Stuart Baird \email{stuartj.e.baird@@gmail.com}, Sophie Gryseels \email{sophie.gryseels@@uantwerpen.be}, Benny Borremans \email{benny.borremans@@uantwerpen.be}
#' 
#' @docType package
#' @name Rmodest
NULL
