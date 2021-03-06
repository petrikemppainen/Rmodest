#' Imports STRUCTURE data
#'
#' Imports STRUCTURE data from files given a path to folder
#' 
#' Given a path to a folder that contains f-files from replicate STRUCTURE runs spanning different Ks (number of assumed populations), \code{importData} extracts all necessary information for \code{\link{Rmodest}}.
#'
#' @param input.data.folder Path to folder that contains f-files generated by STRUCTURE. Each f-file must contain a single unique numeric indentifier but may also contain any other characters.
#' @keywords importData, STRUCTURE, Rmodest
#' @seealso \code{\link{getQdist}} and \code{\link{Qdist2KDAG}}
#' @return 
#' Returns a list with the following objects:
#' \itemize{
#'  \item \code{Qcols}: A matrix that contains all Q-columns ordered by increasing K (individuals as rows in the same order as in the original data)
#'  \item \code{allele_counts}: A matrix that contains allelel counts for each locus in the data set (weighted according to missing data) for each Q-column in \code{Qcols}
#'  \item \code{K}:  K (asssumed number of populations) for each Q-column in \code{Qcols} and \code{allele_counts}
#'  \item \code{NinK}:  Index of each Q-column within K
#'  \item \code{NinR}:  Index of each Q-column within a STRUCTURE run/replicate
#'  \item \code{f_file}:  The name of the file from which the Q-column was extracted from
#'  \item \code{file_nr}:  The unique number identifier for each Qcolumn extracted from \code{f_file}
#' }
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @export
#' @examples
#' \dontrun{
#' modest <- importData("path/to/folder/with/f_files")}


importData <- function(input.data.folder){
  ### check and prepare format ###
  
  if(read.fwf(textConnection(input.data.folder), widths=c(nchar(input.data.folder)-1, nchar(input.data.folder)))[2] != "/"){
    input.data.folder <- paste(input.data.folder, "/", sep="")
  }
  # end check format
  input.files <- paste(input.data.folder, list.files(input.data.folder),sep="")
  
  
  ### get allele counts which are used for calculating Qdist ###
  out <- list()
  inp <- suppressWarnings(readLines(input.files[1]))
  Nind <- as.numeric(gsub(" individuals", "", inp[grep(" individuals",inp)][1]))
  rstarts <- grep("Locus",inp)+1
  Nalleles <- as.numeric(gsub(" alleles", "",inp[rstarts]))
  one_pop <- NULL
  
  ### get allele frequencies
  for(i in 1:length(input.files)){
    inp <- suppressWarnings(readLines(input.files[i]))
    rstarts <- grep("Locus",inp)+1
    Npop <- as.numeric(gsub(" populations assumed", "", inp[grep(" populations assumed",inp)]))
    if(Npop>1){
      out.temp <- as.list(rep(NA, length(rstarts)))    
      for(j in 1:length(rstarts)){
        temp <- do.call('rbind', strsplit(inp[(rstarts[j]+2):(rstarts[j]+1+Nalleles[j])], " ", fixed=TRUE))
        temp <- temp[,(ncol(temp)-Npop+1):ncol(temp)]
        if(is.vector(temp)){
          temp <- as.numeric(temp)
        }else{temp <- apply(temp, 2, as.numeric)}
        out.temp[[j]] <- list(as.numeric(gsub("% missing data", "", inp[rstarts[j]+1])), temp)
      }
      out[[i]] <- out.temp    
    }else{
      out[[i]] <- "one_pop"
      one_pop <- c(one_pop, i)
    }
  }
  
  ###  Get counts
  if(is.null(one_pop)){
    allele_counts <- do.call('cbind', lapply(out, function(y) do.call('rbind', lapply(y, function(x) x[[2]]*Nind*(1-x[[1]])*2))))
  }else{
    allele_counts <- do.call('cbind', lapply(out[-one_pop], function(y) do.call('rbind', lapply(y, function(x) x[[2]]*Nind*(1-x[[1]])*2))))
  }
  colnames(allele_counts) <- c(1:ncol(allele_counts))
  rownames(allele_counts) <- rep(1:length(Nalleles), Nalleles)
  
  
  ### get Qcols.raw, the core of this code is written by Benny Borremans ###
  Qcols.raw <- list()
  for(i in 1:length(input.files)){
    inp <- suppressWarnings(readLines(input.files[i]))
    N <- as.numeric(gsub(" |individuals","",inp[grep("individuals",inp)[which(grep("individuals",inp)>grep("Run parameters",inp))][1]]))
    rstart <- grep("(%Miss)",inp)+1
    rstop <- grep("(%Miss)",inp)+N 
    inp2 <- gsub(",|\\(|\\)"," ",inp[rstart:rstop])
    dat <- read.table(text=inp2)
    if(grepl("Label",inp[rstart-1])){idx.id <- 2}else{idx.id <- 1}    # check if there are individual labels, which changes the columns that need to be extracted
    idx.first.cluster.column <- which(dat[1,]==":")+1
    dat.export <- dat[,c(idx.id,idx.first.cluster.column:ncol(dat))]
    colnames(dat.export)<- c("Label",paste(c(rep("C",(ncol(dat.export)-1)/3),rep(c("LCI","UCI"),(ncol(dat.export)-1)/3)),c(1:((ncol(dat.export)-1)/3),rep(1:((ncol(dat.export)-1)/3),each=2)),sep=""));colnames(dat.export)
    dat.export$Label <- as.character(dat.export$Label)
    Qcols.raw[[i]]<- dat.export
  }
  
  
  
  ### get Qcols and additional information necessary for modoest analyses ###
  fColnames <- lapply(Qcols.raw, function(x) colnames(x))
  extractIndex <- lapply(fColnames, function(y) which(sapply(strsplit(y, "", fixed=T), function(x) x[1])=="C"))
  temp <- lapply(1:length(extractIndex), function(x) Qcols.raw[[x]][,extractIndex[[x]]])
  Qcols <- do.call('cbind', temp)
  temp2 <- sapply(temp, ncol)
  for(i in 1:length(temp2)){
    if(is.null(temp2[[i]])) temp2[[i]] <- 1
  }
  K <- rep(sapply(extractIndex, length), sapply(extractIndex, length))
  f_file <- do.call('rbind', strsplit(input.files, "/"))
  f_file <- f_file[,ncol(f_file)]
  f_file <- rep(f_file, sapply(extractIndex, length))
  f_file <- f_file[K!=1]
  Qcols <- Qcols[,K!=1]
  Qcols <- as.matrix(Qcols)
  K <- K[K!=1]
  
  ### order, filenames must contain unique numbers... ###
  temp <- gsub("[[:alpha:]]", "", f_file)
  file_nr <- suppressWarnings(as.numeric(gsub("[[:punct:]]", "", temp)))
  new.order <- order(K, file_nr)
  Qcols <- Qcols[,new.order]
  colnames(Qcols) <- 1:ncol(Qcols)
  allele_counts <- allele_counts[,new.order]
  f_file <- f_file[new.order]
  file_nr <- file_nr[new.order]
  K <- K[new.order]
  
  #### prepare output ###
  NinR <- unlist(sapply(unique(file_nr), function(x) c(1:length(file_nr[file_nr==x]))))
  NinK <- unlist(sapply(unique(K), function(x) c(1:length(K[K==x]))))
  out <- list(Qcols, allele_counts, K, NinR, NinK, f_file, file_nr)
  names(out) <-  c("Qcols", "allele_counts", "K", "NinR", "NinK", "f_file", "file_nr")
  out
}


