#' convert DNA motif into RNA motif
#'
#' convert DNA motif into RNA motif
#' @param pfm An object of "pcm" or "pfm"
#' @return An object of "pcm" or "pfm" of RNA motif
#' @export
#' @examples 
#' motifs<-importMatrix(dir(file.path(find.package("motifStack"),
#'                          "extdata"),"pcm$", full.names = TRUE))
#' rnaMotifs <- DNAmotifToRNAmotif(motifs)
DNAmotifToRNAmotif <- function(pfm){
  if(is.list(pfm)){
    return(lapply(pfm, DNAmotifToRNAmotif))
  }
  stopifnot(inherits(pfm, c("pcm", "pfm")))
  pfm@alphabet <- "RNA"
  rn <- rownames(pfm@mat)
  rn[rn=="T"] <- "U"
  rownames(pfm@mat) <- rn
  pfm@background["U"] <- pfm@background["T"]
  pfm@background <- pfm@background[rownames(pfm@mat)]
  pfm@color["U"] <- pfm@color["T"]
  pfm@color <- pfm@color[rownames(pfm@mat)]
  pfm
}