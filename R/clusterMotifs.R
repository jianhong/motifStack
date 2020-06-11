#' cluster motifs
#' @description A help function to do matalign and motifHclust in one function.
#' @param motifs A list of pcms of pfms.
#' @param ... parameter to be passed to \link{matalign} function.
#' @return An object of hclust.
#' @export
#' @examples 
#'  if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'   fp <- system.file("extdata", package="motifStack")
#'   fs <- dir(fp, "pcm$")
#'   pcms <- importMatrix(file.path(fp, fs), format="pcm")
#'   hc <- clusterMotifs(pcms)
#'  }
clusterMotifs <- function(motifs, ...){
  d <- matalign(motifs, ...)
  hc <- motifHclust(d, method="average")
  hc
}

pfmList2matrixList <- function(pfms){
  m <- lapply(pfms, pfm2pwm)
  names(m) <- unlist(lapply(pfms, function(.ele) .ele@name))
  m
}

pfm2pcm <- function(pfm, total=1000){
  stopifnot(is(pfm, "pfm"))
  .mat <- pfm@mat
  .mat <- round(total*.mat+.49)
  new("pcm", name=pfm@name, mat=.mat,
      alphabet=pfm@alphabet, color=pfm@color,
      background=pfm@background, tags=pfm@tags,
      markers=pfm@markers)
}