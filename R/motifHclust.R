#' Hierarchical Clustering motifs
#' @description functions to perfom clustering of output of matalign
#' @param align output of matalign, used to generate distance matrix.
#' @param ... parameter to pass to the \link[stats]{hclust}.
#' @return An object of hclust.
#' @importFrom stats hclust dist
#' @export
#' @examples 
#'  if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'   fp <- system.file("extdata", package="motifStack")
#'   fs <- dir(fp, "pcm$")
#'   pcms <- importMatrix(file.path(fp, fs), format="pcm")
#'   align <- matalign(pcms)
#'   hc <- motifHclust(align, method="average")
#'  }
motifHclust <- function(align, ...){
  d <- matalignOut2dist(align)
  d <- dist(d)
  hc <- do.call("hclust", list(d=d, ...))
  return(hc)
}

matalignOut2dist <- function(align){
  stopifnot(
    "align must be output of matalign"=
      all(c("motif1", "motif2", "distance") %in% colnames(align)))
  d <- align[, c("motif1", "motif2", "distance")]
  motifs <- unique(c(align$motif1, align$motif2))
  d <- matrix(NA, nrow = length(motifs), ncol = length(motifs))
  rownames(d) <- colnames(d) <- motifs
  rownames(align) <- paste(align$motif1, align$motif2)
  for(i in motifs){
    for(j in motifs){
      if(i!=j){
        if(paste(i, j) %in% rownames(align)){
          d[i, j] <- align[paste(i, j), "distance"]
        }else{
          d[i, j] <- align[paste(j, i), "distance"]
        }
      }
    }
  }
  return(d)
}