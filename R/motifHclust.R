#' Hierarchical Clustering motifs
#' @description functions to perfom clustering of output of matalign
#' @param align output of matalign, used to generate distance matrix.
#' @param ... parameter to pass to the \link[stats]{hclust}.
#' @return An object of hclust.
#' @importFrom stats hclust dist setNames
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
  # Extract required columns
  d_df <- align[, c("motif1", "motif2", "distance")]
  motifs <- unique(c(align$motif1, align$motif2))
  
  # Initialize matrix
  d <- matrix(NA, nrow = length(motifs), ncol = length(motifs))
  rownames(d) <- colnames(d) <- motifs
  
  # Create lookup for both directions at once
  lookup <- setNames(d_df$distance, paste(d_df$motif1, d_df$motif2))
  lookup_rev <- setNames(d_df$distance, paste(d_df$motif2, d_df$motif1))
  
  # Vectorized assignment
  pairs <- expand.grid(i = motifs, j = motifs, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$i != pairs$j, ]
  keys <- paste(pairs$i, pairs$j)
  
  # Try forward lookup, then reverse
  distances <- lookup[keys]
  missing <- is.na(distances)
  distances[missing] <- lookup_rev[keys[missing]]
  
  # Assign to matrix
  d[cbind(pairs$i, pairs$j)] <- distances
  
  return(d)
}
