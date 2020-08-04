#' align AA motifs
#' 
#' align AA motifs for plotting motifs stack
#' 
#' 
#' @param pcms a list of position frequency matrices, pfms must be a list of
#' class pcm
#' @param threshold information content cutoff threshold for useful postions
#' @param minimalConsensus minimal length of consensus for alignment
#' @return a list of aligned motifs
#' @export
#' @examples
#' 
#' pcms<-importMatrix(system.file("extdata", "prot.meme", package="motifStack"),
#'                    format="meme", to="pfm")
#' motifs<-AAmotifAlignment(pcms)
#'
AAmotifAlignment <- function(pcms, threshold=0.4, minimalConsensus = 0){
  if(length(pcms)<2) stop("less than 2 motifs")
  pcms <- lapply(pcms,function(.ele){
    if(!inherits(.ele, c("pfm", "pcm"))) {
      stop("pcms must be a list of pcm objects.")
    }
    if(.ele@alphabet!="AA") {
      stop("the alphabet of pcm must be AA")
    }
    if(is(.ele, 'pfm')){
      .ele <- pfm2pcm(.ele)
    }
    .ele
  })
  pcmcopy<-list(length=length(pcms))
  pcmcopy[[1]]<-pcms[[1]]
  for(i in 2:length(pcms)){
    pcmcopy[[i]]<-pcms[[i]]
    pcmcopy<-UngappedAlignmentAA(pcmcopy, i, 
                                 threshold=threshold,
                                 minimalConsensus = minimalConsensus)
  }
  l<-unlist(lapply(pcmcopy,function(.ele) ncol(.ele@mat)))
  l<-max(l)-l
  for(i in 1:length(pcmcopy)){
    pcmcopy[[i]]<-addBlank(pcmcopy[[i]],l[i],TRUE)
  }
  pcmcopy
}