#' align DNA motifs
#' 
#' align DNA motifs for plotting motifs stack
#' 
#' 
#' @param pfms a list of position frequency matrices, pfms must be a list of
#' class pfm or psam
#' @param threshold information content cutoff threshold for useful postions
#' @param minimalConsensus minimal length of consensus for alignment
#' @param rcpostfix the postfix for reverse complements
#' @param revcomp a logical vector to indicates whether the reverse complemet
#' should be involved into alignment
#' @return a list of aligned motifs
#' @export
#' @examples
#' 
#' pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
#' motifs<-lapply(pcms,pcm2pfm)
#' motifs<-DNAmotifAlignment(motifs)
#' 
DNAmotifAlignment<-function(pfms, threshold=0.4, minimalConsensus=0, 
                            rcpostfix="(RC)", 
                            revcomp=rep(TRUE, length(pfms))){
  if(length(pfms)<2) stop("less than 2 motifs")
  if(length(revcomp)!=length(pfms)){
    stop("length of revcomp and pfms is not identical")
  }
  if(all(sapply(pfms, function(.ele) is(.ele, "psam")))){
    return(alignmentByPCC(pfms, minimalConsensus, rcpostfix, revcomp))
  }
  lapply(pfms,function(.ele){
    if(!inherits(.ele, c("pfm", "pcm"))) {
      stop("pfms must be a list of pfm objects.")
    }
    if(!.ele@alphabet %in% c("DNA", 'RNA')) {
      stop("the alphabet of pfm must be DNA or RNA")
    }
  })
  alphab <- vapply(pfms, function(.ele) .ele@alphabet=="RNA", FUN.VALUE = TRUE)
  if(any(alphab & revcomp)){
    message("The input is RNA motif. 
            The reverse complementation alignment will be turn off.
            To avoid this message, please set revcomp=rep(FALSE, length(pfms))")
    revcomp[alphab] <- FALSE
  }
  pfmcopy<-list(length=length(pfms))
  pfmcopy[[1]]<-pfms[[1]]
  for(i in 2:length(pfms)){
    pfmcopy[[i]]<-pfms[[i]]
    pfmcopy<-UngappedAlignment(pfmcopy, i, threshold, 
                               minimalConsensus, rcpostfix, revcomp[i])
  }
  l<-unlist(lapply(pfmcopy,function(.ele) ncol(.ele@mat)))
  l<-max(l)-l
  for(i in 1:length(pfmcopy)){
    pfmcopy[[i]]<-addBlank(pfmcopy[[i]],l[i],TRUE)
  }
  pfmcopy
}
