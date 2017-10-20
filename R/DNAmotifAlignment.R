DNAmotifAlignment<-function(pfms, threshold=0.4, minimalConsensus=0, rcpostfix="(RC)", revcomp=rep(TRUE, length(pfms))){
  if(length(pfms)<2) stop("less than 2 motifs")
  if(length(revcomp)!=length(pfms)) stop("length of revcomp and pfms is not identical")
  if(all(sapply(pfms, class)=="psam")){
    return(alignmentByPCC(pfms, minimalConsensus, rcpostfix, revcomp))
  }
  lapply(pfms,function(.ele){
    if(class(.ele)!="pfm" && class(.ele)!="pcm") stop("pfms must be a list of class pfm")
    if(.ele@alphabet!="DNA") stop("the alphabet of pfm must be DNA")
  })
  pfmcopy<-list(length=length(pfms))
  pfmcopy[[1]]<-pfms[[1]]
  for(i in 2:length(pfms)){
    pfmcopy[[i]]<-pfms[[i]]
    pfmcopy<-UngappedAlignment(pfmcopy, i, threshold, minimalConsensus, rcpostfix, revcomp[i])
  }
  l<-unlist(lapply(pfmcopy,function(.ele) ncol(.ele@mat)))
  l<-max(l)-l
  for(i in 1:length(pfmcopy)){
    pfmcopy[[i]]<-addBlank(pfmcopy[[i]],l[i],TRUE)
  }
  pfmcopy
}
