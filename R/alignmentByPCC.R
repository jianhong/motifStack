alignmentByPCC <- function(motifs, minimalConsensus=0, rcpostfix="(RC)", revcomp=rep(TRUE, length(motifs))){
  if(length(motifs)<2) stop("less than 2 motifs")
  if(length(revcomp)!=length(motifs)) stop("length of revcomp and motifs is not identical")
  lapply(motifs,function(.ele){
    if(!inherits(.ele, c("pfm", "psam"))) {
      stop("motifs must be a list of pfm or psam objects.")
    }
    if(.ele@alphabet!="DNA") stop("the alphabet of motif must be DNA")
  })
  motifcopy<-list(length=length(motifs))
  motifcopy[[1]]<-motifs[[1]]
  for(i in 2:length(motifs)){
    motifcopy[[i]]<-motifs[[i]]
    motifcopy<-UngappedPCCAlignment(motifcopy, i, minimalConsensus, rcpostfix, revcomp[i])
  }
  l<-unlist(lapply(motifcopy,function(.ele) ncol(.ele@mat)))
  l<-max(l)-l
  for(i in 1:length(motifcopy)){
    motifcopy[[i]]<-addBlank(motifcopy[[i]],l[i],TRUE)
  }
  motifcopy
}

UngappedPCCAlignment <- function(motifs, i, minimalConsensus=0, rcpostfix="(RC)", revcomp=TRUE){
  res<-getAlignedPCCWithoutGap(motifs[[i-1]], motifs[[i]], revcomp)
  if(res$max>=minimalConsensus){
    if(res$rev){
      motifs[[i]]<-matrixReverseComplement(motifs[[i]])
      motifs[[i]]@name<-paste(motifs[[i]]@name, rcpostfix, sep="")
    }
    if(res$offset>0){
      motifs[[i]]<-addBlank(motifs[[i]], res$offset, FALSE)
    }else{
      if(res[1]<0){
        motifs[1:(i-1)]<-lapply(motifs[1:(i-1)],function(.ele) addBlank(.ele, -1*res$offset, FALSE))
      }
    }
  }
  motifs
}

getAlignedPCCWithoutGap<-function(motif1, motif2, revcomp=TRUE){
  offset1<-getoffsetPosByPCC(motif1@mat, motif2@mat)
  if(revcomp){
    motif3<-matrixReverseComplement(motif2)
    offset2<-getoffsetPosByPCC(motif1@mat, motif3@mat)
  }else{
    offset2<-list(k=0,max=0)
  }
  offset<-0
  rev<-FALSE
  if((offset1$max < offset2$max)){
    rev<-TRUE
    max<-offset2$max
    offset<-offset2$k
  }else{
    max<-offset1$max
    offset<-offset1$k
  }
  list(offset=offset,rev=rev,max=max)    
}

#' @importFrom stats cor.test
getoffsetPosByPCC<-function(query, subject){
  if(!is(query, "matrix") || !is(subject, "matrix")){
    stop("query and subject must be numeric matrix")
  }
  score1<-rep(0, ncol(query))
  score2<-rep(0, ncol(subject))
  pvalue1<-rep(0, ncol(query))
  pvalue2<-rep(0, ncol(subject))
  for(i in seq.int(ncol(query))){
    J<-ifelse(ncol(query)-i>ncol(subject)-1,ncol(subject)-1,ncol(query)-i)
    corr <- cor.test(query[, i:(J+i)], subject[, 1:(J+1)], method="pearson")
    score1[i] <- J+1
    pvalue1[i] <- corr$p.value
  }
  for(i in 1:ncol(subject)){
    J<-ifelse(ncol(subject)-i>ncol(query)-1,ncol(query)-1,ncol(subject)-i)
    corr <- cor.test(query[, 1:(J+1)], subject[, i:(J+i)], method="pearson")
    score2[i] <- J+1
    pvalue2[i] <- corr$p.value
  }
  k1<-match(min(pvalue1),pvalue1)
  k2<-match(min(pvalue2),pvalue2)
  sc <- pvalue1[k1] > pvalue2[k2]
  k <- ifelse(sc, -1*(k2-1), (k1-1))
  max <- ifelse(sc, score2[k2], score1[k1])
  list(k=k, max=max)
}
