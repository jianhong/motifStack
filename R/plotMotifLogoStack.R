plotMotifLogoStack<-function(pfms, ...){
  n<-length(pfms)
  lapply(pfms,function(.ele){
    if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
  })
  opar<-par(mfrow=c(n,1),mar=c(3.5,3.5,1.5,0.5))
  assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
  for(i in 1:(n-1)){
    plot(pfms[[n-i+1]],xlab=NA, ...)
  }
  plot(pfms[[1]], ...)
  rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
  par(opar)
}
