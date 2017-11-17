plotMotifLogoStack<-function(pfms, ...){
  n<-length(pfms)
  if(all(sapply(pfms, class)=="psam")){
    grid.newpage()
    ht <- 1/n
    y0 <- .5 * ht
    for(i in seq.int(n)){
      pushViewport(viewport(y=y0, height=ht))
      ht.title <- convertUnit(unit(1.5, "lines"), unitTo = "npc", valueOnly = TRUE)
      ht.body <- 1 - ht.title
      grid.text(label=pfms[[i]]@name, y=1-.5*ht.title)
      pushViewport(viewport(y=ht.body*.5, height=ht.body))
      plotAffinityLogo(pfms[[i]], newpage=FALSE)
      popViewport()
      popViewport()
      y0 <- y0 + ht
    }
    return()
  }
  lapply(pfms,function(.ele){
    if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
  })
  opar<-par(mfrow=c(n,1),mar=c(3.5,3.5,1.5,0.5))
  assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
  for(i in 1:(n-1)){
    plot(pfms[[n-i+1]],xlab=NA, ...)
  }
  plot(pfms[[1]], ...)
  rm(list="tmp_motifStack_symbolsCache", envir=.globals)
  par(opar)
}
