plotMotifLogoStack<-function(pfms, ...){
  n<-length(pfms)
  if(all(sapply(pfms, function(.ele) is(.ele, "psam")))){
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
    if(!is(.ele, "pfm")) stop("pfms must be a list of pfm objects.")
  })
  assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
  grid.newpage()
  ht <- 1/n
  y0 <- .5 * ht
  for(i in seq.int(n)){
    pushViewport(viewport(y=y0, height=ht))
    plotMotifLogo(pfms[[i]], motifName=pfms[[i]]@name, 
                  p=pfms[[i]]@background, colset = pfms[[i]]@color,
                  xlab=NA, newpage=FALSE, 
                  margins=c(1.5, 4.1, 1.1, .1), ...)
    popViewport()
    y0 <- y0 + ht
  }
  rm(list="tmp_motifStack_symbolsCache", envir=.globals)
  return()
}
