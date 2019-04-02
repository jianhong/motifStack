plotMotifLogo<-function(pfm, motifName, p=rep(0.25, 4), font="Helvetica-Bold", 
                        colset=c("#00811B","#2000C7","#FFB32C","#D00001"), 
                        xaxis=TRUE,yaxis=TRUE,xlab="position",ylab="bits",
                        xlcex=1.2, ylcex=1.2, ncex=1.2, ic.scale=TRUE,
                        newpage=TRUE, margins=c(4.1, 4.1, 2.1, .1)){
  if (class(pfm) == "data.frame"){
    pfm <- as.matrix(pfm)
  }else{
    if(class(pfm) == "pcm"){
      pfm <- pcm2pfm(pfm)
    }
    if(class(pfm) == "pfm"){
      pfm <- pfm@mat
    }
  } 
  if (class(pfm) != "matrix"){
    stop("pfm must be of class matrix or data.frame")
  }
  if(length(p)<nrow(pfm)){
    warning("background length is shorter than number of rows. Will use default background")
    p <- rep(1/nrow(pfm), nrow(pfm))
  }
  if (any(abs(1 - apply(pfm,2,sum)) > 0.01))
    stop("Columns of pfm must add up to 1.0")
  if(class(colset)!="character")
    stop("colset must be of class character")
  if (length(colset)!=nrow(pfm))
    stop(paste("colset length and pfm row numbers different",length(colset),"!=",nrow(pfm)))
  rname<-rownames(pfm)
  if(is.null(rname))
    stop("pfm rowname is empty")
  
  npos<-ncol(pfm)
  ncha<-nrow(pfm)
  key<-paste("x", ncha, font, paste(colset, collapse=""), paste(rname, collapse=""), sep="_")
  symbolsCache <- if(exists("tmp_motifStack_symbolsCache", envir=.globals)) get("tmp_motifStack_symbolsCache", envir=.globals) else list()
  if(!is.null(symbolsCache[[key]])) symbols<-symbolsCache[[key]]
  else {
    symbols<-coloredSymbols(ncha, font, colset, rname)
    symbolsCache[[key]]<-symbols
    assign("tmp_motifStack_symbolsCache", symbolsCache, envir=.globals)
  }
  #calculate postion of each symbol and plot
  if(newpage) grid.newpage()
  
  vp <- plotViewport(margins = margins)
  pushViewport(vp)
  ic<-getIC(pfm, p)
  if(ic.scale){
    ie<-getIE(pfm)
    ie <- max(c(ie, ic))
  }else{
    ie <- 1
  }
  dw<-1/npos
  x.pos<-0
  for(j in 1:npos){
    column<-pfm[,j]
    if(ic.scale){
      heights<-column*ic[j]/ie
    }else{
      heights <- column
    }
    id<-order(heights)
    
    y.pos<-0
    for(i in 1:ncha){
      h<-heights[id[i]]
      if(h>0 && ic[j]>0) suppressWarnings(grid.draw(pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h, just=c(0,0),distort=TRUE)))
      y.pos<-y.pos+h
    }
    x.pos<-x.pos+dw
  }
  if(xaxis) plotXaxis(pfm, p)
  if(yaxis) plotYaxis(ie)
  popViewport()
  if(!is.na(xlab)) grid.text(xlab, y=unit(1, units = "lines"), gp=gpar(cex=xlcex))
  if(!is.na(ylab)) grid.text(ylab, x=unit(1, units = "lines"), gp=gpar(cex=ylcex), rot=90)
  if(!missing(motifName)) grid.text(motifName,y=unit(1, "npc")-unit(1, units = "lines"), gp=gpar(cex=ncex))
}

plotXaxis<-function(pfm, p=rep(0.25, 4)){
  npos<-ncol(pfm)
  ic<-getIC(pfm,p)
  dw<-1/npos
  at<-c()
  label<-c()
  j=1
  for(i in 1:npos){
    if(ic[i]>0){
      at<-c(at,dw*(i-0.5))
      label<-c(label,j)
      j<-j+1
    }
  }
  grob <- xaxisGrob(at=at,label=label, gp=gpar(lwd=1, lex=1, lineheight=1))
  grob$children$labels$y <- unit(-1.1, "lines")
  grid.draw(grob)
}

plotYaxis<-function(ymax){
  ie <- ymax
  majorat<-seq(0,floor(ie),length.out = 5)
  majorat <- majorat/ie
  majorlab<-seq(0,floor(ie),length.out = 5)
  grid.yaxis(at=majorat,label=majorlab, gp=gpar(lwd=1, lex=1, lineheight=1))
  
  if(convertUnit(unit(1, "npc"), unitTo = "lines", valueOnly = TRUE) > 5){
    minorat<-seq(0, floor(ie), length.out = 25)
    minorat <- seq(0, ie, by=diff(minorat)[1])
    minorat <- minorat/ie
    minorat<-minorat[!(minorat %in% majorat)]
    grid.yaxis(at=minorat,label=FALSE,gp=gpar(lwd=1, lineheight=.5))
    at <- c(min(c(minorat, majorat)), max(c(minorat, majorat)))
    grid.yaxis(at=at, label=FALSE, gp=gpar(lwd=1, lineheight=.5))
  }
}

###############################################################################
######## plot motif logo without plot.new
######## to be used to create a better view of stack, eg. radial sty,
###############################################################################
plotMotifLogoA<-function(pfm, font="Helvetica-Bold", ic.scale=TRUE){
  if (class(pfm) != "pfm"){
    stop("pfms must be a list of class pfm")
  }
  rname<-rownames(pfm@mat)
  if(is.null(rname))
    stop("pfm rowname is empty")
  npos<-ncol(pfm@mat)
  ncha<-nrow(pfm@mat) 
  key<-paste("x", ncha, font, paste(pfm@color, collapse=""), paste(rname, collapse=""), sep="_")
  symbolsCache <- if(exists("tmp_motifStack_symbolsCache", envir=.globals)) get("tmp_motifStack_symbolsCache", envir=.globals) else list()
  if(!is.null(symbolsCache[[key]])) symbols<-symbolsCache[[key]]
  else {
    symbols<-coloredSymbols(ncha, font, pfm@color, rname)
    symbolsCache[[key]]<-symbols
    assign("tmp_motifStack_symbolsCache", symbolsCache, envir=.globals)
  }
  
  #calculate postion of each symbol and plot   
  ic<-getIC(pfm)
  ie<-getIE(pfm@mat)
  ie <- max(c(ie, ic))
  dw<-1/npos
  x.pos<-0
  for(j in 1:npos){
    column<-pfm@mat[,j]
    if(ic.scale){
      heights<-column*ic[j]/ie
    }else{
      heights <- column
    }
    id<-order(heights)
    
    y.pos<-0
    for(i in 1:ncha){
      h<-heights[id[i]]
      if(h>0 && ic[j]>0) suppressWarnings(grid.draw(pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h,just=c(0,0),distort=TRUE)))
      y.pos<-y.pos+h
    }
    x.pos<-x.pos+dw
  }
}
