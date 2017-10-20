plotMotifLogo<-function(pfm, motifName, p=rep(0.25, 4), font="Helvetica-Bold", 
                        colset=c("#00811B","#2000C7","#FFB32C","#D00001"), 
                        xaxis=TRUE,yaxis=TRUE,xlab="position",ylab="bits",
                        xlcex=1.2, ylcex=1.2, ncex=1.2, ic.scale=TRUE, 
                        fontsize=12){
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
  # get ghostscript path
  gscmd <- Sys.getenv("R_GSCMD")
  if(is.null(gscmd) || !nzchar(gscmd)) {
    gscmd <- switch(.Platform$OS.type,
                    unix = "gs",
                    windows = {
                      poss <- Sys.which(c("gswin64c.exe",
                                          "gswin32c.exe"))
                      poss <- poss[nzchar(poss)]
                      if (length(poss)) gscmd <- poss
                    })
  }
  if(is.null(gscmd) || !nzchar(gscmd)){
    plot.new()
    text(.5,.5,"Can not locate ghostscript.\nFor windows, please install ghostscript first and then try:\nSys.setenv(R_GSCMD=\"\\\"C:\\\\Path\\\\to\\\\gsbin\\\\gswin32c.exe\\\"\")")
  }else{
    npos<-ncol(pfm)
    ncha<-nrow(pfm)
    key<-paste("x", ncha, font, paste(colset, collapse=""), paste(rname, collapse=""), sep="_")
    symbolsCache <- if(exists("tmp_motifStack_symbolsCache", where=".GlobalEnv")) get("tmp_motifStack_symbolsCache", pos=".GlobalEnv") else list()
    if(!is.null(symbolsCache[[key]])) symbols<-symbolsCache[[key]]
    else {
      symbols<-coloredSymbols(ncha, font, colset, rname, fontsize)
      symbolsCache[[key]]<-symbols
      assign("tmp_motifStack_symbolsCache", symbolsCache, pos=".GlobalEnv")
    }
    #calculate postion of each symbol and plot
    plot.new()
    
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
        if(h>0 && ic[j]>0) picture(symbols[[id[i]]],x.pos,y.pos,x.pos+dw,y.pos+h)
        y.pos<-y.pos+h
      }
      x.pos<-x.pos+dw
    }
    if(xaxis) plotXaxis(pfm, p)
    if(yaxis) plotYaxis(ie)
    if(!is.na(xlab)) mtext(xlab,1,line=2,cex=xlcex)
    if(!is.na(ylab)) mtext(ylab,2,line=2,cex=ylcex)
    if(!missing(motifName)) mtext(motifName,3,line=0,cex=ncex)
  }
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
  axis(1,at=at,labels=label)
}

plotYaxis<-function(ymax){
  ie <- ymax
  majorat<-seq(0,floor(ie)/ceiling(ie),by=1/ceiling(ie))
  majorlab<-0:floor(ie)
  axis(2,at=majorat,labels=majorlab, lwd=0, lwd.ticks=1)
  minorat<-seq(0,ie,by=1/(ceiling(ie)*5))
  minorat<-minorat[!(minorat %in% majorat)]
  axis(2,at=minorat,tcl=par("tcl")*0.5,labels=FALSE,lwd=0,lwd.ticks=1)
  at <- c(min(c(minorat, majorat)), max(c(minorat, majorat)))
  axis(2, at=at, lwd=1, lwd.ticks=0, labels=FALSE)
}

###############################################################################
######## plot motif logo without plot.new
######## to be used to create a better view of stack, eg. radial sty,
###############################################################################
plotMotifLogoA<-function(pfm, font="Helvetica-Bold", ic.scale=TRUE, 
                         fontsize=12){
  if (class(pfm) != "pfm"){
    stop("pfms must be a list of class pfm")
  }
  rname<-rownames(pfm@mat)
  if(is.null(rname))
    stop("pfm rowname is empty")
  # get ghostscript path
  gscmd <- Sys.getenv("R_GSCMD")
  npos<-ncol(pfm@mat)
  ncha<-nrow(pfm@mat) 
  key<-paste("x", ncha, font, paste(pfm@color, collapse=""), paste(rname, collapse=""), sep="_")
  symbolsCache <- if(exists("tmp_motifStack_symbolsCache", where=".GlobalEnv")) get("tmp_motifStack_symbolsCache", pos=".GlobalEnv") else list()
  if(!is.null(symbolsCache[[key]])) symbols<-symbolsCache[[key]]
  else {
    symbols<-coloredSymbols(ncha, font, pfm@color, rname, fontsize)
    symbolsCache[[key]]<-symbols
    assign("tmp_motifStack_symbolsCache", symbolsCache, pos=".GlobalEnv")
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
      if(h>0 && ic[j]>0) grid.draw(pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h,just=c(0,0),distort=TRUE))
      y.pos<-y.pos+h
    }
    x.pos<-x.pos+dw
  }
}
