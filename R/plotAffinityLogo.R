plotAffinityLogo <- function(psam, motifName, font="Helvetica-Bold",
                             colset=c("#00811B","#2000C7","#FFB32C","#D00001"),
                             alpha=0.5, newpage=TRUE, fontsize=12){
  if(class(psam)=="data.frame"){
    psam <- as.matrix(psam)
  }else{
    if(class(psam) == "psam"){
      psam <- psam@mat
    }
  }
  if (class(psam) != "matrix"){
    stop("psam must be of class matrix or data.frame")
  }
  if(class(colset)!="character")
    stop("colset must be of class character")
  if (length(colset)!=nrow(psam))
    stop(paste("colset length and psam row numbers different",length(colset),"!=",nrow(psam)))
  rname<-rownames(psam)
  if(is.null(rname))
    stop("psam rowname is empty")
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
    if(plot.new) plot.new()
    text(.5,.5,"Can not locate ghostscript.\nFor windows, please install ghostscript first and then try:\nSys.setenv(R_GSCMD=\"\\\"C:\\\\Path\\\\to\\\\gsbin\\\\gswin32c.exe\\\"\")")
  }else{
    npos<-ncol(psam)
    ncha<-nrow(psam)
    key<-paste("x", ncha, font, paste(colset, collapse=""), paste(rname, collapse=""), sep="_")
    symbolsCache <- if(exists("tmp_motifStack_symbolsCache", envir=.globals)) get("tmp_motifStack_symbolsCache", envir=.globals) else list()
    if(!is.null(symbolsCache[[key]])){
      symbols<-symbolsCache[[key]]
    } else {
      symbols<-coloredSymbols(ncha, font, colset, rname, fontsize)
      symbolsCache[[key]]<-symbols
      assign("tmp_motifStack_symbolsCache", symbolsCache, envir=.globals)
    }
    #calculate postion of each symbol and plot
    if(newpage) grid.newpage()
    
    ddG <- log(psam+1e-2)
    ddG.mu<-colMeans(ddG)
    ddG.height <- t(t(ddG) - ddG.mu)
    ddG.h.pos <- ddG.h.neg <- ddG.height
    ddG.h.pos[ddG.h.pos<0] <- 0
    ddG.h.neg[ddG.h.neg>0] <- 0
    ie <- max(abs(min(colSums(ddG.h.neg))), max(colSums(ddG.h.pos)))
    ddG.height <- ddG.height/(2*ie)
    dw<-1/npos
    x.pos<-0
    grid.lines(y=unit(c(.5, .5), "npc"), gp=gpar(lty=3))
    for(j in seq.int(npos)){
      heights<-ddG.height[,j]
      ## less than 0
      id<-order(heights, decreasing = TRUE)
      y.pos<-0.5
      for(i in seq.int(ncha)){
        h<-heights[id[i]]
        if(h<0) {
          y.pos<-y.pos+h
          grid.draw(pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,-h,just=c(0,0),distort=TRUE, gp=gpar(alpha=alpha)))
        }
      }
      ## greater than 0
      id<-order(heights, decreasing = FALSE)
      y.pos<-0.5
      for(i in seq.int(ncha)){
        h<-heights[id[i]]
        if(h>0) {
          grid.draw(pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h,just=c(0,0),distort=TRUE))
          y.pos<-y.pos+h
        }
      }
      x.pos<-x.pos+dw
    }
  }
}

