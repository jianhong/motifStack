plotAffinityLogo <- function(psam, motifName, font="Helvetica-Bold", fontface="bold",
                             colset=c("#00811B","#2000C7","#FFB32C","#D00001"),
                             alpha=0.5, newpage=TRUE, draw=TRUE){
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
  npos<-ncol(psam)
  ncha<-nrow(psam)
  key<-paste("x", ncha, font, paste(colset, collapse=""), paste(rname, collapse=""), sep="_")
  symbolsCache <- if(exists("tmp_motifStack_symbolsCache", envir=.globals)) get("tmp_motifStack_symbolsCache", envir=.globals) else list()
  symbols<-coloredSymbols(ncha, font, colset, rname, alpha, fontface=fontface)
  symbolsCache[[key]]<-symbols
  assign("tmp_motifStack_symbolsCache", symbolsCache, envir=.globals)
  
  #calculate postion of each symbol and plot
  plot <- gList()
  
  
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
  plot <- gList(plot, linesGrob(y=unit(c(.5, .5), "npc"), gp=gpar(lty=3)))
  for(j in seq.int(npos)){
    heights<-ddG.height[,j]
    ## less than 0
    id<-order(heights, decreasing = TRUE)
    y.pos<-0.5
    for(i in seq.int(ncha)){
      h<-heights[id[i]]
      if(h<0) {
        y.pos<-y.pos+h
        plot <- gList(plot, pictureGrob(symbols[[paste0(id[i], "_", alpha)]],x.pos,y.pos,dw,-h,just=c(0,0),distort=TRUE))
      }
    }
    ## greater than 0
    id<-order(heights, decreasing = FALSE)
    y.pos<-0.5
    for(i in seq.int(ncha)){
      h<-heights[id[i]]
      if(h>0) {
        plot <- gList(plot, pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h,just=c(0,0),distort=TRUE))
        y.pos<-y.pos+h
      }
    }
    x.pos<-x.pos+dw
  }
  
  if(draw){
    if(newpage) grid.newpage()
    suppressWarnings(grid.draw(plot))
  }else{
    plot
  }
}

