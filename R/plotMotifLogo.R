plotMotifLogo<-function(pfm, motifName, p=rep(0.25, 4), font="Helvetica-Bold", fontface="bold", 
                        colset=c("#00811B","#2000C7","#FFB32C","#D00001"), 
                        xaxis=TRUE,yaxis=TRUE,xlab="position",ylab="bits",
                        xlcex=1.2, ylcex=1.2, ncex=1.2, ic.scale=TRUE,
                        newpage=TRUE, margins=c(4.1, 4.1, 2.1, .1), draw=TRUE){
  markers <- NULL
  if (class(pfm) == "data.frame"){
    pfm <- as.matrix(pfm)
  }else{
    if(class(pfm) == "pcm"){
      pfm <- pcm2pfm(pfm)
    }
    if(class(pfm) == "pfm"){
      markers <- pfm@markers
      if(missing(motifName)) motifName = pfm@name
      p=pfm@background[rownames(pfm@mat)]
      colset=pfm@color[rownames(pfm@mat)]
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
    symbols<-coloredSymbols(ncha, font, colset, rname, fontface=fontface)
    symbolsCache[[key]]<-symbols
    assign("tmp_motifStack_symbolsCache", symbolsCache, envir=.globals)
  }
  #calculate postion of each symbol and plot
  plot <- gList()
  ic<-getIC(pfm, p)
  if(ic.scale){
    ie<-getIE(pfm)
    ie <- max(c(ie, ic))
  }else{
    ie <- 1
  }
  dw<-1/npos
  x.pos<-0
  y.poss <- numeric(length=npos)
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
      if(h>0 && ic[j]>0) plot <- gList(plot, pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h, just=c(0,0),distort=TRUE))
      y.pos<-y.pos+h
    }
    y.poss[j] <- y.pos
    x.pos<-x.pos+dw
  }
  if(length(markers)>0){
    plot <- gList(plot, plotMarkers(markers, dw, y.poss))
  }
  if(xaxis) plot <- gList(plot, plotXaxis(pfm, p))
  if(yaxis) plot <- gList(plot, plotYaxis(ie))
  plot <- gTree(children=plot, vp= plotViewport(margins = margins))
  if(!is.na(xlab)) plot <- gList(plot, textGrob(xlab, y=unit(1, units = "lines"), gp=gpar(cex=xlcex), name="xlab"))
  if(!is.na(ylab)) plot <- gList(plot, textGrob(ylab, x=unit(1, units = "lines"), gp=gpar(cex=ylcex), rot=90, name="ylab"))
  if(!missing(motifName)) plot <- gList(plot, textGrob(motifName,y=unit(1, "npc")-unit(.5, units = "lines"), gp=gpar(cex=ncex), name="motifName"))
  if(draw){
    if(newpage) grid.newpage()
    suppressWarnings(grid.draw(plot))
  }else{
    plot
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
  grob <- xaxisGrob(at=at,label=label, gp=gpar(lwd=1, lex=1, lineheight=1))
  grob$children$labels$y <- unit(-1.1, "lines")
  grob
}


plotYaxis<-function(ymax){
  ie <- ymax
  majorat<-seq(0,floor(ie),length.out = 5)
  majorat <- majorat/ie
  majorlab<-seq(0,floor(ie),length.out = 5)
  majorY <- yaxisGrob(at=majorat,label=majorlab,
                      name = "majorY",
                      gp=gpar(lwd=1, lex=1, lineheight=1))
  YgList <- gList(majorY)
  if(convertUnit(unit(1, "npc"), unitTo = "lines", valueOnly = TRUE) > 10){
    minorat<-seq(0, floor(ie), length.out = 25)
    minorat <- seq(0, ie, by=diff(minorat)[1])
    minorat <- minorat/ie
    minorat<-minorat[!(minorat %in% majorat)]
    minorY <- yaxisGrob(at=minorat,label=FALSE,name="minorY",
                        gp=gpar(lwd=1, lineheight=.5))
    at <- c(min(c(minorat, majorat)), max(c(minorat, majorat)))
    terminY <- yaxisGrob(at=at, label=FALSE, name="endY",
                         gp=gpar(lwd=1, lineheight=.5))
    YgList <- gList(YgList, minorY, terminY)
  }
  YgList
}


###############################################################################
######## plot motif logo without plot.new
######## to be used to create a better view of stack, eg. radial sty,
###############################################################################
plotMotifLogoA<-function(pfm, font="Helvetica-Bold", fontface="bold", ic.scale=TRUE, draw=TRUE){
  if (class(pfm) != "pfm"){
    stop("pfms must be a pfm object")
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
    symbols<-coloredSymbols(ncha, font, pfm@color, rname, fontface=fontface)
    symbolsCache[[key]]<-symbols
    assign("tmp_motifStack_symbolsCache", symbolsCache, envir=.globals)
  }
  
  #calculate postion of each symbol and plot   
  plot <- gList()
  ic<-getIC(pfm)
  ie<-getIE(pfm@mat)
  ie <- max(c(ie, ic))
  dw<-1/npos
  x.pos<-0
  y.poss <- numeric(length=npos)
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
      if(h>0 && ic[j]>0) plot <- gList(plot, pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h,just=c(0,0),distort=TRUE))
      y.pos<-y.pos+h
    }
    y.poss[j] <- y.pos
    x.pos<-x.pos+dw
  }
  markers <- pfm@markers
  if(length(markers)>0){
    plot <- gList(plot, plotMarkers(markers, dw, y.poss))
  }
  if(draw){
    suppressWarnings(grid.draw(plot))
  }else{
    plot
  }
}

plotMarkers <- function(markers, dw, h, lo=NULL){
  do.call(gList, lapply(markers, function(m){
    switch(m@type,
           "rect"={
             pos <- mapply(seq, m@start, m@stop, SIMPLIFY = FALSE)
             if(length(lo)>0 && length(lo)==length(h)){
               height <- sapply(pos, function(.ele) max(h[.ele])-min(lo[.ele]))
               y <- sapply(pos, function(.ele) (min(lo[.ele])+max(h[.ele]))/2)
             }else{
               height <- sapply(pos, function(.ele) max(h[.ele]))
               y <- height/2
             }
             rectGrob(x= (m@start + m@stop-1)*dw/2,
                      y = y,
                      width = dw*(m@stop - m@start + 1),
                      height = height, 
                      gp = m@gp)
             },
           "text"={
             pos <- mapply(seq, m@start, m@stop, SIMPLIFY = FALSE)
             label <- rep(m@label, lengths(pos))
             pos <- unlist(pos)
             textGrob(label=label,
                      x = (pos-.5)*dw, 
                      y = h[pos],
                      vjust = -0.1,
                      gp = m@gp)
             },
           "line"={
             pos <- mapply(seq, m@start, m@stop, SIMPLIFY = FALSE)
             height <- sapply(pos, function(.ele) max(h[.ele]))
             linesGrob(x = as.numeric(rbind(m@start-1, m@stop)*dw), 
                       y = height,
                       gp = m@gp)
             })
  }))
}

################## for ggplot2 ###################
motifGrob <- function(pfm, x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                      width = unit(1, "npc"), height = unit(1, "npc"),
                      angle = 0, ic.scale=TRUE, default.units = "native", name=NULL, 
                      gp = gpar(fontfamily="Helvetica-Bold",
                                fontface="bold")){
  vp <- viewport(x=x, y=y, width = width, height = height, 
                 default.units = default.units, angle = angle, name = "motifVP")
  fontfamily <- ifelse(length(gp$fontfamily)>0, gp$fontfamily[1], "Helvetica-Bold")
  fontface <- ifelse(length(gp$fontface)>0, gp$fontface[1], "bold")
  gTree(children=plotMotifLogoA(pfm, font = fontfamily, fontface = fontface, ic.scale = ic.scale, draw=FALSE),
        name=name, vp=vp, cl="motif")
}

geom_motif <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       ...,
                       ic.scale=TRUE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomMotif,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      ic.scale = ic.scale,
      ...
    )
  )
}

GeomMotif <- 
  ggproto("GeoMotif", Geom, 
          required_aes = c("xmin", "ymin", "xmax", "ymax", "motif"),
          default_aes = aes(angle=0, fontfamily="Helvetica-Bold", 
                            fontface="bold"),
          draw_panel = function(data, panel_params, coord, ic.scale=TRUE){
            coords <- coord$transform(data, panel_params)
            n <- nrow(coords)
            do.call(gList, lapply(seq.int(n), function(.id){
              motifGrob(pfm=coords$motif[[.id]],
                        x=mean(c(coords$xmin[.id], coords$xmax[.id])),
                        y=mean(c(coords$ymin[.id], coords$ymax[.id])),
                        width=coords$xmax[.id] - coords$xmin[.id], 
                        height=coords$ymax[.id] - coords$ymin[.id],
                        angle=coords$angle[.id], ic.scale=ic.scale,
                        default.units = "native",
                        gp = gpar(fontfamily=coords$fontfamily[.id],
                                  fontface=coords$fontface[.id]),
                        name = paste0("geom_motif_", .id))
            }))
          })

