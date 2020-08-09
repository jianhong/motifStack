#' plot sequence logo
#' 
#' plot amino acid or DNA sequence logo
#' 
#' 
#' @param pfm a position frequency matrices
#' @param motifName motif name
#' @param p background possibility
#' @param font font of logo
#' @param fontface fontface of logo
#' @param colset color setting for each logo letter
#' @param xaxis draw x-axis or not
#' @param yaxis draw y-axis or not
#' @param xlab x-label, do nothing if set xlab as NA
#' @param ylab y-label, do nothing if set ylab as NA
#' @param xlcex cex value for x-label
#' @param ylcex cex value for y-label
#' @param ncex cex value for motif name
#' @param ic.scale logical If TRUE, the height of each column is proportional
#' to its information content. Otherwise, all columns have the same height.
#' @param newpage logical If TRUE, plot it in a new page.
#' @param margins A numeric vector interpreted in the same way as par(mar) in
#' base graphics.
#' @param draw Vector (logical(1)). TRUE to plot. FALSE, return a gList
#' @return none
#' @export
#' @importFrom grid gList gTree plotViewport textGrob unit gpar grid.newpage
#' grid.draw gList grid.draw
#' @examples
#' 
#' pcm<-matrix(runif(40,0,100),nrow=4,ncol=10)
#' pfm<-pcm2pfm(pcm)
#' rownames(pfm)<-c("A","C","G","T")
#' plotMotifLogo(pfm)
#' 
plotMotifLogo<-function(pfm, motifName, p=rep(0.25, 4), 
                        font="Helvetica-Bold", fontface="bold", 
                        colset=c("#00811B","#2000C7","#FFB32C","#D00001"), 
                        xaxis=TRUE,yaxis=TRUE,xlab="position",ylab="bits",
                        xlcex=1.2, ylcex=1.2, ncex=1.2, ic.scale=TRUE,
                        newpage=TRUE, margins=c(4.1, 4.1, 2.1, .1), draw=TRUE){
  markers <- NULL
  if (is(pfm, "data.frame")){
    pfm <- as.matrix(pfm)
  }else{
    if(is(pfm, "pcm")){
      pfm <- pcm2pfm(pfm)
    }
    if(is(pfm, "pfm")){
      markers <- pfm@markers
      if(missing(motifName)) motifName = pfm@name
      p=pfm@background[rownames(pfm@mat)]
      p <- p/sum(p)
      colset=pfm@color[rownames(pfm@mat)]
      pfm <- pfm@mat
    }
  } 
  if (!is(pfm, "matrix")){
    stop("pfm must be a matrix")
  }
  if(length(p)<nrow(pfm)){
    warning("background length is shorter than number of rows. Will use default background")
    p <- rep(1/nrow(pfm), nrow(pfm))
  }
  if (any(abs(1 - apply(pfm,2,sum)) > 0.01))
    stop("Columns of pfm must add up to 1.0")
  if(!is(colset, "character"))
    stop("colset must be character")
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
  pictureGrob <- get("pictureGrob", envir = .globals)
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
      if(h>0 && ic[j]>0) plot <- 
          gList(plot, 
                pictureGrob(picture=symbols[[id[i]]],
                            x=x.pos,
                            y=y.pos,
                            width = dw,
                            height = h, 
                            just=c(0,0),
                            distort=TRUE))
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



#' plot x-axis
#' 
#' plot x-axis for the sequence logo
#' 
#' 
#' @param pfm position frequency matrices
#' @param p background possibility
#' @return none
#' @importFrom grid xaxisGrob gpar unit 
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




#' plot y-axis
#' 
#' plot y-axis for the sequence logo
#' 
#' 
#' @param ymax max value of y axix
#' @return none
#' @importFrom grid yaxisGrob gpar unit gList convertUnit
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


#' plot sequence logo without plot.new
#' 
#' plot amino acid or DNA sequence logo in a given canvas
#' 
#' 
#' @param pfm an object of pfm
#' @param font font of logo
#' @param fontface fontface of logo
#' @param ic.scale logical If TRUE, the height of each column is proportional
#' to its information content. Otherwise, all columns have the same height.
#' @param draw Vector (logical(1)). TRUE to plot. FALSE, return a gList
#' @return none
#' @export
#' @examples
#' 
#' pcm<-matrix(runif(40,0,100),nrow=4,ncol=10)
#' pfm<-pcm2pfm(pcm)
#' rownames(pfm)<-c("A","C","G","T")
#' motif <- new("pfm", mat=pfm, name="bin_SOLEXA")
#' plotMotifLogoA(motif)
#' 
plotMotifLogoA<-function(pfm, font="Helvetica-Bold", fontface="bold", ic.scale=TRUE, draw=TRUE){
  if (is(pfm, "pcm")){
    pfm <- pcm2pfm(pfm)
  }
  if (!is(pfm, "pfm")){
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
  pictureGrob <- get("pictureGrob", envir = .globals)
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
      if(h>0 && ic[j]>0) plot <- 
          gList(plot, 
                pictureGrob(picture = symbols[[id[i]]],
                            x = x.pos,
                            y = y.pos,
                            width = dw,
                            height = h,
                            just=c(0,0),
                            distort=TRUE))
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

#' @importFrom grid gList rectGrob textGrob linesGrob
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
             res <- rectGrob(x= (m@start + m@stop-1)*dw/2,
                             y = y,
                             width = dw*(m@stop - m@start + 1),
                             height = height, 
                             gp = m@gp)
             if(any(nchar(m@label)>0)){
               tG <- textGrob(label=m@label,
                              x = (m@start+m@stop-1)*dw/2,
                              y = height,
                              vjust = -.5,
                              gp = m@gp)
               res <- gList(res, tG)
             }
             res
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
             res <- segmentsGrob(x0 = (m@start-1)*dw, 
                                 x1 = m@stop*dw,
                                 y0 = height,
                                 y1 = height,
                                 gp = m@gp)
             if(any(nchar(m@label)>0)){
               tG <- textGrob(label=m@label,
                              x = (m@start+m@stop-1)*dw/2,
                              y = height,
                              vjust = -.5,
                              gp = m@gp)
               res <- gList(res, tG)
             }
             res
             })
  }))
}

################## for ggplot2 ###################


#' Motif Grob
#' 
#' This function create a motif grob.
#' 
#' 
#' @param pfm an object of pfm
#' @param x A numeric vector or unit object specifying x-values.
#' @param y A numeric vector or unit object specifying y-values.
#' @param width A numeric vector or unit object specifying width.
#' @param height A numeric vector or unit object specifying height.
#' @param angle A numeric value indicating the angle of rotation of the motif.
#' Positive values indicate the amount of rotation, in degrees, anticlockwise
#' from the positive x-axis.
#' @param ic.scale logical If TRUE, the height of each column is proportional
#' to its information content. Otherwise, all columns have the same height.
#' @param default.units A string indicating the default units to use if x, y,
#' width, or height are only given as numeric vectors.
#' @param name A character value to uniquely identify the motifGrob once it has
#' been pushed onto the grob tree.
#' @param gp A gpar object, typically the output from a call to the function
#' gpar. The list will be used as parameter of plotMotifLogoA.
#' @return An gTree object.
#' @author Jianhong Ou
#' @export
#' @importFrom grid unit gpar viewport gTree
#' @examples
#' 
#' pcm<-matrix(runif(40,0,100),nrow=4,ncol=10)
#' pfm<-pcm2pfm(pcm)
#' rownames(pfm)<-c("A","C","G","T")
#' motif <- new("pfm", mat=pfm, name="bin_SOLEXA")
#' motifGrob(motif)
#' 
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



#' geom_motif
#' 
#' geom_motif uses the locations of the four corners (xmin, xmax, ymin and
#' ymax) to plot motifs.
#' 
#' 
#' @param mapping Set of aesthetic mappings created by aes() or aes_(). If
#' specified and inherit.aes = TRUE (the default), it is combined with the
#' default mapping at the top level of the plot.  You must supply mapping if
#' there is no plot mapping.
#' @param data The data to be displayed in this layer.
#' @param stat The statistical transformation to use on the data for this
#' layer, as a string.
#' @param position Position adjustment, either as a string, or the result of a
#' call to a position adjustment function.
#' @param \dots Other arguments passed on to layer().
#' @param ic.scale logical If TRUE, the height of each column is proportional
#' to its information content.  Otherwise, all columns have the same height.
#' @param use.xy logical If TRUE, the required aesthethics will be x, y, width,
#' height, and motif.  Otherwise, xmin, ymin, xmax, ymax and motif.
#' @param show.legend Not used.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than
#' combining with them.
#' @return a layer that contains GeomMotif object.
#' @section Aesthetics: geom_motif() understands the following aesthetics
#' (required aesthetics are in bold): \itemize{ \item \strong{\code{xmin}}
#' \item \strong{\code{xmax}} \item \strong{\code{ymin}} \item
#' \strong{\code{ymax}} \item \strong{\code{motif}} \item \code{angle} \item
#' \code{fontfamily} \item \code{fontface} }
#' 
#' OR
#' 
#' \itemize{ \item \strong{\code{x}} \item \strong{\code{y}} \item
#' \strong{\code{width}} \item \strong{\code{height}} \item
#' \strong{\code{motif}} \item \code{angle} \item \code{fontfamily} \item
#' \code{fontface} }
#' @author Jianhong Ou
#' @export
#' @importFrom ggplot2 layer
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), 
#'                             "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
#' 
#' df <- data.frame(xmin=c(.25, .25), ymin=c(.25, .75), xmax=c(.75, .75), ymax=c(.5, 1))
#' df$motif <- list(pcm2pfm(motif), pcm2pfm(motif))
#' 
#' library(ggplot2)
#' ggplot(df, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, motif=motif)) + 
#' geom_motif() + theme_bw() + ylim(0, 1) + xlim(0, 1)
#' 
geom_motif <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       ...,
                       ic.scale=TRUE,
                       use.xy=FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = if(use.xy) GeomMotifA else GeomMotif,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      ic.scale = ic.scale,
      ...
    )
  )
}

#' GeomMotif object
#' 
#' GeomMotif object is a ggproto object.
#' 
#' 
#' @name GeomMotif
#' @format The format is: Classes 'GeoMotif', 'Geom', 'ggproto', 'gg' <ggproto
#' object: Class GeoMotif, Geom, gg> aesthetics: function default_aes: uneval
#' draw_group: function draw_key: function draw_layer: function draw_panel:
#' function extra_params: na.rm handle_na: function non_missing_aes:
#' optional_aes: parameters: function required_aes: xmin ymin xmax ymax motif
#' setup_data: function use_defaults: function super: <ggproto object: Class
#' Geom, gg>
#' @seealso geom_motif
#' @export
#' @importFrom ggplot2 ggproto aes Geom
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), 
#'                             "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
#' 
#' df <- data.frame(xmin=c(.25, .25), ymin=c(.25, .75), xmax=c(.75, .75), ymax=c(.5, 1))
#' df$motif <- list(pcm2pfm(motif), pcm2pfm(motif))
#' 
#' library(ggplot2)
#' 
#' ggplot(df, aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, motif=motif)) + 
#' geom_motif() + theme_bw() + ylim(0, 1) + xlim(0, 1)
#' 
#' 

GeomMotif <- 
  ggproto("GeomMotif", Geom, 
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

GeomMotifA <- 
  ggproto("GeomMotif", Geom, 
          required_aes = c("x", "y", "width", "height", "motif"),
          default_aes = aes(angle=0, fontfamily="Helvetica-Bold", 
                            fontface="bold"),
          draw_panel = function(data, panel_params, coord, ic.scale=TRUE){
            if(is.unit(data$x)) data$x <- convertUnit(data$x, unitTo = "native", valueOnly = TRUE)
            if(is.unit(data$y)) data$y <- convertUnit(data$x, unitTo = "native", valueOnly = TRUE)
            if(is.unit(data$width)) data$width <- convertUnit(data$width, unitTo = "native", valueOnly = TRUE)
            if(is.unit(data$height)) data$height <- convertUnit(data$height, unitTo = "native", valueOnly = TRUE)
            hw <- data$width/2
            hh <- data$height/2
            data$xmin <- data$x - hw
            data$ymin <- data$y - hh
            data$xmax <- data$x + hw
            data$ymax <- data$y + hh
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
