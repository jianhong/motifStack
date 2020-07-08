#' plot affinity logo
#' 
#' plot affinity logo
#' 
#' 
#' @param psam a position-specific affinity matrix
#' @param motifName motif name
#' @param font font of logo
#' @param fontface fontface of logo
#' @param colset color setting for each logo letter
#' @param alpha Alpha channel for transparency of low affinity letters.
#' @param newpage plot in a new canvas or not.
#' @param draw Vector (logical(1)). TRUE to plot. FALSE, return a gList
#' @return none
#' @references Barrett C. Foat, Alexandre V. Morozov, Harmen J. Bussemaker;
#' Statistical mechanical modeling of genome-wide transcription factor
#' occupancy data by MatrixREDUCE, Bioinformatics, Volume 22, Issue 14, 15 July
#' 2006, Pages e141-e149, https://doi.org/10.1093/bioinformatics/btl223
#' @export
#' @importFrom grid gList linesGrob unit gpar grid.newpage grid.draw 
#' @examples
#' 
#' psam <- importMatrix(file.path(find.package("motifStack"), "extdata", "PSAM.mxr"), 
#'                      format="psam")[[1]]
#' plotAffinityLogo(psam)
#' 
plotAffinityLogo <- function(psam, motifName, font="Helvetica-Bold", fontface="bold",
                             colset=c("#00811B","#2000C7","#FFB32C","#D00001"),
                             alpha=0.5, newpage=TRUE, draw=TRUE){
  markers <- NULL
  pssm <- FALSE
  if(is(psam, "data.frame")){
    psam <- as.matrix(psam)
  }else{
    if(inherits(psam, c("psam", "pssm"))){
      markers <- psam@markers 
      if(missing(motifName)) motifName = psam@name
      colset=psam@color[rownames(psam@mat)]
      pssm <- is(psam, "pssm")
      psam <- psam@mat
    }
  }
  if (!is(psam, "matrix")){
    stop("psam must be a matrix.")
  }
  if(!is(colset, "character"))
    stop("colset must be a character")
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
  pictureGrob <- get("pictureGrob", envir = .globals)
  
  #calculate postion of each symbol and plot
  plot <- gList()
  
  ddG <- if(pssm) psam else log(psam+1e-2)
  ddG.mu<-colMeans(ddG)
  ddG.height <- t(t(ddG) - ddG.mu)
  ddG.h.pos <- ddG.h.neg <- ddG.height
  ddG.h.pos[ddG.h.pos<0] <- 0
  ddG.h.neg[ddG.h.neg>0] <- 0
  ie <- max(abs(min(colSums(ddG.h.neg))), max(colSums(ddG.h.pos)))
  ddG.height <- ddG.height/(2*ie)
  dw<-1/npos
  x.pos<-0
  y.poss <- numeric(length=npos)
  y.low.poss <- numeric(length=npos)
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
        plot <- 
          gList(plot, 
                pictureGrob(picture=symbols[[paste0(id[i], "_", alpha)]],
                            x = x.pos,
                            y = y.pos,
                            width = dw,
                            height = -h,
                            just=c(0,0),
                            distort=TRUE))
      }
    }
    y.low.poss[j] <- y.pos
    ## greater than 0
    id<-order(heights, decreasing = FALSE)
    y.pos<-0.5
    for(i in seq.int(ncha)){
      h<-heights[id[i]]
      if(h>0) {
        plot <- 
          gList(plot, 
                pictureGrob(picture = symbols[[id[i]]],
                            x = x.pos,
                            y = y.pos,
                            width = dw,
                            height =h,
                            just=c(0,0),
                            distort=TRUE))
        y.pos<-y.pos+h
      }
    }
    y.poss[j] <- y.pos
    x.pos<-x.pos+dw
  }
  if(length(markers)>0){
    plot <- gList(plot, plotMarkers(markers, dw, y.poss, y.low.poss))
  }
  
  if(draw){
    if(newpage) grid.newpage()
    suppressWarnings(grid.draw(plot))
  }else{
    plot
  }
}

