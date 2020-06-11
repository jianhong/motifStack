###############################################################################
######## radial style stack
######## 
###############################################################################


#' plot sequence logo stacks with a radial phylogenic tree
#' 
#' plot sequence logo stacks with a radial phylogenic tree
#' 
#' 
#' @param phylog an object of class phylog
#' @param pfms a list of objects of class pfm
#' @param circle a size coefficient for the outer circle of the labels. Please
#' note this is the position of inner.label.cirle.
#' @param circle.motif a size coefficient for the motif circle
#' @param cleaves a character size for plotting the points that represent the
#' leaves, used with par("cex")*cleaves. If zero, no points are drawn
#' @param cnodes a character size for plotting the points that represent the
#' nodes, used with par("cex")*cnodes. If zero, no points are drawn
#' @param labels.leaves a vector of strings of characters for the leaves labels
#' @param clabel.leaves a character size for the leaves labels, used with
#' par("cex")*clabel.leaves
#' @param labels.nodes a vector of strings of characters for the nodes labels
#' @param clabel.nodes a character size for the nodes labels, used with
#' par("cex")*clabel.nodes. If zero, no nodes labels are drawn
#' @param draw.box if TRUE draws a box around the current plot with the
#' function box()
#' @param col.leaves a vector of colors for leaves labels
#' @param col.leaves.bg a vector of colors for background of leaves labels
#' @param col.leaves.bg.alpha alpha value [0, 1] for the colors of backgroud of
#' leaves labels
#' @param col.bg a vector of colors for tree background
#' @param col.bg.alpha a alpha value [0, 1] of colors for tree background
#' @param col.inner.label.circle a vector of colors for inner cirlce of pfms
#' @param inner.label.circle.width width for inner circle of pfms
#' @param col.outer.label.circle a vector of colors for outer circle of pfms
#' @param outer.label.circle.width width for outer circle of pfms
#' @param clockwise a logical value indicating if slices are drawn clockwise or
#' counter clockwise
#' @param init.angle number specifying the starting angle (in degrees) for the
#' slices. Defaults to 0 (i.e., `3 o'clock`) unless clockwise is true where
#' init.angle defaults to 90 (degrees), (i.e., `12 o'clock`)
#' @param angle number specifying the angle (in degrees) for phylogenic tree.
#' Defaults 360
#' @param pfmNameSpliter spliter when name of pfms contain multiple node of
#' labels.leaves
#' @param rcpostfix the postfix for reverse complements
#' @param motifScale the scale of logo size
#' @param ic.scale logical. If TRUE, the height of each column is proportional
#' to its information content. Otherwise, all columns have the same height.
#' @param plotIndex logical. If TRUE, will plot index number in the motifLogo
#' which can help user to describe the motifLogo
#' @param IndexCol The color of the index number when plotIndex is TRUE.
#' @param IndexCex The cex of the index number when plotIndex is TRUE.
#' @param groupDistance show groupDistance on the draw
#' @param groupDistanceLineCol groupDistance line color, default: red
#' @param plotAxis logical. If TRUE, will plot distance axis.
#' @param font font of logo
#' @return none
#' @seealso \link[ade4:plot.phylog]{plot.phylog}
#' @export
#' @importFrom graphics par plot.default strwidth polygon text segments
#' points lines symbols box
#' @importFrom grDevices grey dev.size
#' @importFrom stats median
#' @importFrom grid pushViewport viewport popViewport grid.text gpar grid.xaxis
#' @examples
#' 
#'   if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'     library("MotifDb")
#'     matrix.fly <- query(MotifDb, "Dmelanogaster")
#'     motifs <- as.list(matrix.fly)
#'     motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", 
#'                             names(motifs), fixed=TRUE)]
#'     names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
#'                 gsub("_FBgn[0-9]+$", "", 
#'                   gsub("[^a-zA-Z0-9]","_", 
#'                      gsub("(_[0-9]+)+$", "", names(motifs)))))
#'     motifs <- motifs[unique(names(motifs))]
#'     pfms <- sample(motifs, 50)
#'     hc <- clusterMotifs(pfms)
#'     library(ade4)
#'     phylog <- ade4::hclust2phylog(hc)
#'     leaves <- names(phylog$leaves)
#'     pfms <- pfms[leaves]
#'     pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#'     pfms <- DNAmotifAlignment(pfms, minimalConsensus=3)
#'     library(RColorBrewer)
#'     color <- brewer.pal(12, "Set3")
#'     plotMotifStackWithRadialPhylog(phylog, pfms, circle=0.9, 
#'                  cleaves = 0.5, clabel.leaves = 0.7, 
#'                  col.bg=rep(color, each=5), col.leaves=rep(color, each=5))
#'   }
#' 
plotMotifStackWithRadialPhylog <- 
  function (phylog, pfms=NULL,
            circle=.75, circle.motif=NA, cleaves=1, cnodes=0,
            labels.leaves=names(phylog$leaves), clabel.leaves=1,
            labels.nodes=names(phylog$nodes), clabel.nodes=0,
            draw.box=FALSE,
            col.leaves=rep("black", length(labels.leaves)),
            col.leaves.bg=NULL, col.leaves.bg.alpha=1,
            col.bg=NULL, col.bg.alpha=1,
            col.inner.label.circle=NULL, inner.label.circle.width="default",
            col.outer.label.circle=NULL, outer.label.circle.width="default",
            clockwise =FALSE, init.angle=if(clockwise) 90 else 0,
            angle=360, pfmNameSpliter=";", rcpostfix="(RC)", 
            motifScale=c("linear","logarithmic"), ic.scale=TRUE,
            plotIndex=FALSE, IndexCol="black", IndexCex=.8,
            groupDistance=NA, groupDistanceLineCol="red", 
            plotAxis=FALSE, font="Helvetica-Bold")
{
  if (!inherits(phylog, "phylog"))
    stop("Non convenient data")
  leaves.number <- length(phylog$leaves)
  checkLength <- function(tobechecked){
    !((length(tobechecked)>=leaves.number)||is.null(tobechecked))
  }
  checkNA <- function(tobechecked){
    if(is.null(tobechecked)) return(FALSE)
    return(any(is.na(tobechecked)))
  }
  for(tobechecked in c("col.leaves", "col.leaves.bg", "col.bg", "col.inner.label.circle", "col.outer.label.circle")){
    if(checkLength(eval(as.symbol(tobechecked)))) stop(paste("the length of", tobechecked, "should be same as the length of leaves"))
    if(checkNA(eval(as.symbol(tobechecked)))) stop(paste("contain NA in", tobechecked))
  }
  motifScale <- match.arg(motifScale)
  leaves.names <- names(phylog$leaves)
  nodes.number <- length(phylog$nodes)
  nodes.names <- names(phylog$nodes)
  if (length(labels.leaves) != leaves.number)
    labels.leaves <- names(phylog$leaves)
  if (length(labels.nodes) != nodes.number)
    labels.nodes <- names(phylog$nodes)
  if (circle < 0)
    stop("'circle': non convenient value")
  leaves.car <- gsub("[_]", " ", labels.leaves)
  nodes.car <- gsub("[_]", " ", labels.nodes)
  opar <- par(mar = par("mar"), srt = par("srt"))
  on.exit(par(opar))
  par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow=c(1,1))
  dis <- phylog$droot
  max_Dis <- max(dis)
  axis_pos <- circle
  if(!is.na(groupDistance)){
    groupDistance <- (max_Dis - groupDistance) * circle / max_Dis
  }
  dis <- dis/max_Dis
  rayon <- circle
  dis <- dis * rayon
  dist.leaves <- dis[leaves.names]
  dist.nodes <- dis[nodes.names]
  asp <- c(1, 1)
  if(is.null(pfms)){
    plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", xlim = c(-2, 2), ylim = c(-2, 2),
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
  }else{
    pin <- dev.size("in") #pin <- par("pin")
    if (pin[1L] > pin[2L]) asp <- c(pin[2L]/pin[1L], 1)
    else asp <- c(1, pin[1L]/pin[2L])
    plot.default(0, 0, type = "n", asp=1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
  }
  d.rayon <- rayon/(nodes.number - 1)
  twopi <- if (clockwise) -2 * pi else 2 * pi
  alpha <- twopi * angle * (1:leaves.number)/leaves.number/360 + init.angle * pi/180
  names(alpha) <- leaves.names
  x <- dist.leaves * cos(alpha)
  y <- dist.leaves * sin(alpha)
  xcar <- (rayon + d.rayon) * cos(alpha)
  ycar <- (rayon + d.rayon) * sin(alpha)
  
  rayonWidth <- max(unlist(lapply(leaves.names, strwidth, units="user", cex=clabel.leaves)))
  circle.motif <- ifelse(is.na(circle.motif), rayon + d.rayon + rayonWidth, circle.motif)
  ##for logos position
  maxvpwidth <- 0
  if(!is.null(pfms)){
    beta <- alpha * 180 / pi
    vpheight <- 2 * pi * angle * circle.motif / 360 / leaves.number / 5
    vpheight <- vpheight * asp[2L]
    xm <- circle.motif * cos(alpha) * asp[1L] / 5 + 0.5
    ym <- circle.motif * sin(alpha) * asp[2L] / 5 + 0.5
    pfmNamesLen <- sapply(pfms, function(.ele) 
      length(strsplit(.ele@name, pfmNameSpliter)[[1]]))
    if(motifScale=="linear")
      vph <- 2.5*vpheight*pfmNamesLen
    else vph <- 2.5*vpheight*(1+log2(pfmNamesLen+0.0001))
    maxvpwidth <- max(mapply(function(.ele, f) ncol(.ele@mat)*f, pfms, vph))
  }
  if(inner.label.circle.width=="default") inner.label.circle.width <- rayonWidth/10
  if(outer.label.circle.width=="default") outer.label.circle.width <- rayonWidth/10
  
  ratio <- if(!is.null(pfms)) 2.5/(circle.motif+maxvpwidth+outer.label.circle.width) else 1
  if(ratio < 1){
    x <- x * ratio
    y <- y * ratio
    xcar <- xcar * ratio
    ycar <- ycar * ratio
    dis <- dis * ratio
    groupDistance <- groupDistance * ratio
    xm <- (xm-.5) * ratio + .5
    ym <- (ym-.5) * ratio + .5
    vpheight <- vpheight * ratio
    axis_pos <- axis_pos * ratio
  }
  ##for plot background
  if(!is.null(col.bg)) col.bg <- highlightCol(col.bg, col.bg.alpha)
  if(!is.null(col.leaves.bg)) col.leaves.bg <- highlightCol(col.leaves.bg, col.leaves.bg.alpha)
  gamma <- twopi * angle * ((1:(leaves.number+1))-0.5)/leaves.number/360 + init.angle * pi/180
  n <- max(2, floor(200*360/leaves.number))
  plotBgArc <- function(r,bgcol,inr){
    if(ratio < 1){
      r <- r*ratio
      inr <- inr*ratio
    }
    t2xy <- function(rx,t) list(x=rx*cos(t), y=rx*sin(t))
    oldcol <- bgcol[1]
    start <- 1
    icnt <- 1
    for(i in 1:leaves.number){
      oldcol <- bgcol[i]
      if(i==leaves.number || bgcol[i+1]!=oldcol){
        P <- t2xy(r, seq.int(gamma[start], gamma[start+icnt], length.out=n*icnt))
        polygon(c(P$x, 0), c(P$y, 0), border=bgcol[i], col=bgcol[i])
        start <- i+1
        icnt <- 1
      } else {
        icnt <- icnt + 1
      }
    }
    if(inr!=0){
      P <- t2xy(inr, seq.int(0, twopi, length.out=n*leaves.number))
      polygon(c(P$x, 0), c(P$y, 0), border="white", col="white")
    }
  }
  if (clabel.leaves > 0) {
    if(!is.null(col.outer.label.circle)) ##plot outer.label.circle
      plotBgArc(circle.motif+maxvpwidth+outer.label.circle.width, col.outer.label.circle, circle.motif+maxvpwidth)
    if(!is.null(col.inner.label.circle)) #plot inner.label.circle
      plotBgArc(circle.motif, col.inner.label.circle, circle.motif - inner.label.circle.width)
    if(!is.null(col.leaves.bg)) ##plot leaves bg
      plotBgArc(circle.motif - inner.label.circle.width, col.leaves.bg, rayon+d.rayon)
    if(!is.null(col.bg)) ##plot center bg
      plotBgArc(ifelse(mean(dist.leaves)/max(dis) > .9, mean(dist.leaves), max(dis)-d.rayon), col.bg, 0)##plotBgArc(mean(dist.leaves), col.bg, 0)
    for(i in 1:leaves.number) {
      par(srt = alpha[i] * 180/pi)
      text(xcar[i], ycar[i], leaves.car[i], adj = 0, col=col.leaves[i], cex = par("cex") *
             clabel.leaves)
      segments(xcar[i], ycar[i], x[i], y[i], col = grey(0.7))
    }
    
    assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
    if(!is.null(pfms)){
      ##extract names
      for(metaChar in c("\\","$","*","+",".","?","[","]","^","{","}","|","(",")"))
      {
        rcpostfix <- gsub(metaChar,paste("\\",metaChar,sep=""),rcpostfix,fixed=TRUE)
      }
      pfmNames <- lapply(pfms, function(.ele) .ele@name)
      for(i in 1:length(pfmNames)){
        pfmname <- unlist(strsplit(pfmNames[[i]], pfmNameSpliter))
        pfmname <- gsub(paste(rcpostfix,"$",sep=""),"",pfmname)
        pfmIdx <- which(makeLeaveNames(labels.leaves) %in% makeLeaveNames(pfmname))
        if(length(pfmIdx)==0) 
          pfmIdx <- which(makeLeaveNames(names(phylog$leaves)) 
                          %in% makeLeaveNames(pfmname))
        if(length(pfmIdx)>0){
          vph <- ifelse(motifScale=="linear",
                        vpheight*length(pfmname),
                        vpheight*(1+log2(length(pfmname))))
          vpw <- vph * ncol(pfms[[i]]@mat) / 2
          vpd <- sqrt(vph*vph+vpw*vpw) / 2
          angle <- median(beta[pfmIdx])
          if(length(pfmIdx)%%2==1){
            this.pfmIdx <- which(beta[pfmIdx] == angle)[1]
            vpx <- xm[pfmIdx[this.pfmIdx]] + vpd * cos(alpha[pfmIdx[this.pfmIdx]]) * asp[1L]
            vpy <- ym[pfmIdx[this.pfmIdx]] + vpd * sin(alpha[pfmIdx[this.pfmIdx]]) * asp[2L]
            vpx1 <- xm[pfmIdx[this.pfmIdx]] - inner.label.circle.width * cos(alpha[pfmIdx[this.pfmIdx]]) * asp[1L] *1.1/5
            vpy1 <- ym[pfmIdx[this.pfmIdx]] - inner.label.circle.width * sin(alpha[pfmIdx[this.pfmIdx]]) * asp[2L] *1.1/5
          }else{
            this.pfmIdx <- order(abs(beta[pfmIdx] - angle))[1:2]
            vpx <- median(xm[pfmIdx[this.pfmIdx]]) + vpd * cos(median(alpha[pfmIdx[this.pfmIdx]])) * asp[1L]
            vpy <- median(ym[pfmIdx[this.pfmIdx]]) + vpd * sin(median(alpha[pfmIdx[this.pfmIdx]])) * asp[2L]
            vpx1 <- median(xm[pfmIdx[this.pfmIdx]]) - inner.label.circle.width * cos(median(alpha[pfmIdx[this.pfmIdx]])) * asp[1L] *1.1/5
            vpy1 <- median(ym[pfmIdx[this.pfmIdx]]) - inner.label.circle.width * sin(median(alpha[pfmIdx[this.pfmIdx]])) * asp[2L] *1.1/5
          }
          pushViewport(viewport(x=vpx, y=vpy, width=vpw, height=vph, angle=angle))
          if(!is(pfms[[i]], "psam")){
            plotMotifLogoA(pfms[[i]], font=font, ic.scale=ic.scale)
          }else{
            plotAffinityLogo(pfms[[i]], font=font, newpage=FALSE)
          }
          popViewport()
          if(plotIndex) {
            grid.text(label=i, x=vpx1, 
                      y=vpy1, 
                      gp=gpar(col=IndexCol, cex=IndexCex), rot=angle, just="right")
          }
        }else{
          warning(paste("No leave named as ", paste(pfmname, collapse=", ")), sep="")
        }
      }
    }
    rm(list="tmp_motifStack_symbolsCache", envir=.globals)
  }
  if (cleaves > 0) {
    for (i in 1:leaves.number) points(x[i], y[i], pch = 21, col=col.leaves[i],
                                      bg = col.leaves[i], cex = par("cex") * cleaves)
  }
  ang <- rep(0, length(dist.nodes))
  names(ang) <- names(dist.nodes)
  ang <- c(alpha, ang)
  for (i in 1:length(phylog$parts)) {
    w <- phylog$parts[[i]]
    but <- names(phylog$parts)[i]
    ang[but] <- mean(ang[w])
    b <- range(ang[w])
    a.seq <- c(seq(b[1], b[2], by = pi/180), b[2])
    lines(dis[but] * cos(a.seq), dis[but] * sin(a.seq), col="#222222", lwd=par("cex") * clabel.leaves)
    x1 <- dis[w] * cos(ang[w])
    y1 <- dis[w] * sin(ang[w])
    x2 <- dis[but] * cos(ang[w])
    y2 <- dis[but] * sin(ang[w])
    segments(x1, y1, x2, y2, col="#222222", lwd=par("cex") * clabel.leaves)
  } 
  if(!is.na(groupDistance)){
    if(length(groupDistanceLineCol)!=length(groupDistance)){
      groupDistanceLineCol <- rep(groupDistanceLineCol, ceiling(length(groupDistance)/length(groupDistanceLineCol)))[1:length(groupDistance)]
    }
    for(i in 1:length(groupDistance)){
      symbols(x=0, y=0, circles=groupDistance[i], fg=groupDistanceLineCol[i], lty=2, inches=FALSE, add=TRUE)
    }
  }
  if (cnodes > 0) {
    for (i in 1:length(phylog$parts)) {
      w <- phylog$parts[[i]]
      but <- names(phylog$parts)[i]
      ang[but] <- mean(ang[w])
      points(dis[but] * cos(ang[but]), dis[but] * sin(ang[but]),
             pch = 21, bg = "white", cex = par("cex") * cnodes)
    }
  }
  
  ## draw distance indix  
  if(plotAxis){
    wd <- if(!is.null(pfms)) 5 else 4
    if(clockwise){
      vp <- viewport(x=0.5, y=0.5, width=axis_pos/wd, height=.1, 
                     xscale=c(0, max_Dis), angle=init.angle, just=c(0, 0))
      pushViewport(vp)
      grid.xaxis(gp=gpar(cex = par("cex") * clabel.leaves, col="lightgray"),
                 main=TRUE)
      popViewport()
    }else{
      vp <- viewport(x=0.5, y=0.5, width=axis_pos/wd, height=.1, 
                     xscale=c(0, max_Dis), angle=init.angle, just=c(0, 1))
      pushViewport(vp)
      grid.xaxis(gp=gpar(cex = par("cex") * clabel.leaves, col="lightgray"),
                 main=FALSE)
      popViewport()
    }
  }
  
  points(0, 0, pch = 21, cex = par("cex") * 2 * clabel.leaves, bg = "red")
  if (clabel.nodes > 0) {
    delta <- strwidth(as.character(length(dist.nodes)), cex = par("cex") *
                        clabel.nodes)
    for (j in 1:length(dist.nodes)) {
      i <- names(dist.nodes)[j]
      par(srt = (ang[i] * 360/2/pi + 90))
      x1 <- dis[i] * cos(ang[i])
      y1 <- dis[i] * sin(ang[i])
      symbols(x1, y1, delta, bg = "white", add = TRUE,
              inches = FALSE)
      text(x1, y1, nodes.car[j], adj = 0.5, cex = par("cex") *
             clabel.nodes)
    }
  }
  
  if (draw.box)
    box()
  return(invisible())
}
