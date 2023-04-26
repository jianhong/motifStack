###############################################################################
######## phylog style stack
######## 
###############################################################################



#' plot sequence logo stacks with a ape4-style phylogenic tree
#' 
#' plot sequence logo stacks with a ape4-style phylogenic tree
#' 
#' 
#' @param phylog an object of class phylog
#' @param pfms a list of objects of class pfm
#' @param f.phylog a size coefficient for tree size (a parameter to draw the
#' tree in proportion to leaves label)
#' @param f.logo a size coefficient for the motif
#' @param cleaves a character size for plotting the points that represent the
#' leaves, used with par("cex")*cleaves. If zero, no points are drawn
#' @param cnodes a character size for plotting the points that represent the
#' nodes, used with par("cex")*cnodes. If zero, no points are drawn
#' @param labels.leaves a vector of strings of characters for the leaves labels
#' @param clabel.leaves a character size for the leaves labels, used with
#' par("cex")*clavel.leaves
#' @param labels.nodes a vector of strings of characters for the nodes labels
#' @param clabel.nodes a character size for the nodes labels, used with
#' par("cex")*clabel.nodes. If zero, no nodes labels are drawn
#' @param font font of logo
#' @param ic.scale logical If TRUE, the height of each column is proportional
#' to its information content. Otherwise, all columns have the same height.
#' @return none
#' @seealso \link[ade4:plot.phylog]{plot.phylog}
#' @export
#' @importFrom grid grid.newpage pushViewport viewport grid.text gpar
#' grid.segments popViewport grid.points
#' @importFrom graphics par
#' @importFrom grDevices grey
#' @examples
#' 
#'   if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'     library("MotifDb")
#'     matrix.fly <- query(MotifDb, "Dmelanogaster")
#'     motifs <- as.list(matrix.fly)
#'     motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", names(motifs), fixed=TRUE)]
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
#'     plotMotifStackWithPhylog(phylog, pfms, f.phylog=0.3, 
#'                              cleaves = 0.5, clabel.leaves = 0.7)
#'   }
#' 
plotMotifStackWithPhylog <- function(phylog, pfms=NULL,
                                     f.phylog = 0.3, f.logo = NULL, 
                                     cleaves =1, cnodes =0,
                                     labels.leaves = names(phylog$leaves), 
                                     clabel.leaves=1,
                                     labels.nodes = names(phylog$nodes), 
                                     clabel.nodes = 0, 
                                     font="sans", ic.scale=TRUE
){
  if(!inherits(phylog, "phylog")) stop("phylog must be an object of phylog")
  n<-length(pfms)
  if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) {
    pfms <- lapply(pfms, pcm2pfm)
  }
  lapply(pfms,function(.ele){
    if(!inherits(.ele, c("pfm", "psam"))) stop("pfms must be a list of pfm or psam objects.")
  })
  leaves.number <- length(phylog$leaves)
  leaves.names<- names(phylog$leaves)
  nodes.number <- length(phylog$nodes)
  nodes.names <- names(phylog$nodes)
  if (length(labels.leaves) != leaves.number) labels.leaves <- names(phylog$leaves)
  if (length(labels.nodes) != nodes.number) labels.nodes <- names(phylog$nodes)
  leaves.car <- gsub("[_]"," ",labels.leaves)
  nodes.car <- gsub("[_]"," ",labels.nodes)
  #opar <- par(mar = c(0, 0, 0, 0))
  #on.exit(par(opar))
  
  if (f.phylog < 0.05) f.phylog <- 0.05 
  if (f.phylog > 0.75) f.phylog <- 0.75
  
  maxx <- max(phylog$droot)
  # plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
  #              yaxt = "n", xlim = c(-maxx*0.15, maxx/f.phylog), ylim = c(-0, 1), xaxs = "i", 
  #              yaxs = "i", frame.plot = FALSE)
  plot <- gList()
  vp <- viewport(xscale = c(-maxx*0.15, maxx/f.phylog))

  x.leaves <- phylog$droot[leaves.names]
  x.nodes <- phylog$droot[nodes.names]
  
  y <- (leaves.number:1)/(leaves.number + 1)
  names(y) <- leaves.names
  
  xcar <- maxx*1.05
  xx <- c(x.leaves, x.nodes)
  
  assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
  vpheight <- y[length(y)] * .95
  if(is.null(f.logo)){
    f.logo <- max(unlist(lapply(leaves.names, strwidth, units="figure", cex=clabel.leaves)))
    if(!is.null(pfms)){
      pfms.width <- max(sapply(pfms, function(x) ncol(x@mat)))
      pfms.width <- strwidth(paste(rep("M", pfms.width), collapse=""), units="figure")
      f.logo <- max(f.logo, pfms.width)
    }
    if(clabel.leaves>0) f.logo <- f.phylog+f.logo+.01 else f.logo <- f.phylog+.05
  }else{
    if (f.logo > 0.85) f.logo <- 0.85
    if (f.logo < f.phylog) f.logo <- f.phylog+.05 
  }
  subplot <- gList()
  for (i in seq.int(leaves.number)) {
    if(clabel.leaves>0) 
      subplot <- 
        gList(subplot, 
              textGrob(x=xcar, y=y[i], label=leaves.car[i], just = 0, 
                gp=gpar(cex = par("cex") * clabel.leaves),
                default.units = "native"))
    subplot <- gList(subplot, 
                          segmentsGrob(x0=xcar, y0=y[i], x1=xx[i], y1=y[i], 
                                       gp = gpar(col = grey(0.7)), 
                                       default.units = "native"))
    
    if(!is.null(pfms)){
      vpwidth <- vpheight * ncol(pfms[[i]]@mat) / 2
      if(!is.null(pfms[[i]])) 
        if(!is(pfms[[i]], "psam")){
          subplot <- 
            gList(subplot,
                  gTree(children = plotMotifLogoA(pfms[[i]],
                                                  font=font,
                                                  ic.scale=ic.scale,
                                                  draw = FALSE),
                        vp = viewport(x=f.logo, y=y[i],
                                      width=vpwidth,
                                      height=vpheight,
                                      just=c(0, .5))))
        }else{
          subplot <- 
            gList(subplot,
                  gTree(children = plotAffinityLogo(pfms[[i]],
                                                    font=font,
                                                    newpage=FALSE,
                                                    draw = FALSE),
                        vp = viewport(x=f.logo,
                                      y=y[i],
                                      width=vpwidth,
                                      height=vpheight,
                                      just=c(0, .5))))
        }
    }
  }
  rm(list="tmp_motifStack_symbolsCache", envir=.globals)
  
  yleaves <- y[1:leaves.number]
  xleaves <- xx[1:leaves.number]
  if (cleaves > 0) {
    for (i in seq.int(leaves.number)) {
      subplot <- gList(subplot,
                       pointsGrob(x=xx[i], y=y[i], pch = 21,
                                  gp=gpar(fill=1,
                                          cex = par("cex") * cleaves * .5), 
                                   default.units = "native"))
    }
  }
  yn <- rep(0, nodes.number)
  names(yn) <- nodes.names
  y <- c(y, yn)
  for (i in seq_along(phylog$parts)) {
    w <- phylog$parts[[i]]
    but <- names(phylog$parts)[i]
    y[but] <- mean(y[w])
    b <- range(y[w])
    subplot <- gList(subplot,
                     segmentsGrob(x0=xx[but],
                                  y0=b[1],
                                  x1=xx[but],
                                  y1=b[2],
                                  default.units = "native"))
    x1 <- xx[w]
    y1 <- y[w]
    x2 <- rep(xx[but], length(w))
    subplot <- gList(subplot,
                     segmentsGrob(x0=x1,
                                  y0=y1,
                                  x1=x2,
                                  y1=y1, 
                                  default.units = "native"))
  }
  if (cnodes > 0) {
    for (i in nodes.names) {
      subplot <- gList(subplot,
                       pointsGrob(x=xx[i], y=y[i],
                                  pch = 21,
                                  gp=gpar(fill="white", cex = cnodes)))
    }
  }
  if (clabel.nodes > 0) {
    subplot <- gList(subplot,
                     grid.eti(xx[names(x.nodes)],
                              y[names(x.nodes)],
                              nodes.car,
                              clabel.nodes,
                              draw = FALSE))
  }
  if (cleaves > 0){
    subplot <- gList(subplot,
                     pointsGrob(x=xleaves, y=yleaves, pch = 21, 
                                gp = gpar(fill=1, cex = par("cex") * cleaves *.5),
                                default.units = "native"))
  } 
  plot <- gTree(children = subplot, vp = vp)
  suppressWarnings(grid.draw(plot))
  return(invisible(plot))
}

