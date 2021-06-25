#' plot sequence logos stack with hierarchical cluster tree
#' 
#' plot sequence logos stack with hierarchical cluster tree
#' 
#' 
#' @param pfms a list of position frequency matrices, pfms must be a list of
#' class pfm
#' @param hc an object of the type produced by stats::hclust
#' @param treewidth the width to show tree
#' @param trueDist logical flags to use hclust height or not.
#' @param ... other parameters can be passed to plotMotifLogo function
#' @return none
#' @export
#' @importFrom grid grid.newpage grid.lines viewport pushViewport convertUnit 
#' unit grid.text popViewport
#' @examples
#' 
#'   #####Input#####
#'   pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
#' 
#'   #####Clustering#####
#'   hc <- clusterMotifs(pcms)
#' 
#'   ##reorder the motifs for plotMotifLogoStack
#'   motifs<-pcms[hc$order]
#'   motifs <- lapply(motifs, pcm2pfm)
#'   ##do alignment
#'   motifs<-DNAmotifAlignment(motifs)
#'   ##plot stacks
#'   plotMotifLogoStack(motifs, ncex=1.0)
#'   plotMotifLogoStackWithTree(motifs, hc=hc)
#' 
plotMotifLogoStackWithTree<-function(pfms, hc, treewidth=1/8, trueDist=FALSE, ...){
  n<-length(pfms)
  if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) {
    pfms <- lapply(pfms, pcm2pfm)
  }
  lapply(pfms,function(.ele){
    if(!inherits(.ele, c("pfm", "psam"))) stop("pfms must be a list of pfm or psam objects")
  })
  if(!is(hc, "hclust")) stop("hc must be hclust object")
  if(treewidth>0.5) stop("treewidth can not greater than 0.5")
  #opar<-par(mar=c(0,0,0,0), mfrow=par("mfrow"))
  #layout(matrix(c(rep(1,n),rep(2:(n+1),ceiling(1/treewidth)-1)),nrow=n,ncol=ceiling(1/treewidth)))
  #plot tree
  
  dots <- list(...)
  if(is.null(dots$draw)) dots$draw <- TRUE
  if(is.null(dots$newpage)) dots$newpage <- TRUE
  
  if(trueDist) h <- hc$height / max(hc$height) / 1.05
  else h<- seq(0.01, 0.95, length.out=length(hc$height))[order(hc$height)]
  m <- hc$merge
  o <- hc$order
  
  m[m > 0] <- n + m[m > 0] 
  m[m < 0] <- abs(m[m < 0])
  
  dist <- matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
  dist[1:n, 2] <- 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
  
  for(i in 1:nrow(m)){
    dist[n + i, 2] <- (dist[m[i, 1], 2] + dist[m[i, 2], 2]) / 2
    dist[n + i, 1] <- h[i]
  }
  
  trx <- function(x) 0.9*treewidth*(1 - x)+ 0.01
  for(i in 1:nrow(dist)){
    dist[i, 1] <- trx(dist[i, 1])
    #dist[i, 2] <- dist[i, 2] + 1/(4*n)
  }
  
  draw_connection <- function(x1, x2, y1, y2, x){
    gList(
      linesGrob(x = c(x1, x), y = c(y1, y1)),
      linesGrob(x = c(x2, x), y = c(y2, y2)),
      linesGrob(x = c(x, x), y = c(y1, y2))
    )
  }
  plot <- gList()
  for(i in 1:nrow(m)){
    plot <- gList(plot, 
                  draw_connection(dist[m[i, 1], 1],
                                  dist[m[i, 2], 1],
                                  dist[m[i, 1], 2],
                                  dist[m[i, 2], 2], trx(h[i])))
  }
  
  #plot logo
  subplots <- gList()
  if(all(sapply(pfms, function(.ele) is(.ele, "psam")))){
    .local <- function(..., newpage, draw){
      plotAffinityLogo(..., newpage = FALSE, draw = FALSE)
    }
    ht <- 1/n
    y0 <- .5 * ht
    for(i in seq.int(n)){
      ht.title <- convertUnit(unit(1.5, "lines"), unitTo = "npc", valueOnly = TRUE)
      ht.body <- 1 - ht.title
      plotsub <- gTree(children = .local(psam=pfms[[i]], ...), 
                       vp = viewport(y=ht.body*.5, height=ht.body))
      subplots <- gList(subplots, 
                    gTree(children = gList(textGrob(label=pfms[[i]]@name,
                                                     y=1-.5*ht.title),
                                           plotsub),
                          vp = viewport(y=y0, height=ht)))
      y0 <- y0 + ht
    }
  }else{
    lapply(pfms,function(.ele){
      if(!is(.ele, "pfm")) stop("pfms must be a list of pfm objects.")
    })
    assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
    .local <- function(..., newpage, draw){
      plotMotifLogo(..., newpage = FALSE, draw = FALSE)
    }
    ht <- 1/n
    y0 <- .5 * ht
    for(i in seq.int(n)){
      subplots <- 
        gList(subplots, gTree(children=
                            .local(pfm=pfms[[i]], motifName=pfms[[i]]@name, 
                                   p=pfms[[i]]@background,
                                   colset = pfms[[i]]@color,
                                   xlab=NA, newpage=FALSE,
                                   margins=c(1.5, 4.1, 1.1, .1), ...),
                          vp= viewport(y=y0, height=ht)))
      y0 <- y0 + ht
    }
    rm(list="tmp_motifStack_symbolsCache", envir=.globals)
  }
  plot <- 
    gList(plot, gTree(children = subplots,
                      vp = viewport(x = .5 + treewidth/2, y = .5,
                                    width = 1-treewidth, height = 1)))
  if(dots$draw){
    if(dots$newpage) grid.newpage()
    suppressWarnings(grid.draw(plot))
  }else{
    plot
  }
}
