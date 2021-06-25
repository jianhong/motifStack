#' plot sequence logos stack
#' 
#' plot sequence logos stack
#' 
#' 
#' @param pfms a list of position frequency matrices, pfms must be a list of
#' class pfm
#' @param ... other parameters can be passed to plotMotifLogo function
#' @return none
#' @export
#' @importFrom grid grid.newpage pushViewport viewport convertUnit unit
#' grid.text popViewport
#' @examples
#' 
#'   pcm1<-matrix(c(0,50,0,50,
#'            100,0,0,0,
#'            0,100,0,0,
#'            0,0,100,0,
#'            0,0,0,100,
#'            50,50,0,0,
#'            0,0,50,50), nrow=4)
#'   pcm2<-matrix(c(50,50,0,0,
#'            0,100,0,0,
#'            0,50,50,0,
#'            0,0,0,100,
#'            50,50,0,0,
#'            0,0,50,50), nrow=4)
#'   rownames(pcm1)<-c("A","C","G","T")
#'   rownames(pcm2)<-c("A","C","G","T")
#'   pfms<-list(p1=new("pfm",mat=pcm2pfm(pcm1),name="m1"),
#'          p2=new("pfm",mat=pcm2pfm(pcm2),name="m2"))
#'   pfms<-DNAmotifAlignment(pfms)
#'   plotMotifLogoStack(pfms)
#' 
plotMotifLogoStack<-function(pfms, ...){
  n<-length(pfms)
  dots <- list(...)
  if(is.null(dots$draw)) dots$draw <- TRUE
  if(is.null(dots$newpage)) dots$newpage <- TRUE
  if(all(sapply(pfms, function(.ele) is(.ele, "psam")))){
    .local <- function(..., newpage, draw){
      plotAffinityLogo(..., newpage = FALSE, draw = FALSE)
    }
    plot <- gList()
    ht <- 1/n
    y0 <- .5 * ht
    for(i in seq.int(n)){
      ht.title <- convertUnit(unit(1.5, "lines"), unitTo = "npc", valueOnly = TRUE)
      ht.body <- 1 - ht.title
      plotsub <- gTree(children = .local(psam=pfms[[i]], ...), 
                       vp = viewport(y=ht.body*.5, height=ht.body))
      plot <- gList(plot, 
                    gTree(children = gList(textGrob(label=pfms[[i]]@name,
                                                     y=1-.5*ht.title),
                                           plotsub),
                          vp = viewport(y=y0, height=ht)))
      y0 <- y0 + ht
    }
    if(dots$draw){
      if(dots$newpage) grid.newpage()
      suppressWarnings(grid.draw(plot))
    }else{
      return(plot)
    }
    return()
  }
  if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) {
    pfms <- lapply(pfms, pcm2pfm)
  }
  lapply(pfms,function(.ele){
    if(!is(.ele, "pfm")) stop("pfms must be a list of pfm objects.")
  })
  assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
  .local <- function(..., newpage, draw){
    plotMotifLogo(..., newpage = FALSE, draw = FALSE)
  }
  plot <- gList()
  ht <- 1/n
  y0 <- .5 * ht
  for(i in seq.int(n)){
    plot <- 
      gList(plot, gTree(children=
                          .local(pfm=pfms[[i]], motifName=pfms[[i]]@name, 
                                 p=pfms[[i]]@background,
                                 colset = pfms[[i]]@color,
                                 xlab=NA, newpage=FALSE,
                                 margins=c(1.5, 4.1, 1.1, .1), ...),
                        vp= viewport(y=y0, height=ht)))
    y0 <- y0 + ht
  }
  rm(list="tmp_motifStack_symbolsCache", envir=.globals)
  if(dots$draw){
    if(dots$newpage) grid.newpage()
    suppressWarnings(grid.draw(plot))
  }else{
    plot
  }
}
