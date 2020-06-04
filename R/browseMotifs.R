#' browse motifs
#' 
#' browse motifs in a web browser
#' 
#' 
#' @param pfms a list of \link{pfm}
#' @param phylog layout type. see \link[Rgraphviz]{GraphvizLayouts}
#' @param layout layout type. Could be tree, cluster or radialPhylog.
#' @param nodeRadius node radius, default 2.5px.
#' @param baseWidth,baseHeight width and height of each alphabet of the motif
#' logo.
#' @param xaxis,yaxis plot x-axis or y-axis or not in the motifs.
#' @param width width of the figure
#' @param height height of the figure
#' @param ... parameters not used
#' @return An object of class htmlwidget that will intelligently print itself
#' into HTML in a variety of contexts including the R console, within R
#' Markdown documents, and within Shiny output bindings.
#' @keywords plot
#' @export
#' @importFrom ade4 hclust2phylog
#' @importFrom htmlwidgets createWidget
#' @examples
#' 
#' library("MotifDb")
#' matrix.fly <- query(MotifDb, "Dmelanogaster")
#' motifs <- as.list(matrix.fly)
#' motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", 
#'                        names(motifs), fixed=TRUE)]
#' names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
#'                      gsub("_FBgn[0-9]+$", "", 
#'                           gsub("[^a-zA-Z0-9]","_", 
#'                                gsub("(_[0-9]+)+$", "", names(motifs)))))
#' motifs <- motifs[unique(names(motifs))]
#' pfms <- sample(motifs, 10)
#' pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#' browseMotifs(pfms)
#' 
browseMotifs <- function(pfms, phylog,
                      layout=c("tree", "cluster", "radialPhylog"),
                      nodeRadius=2.5, baseWidth=12, baseHeight=30,
                      xaxis=TRUE, yaxis=TRUE,
                      width=NULL, height=NULL,
                      ...){
  if(!is.list(pfms)){
    pfms <- list(pfms)
  }
  layout <- match.arg(layout)
  if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) {
    pfms <- lapply(pfms, pcm2pfm)
  }
  if (any(unlist(lapply(pfms, function(.ele) !is(.ele, "pfm"))))) 
    stop("pfms must be a list of pfm objects")
  maxW <- 0
  if(length(pfms)<=2){
    elements <- lapply(pfms, function(pfm){
      w <- (ncol(pfm$mat)+1)*baseWidth
      if(w > maxW) maxW <<- w
      list(name="Root",
           motif=list(t(pfm$mat)),
           letters=rownames(pfm$mat),
           color=as.list(pfm$color),
           background=as.list(pfm$background),
           width=w,
           height=baseHeight,
           xaxis=xaxis,
           yaxis=yaxis)
    })
    names(elements) <- NULL
  }else{
    if(missing(phylog)){
      hc <- clusterMotifs(pfms)
      pfms <- pfms[hc$order]
      #pfms <- DNAmotifAlignment(pfms)
      phylog <- hclust2phylog(hc)
    }
    stopifnot(inherits(phylog, "phylog"))
    names(pfms) <- names(phylog$leaves)
    ## list to recursive list
    flatListToRecursiveList <- function(flatList, nodeName="Root"){
      if(nodeName %in% names(flatList)){## nodes
        out <- lapply(flatList[[nodeName]], 
                      flatListToRecursiveList, 
                      flatList=flatList)
        names(out) <- NULL
        ele <- list(name=nodeName, children=out)
      }else{## leaves
        pfm <- pfms[[nodeName]]
        w <- (ncol(pfm$mat)+1)*baseWidth
        if(w > maxW) maxW <<- w
        ele <- list(name=nodeName,
                    motif=list(t(pfm$mat)),
                    letters=rownames(pfm$mat),
                    color=as.list(pfm$color),
                    background=as.list(pfm$background),
                    width=w,
                    height=baseHeight,
                    xaxis=xaxis,
                    yaxis=yaxis)
      }
      return(ele)
    }
    elements <- flatListToRecursiveList(phylog$parts)
  }
  
  x <- list(
    elements = elements,
    nodeRadius = nodeRadius,
    maxH = baseHeight,
    maxW = maxW,
    layout = layout
  )
  
  htmlwidgets::createWidget(
    name = 'browseMotifs',
    x = x,
    width = width,
    height = height,
    package = getPackageName()
  )
}

#' Shiny bindings for browseMotifs
#'
#' Output and render functions for using browseMotifs within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a browseMotifs
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name browseMotifs-shiny
#'
#' @export
#' @importFrom htmlwidgets shinyWidgetOutput
browseMotifsOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'browseMotifs', width, height, 
                                 package = 'motifStack')
}

#' @rdname browseMotifs-shiny
#' @export
#' @importFrom htmlwidgets shinyRenderWidget
renderbrowseMotifs <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, browseMotifsOutput, env, quoted = TRUE)
}
