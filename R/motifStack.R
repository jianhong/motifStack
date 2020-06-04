###############################################################################
####### motifstack


#' plot a DNA sequence logo stack
#' 
#' Plot a DNA sequence logo stack
#' 
#' 
#' @param pfms a list of objects of class \link{pfm}
#' @param layout layout of the logo stack, stack, treeview or radialPhylog
#' @param \dots any parameters could to pass to \link{plotMotifLogoStack},
#' \link{plotMotifLogoStackWithTree}, \link{plotMotifStackWithPhylog} or
#' \link{plotMotifStackWithRadialPhylog}
#' @return return a list contains pfms and phylog
#' @export
#' @importFrom ade4 hclust2phylog
#' @examples
#' 
#'   if(interactive()){
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
#'     pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#'     motifStack(pfms, "radialPhylog")
#'   }
#' 
motifStack <-function(pfms, 
                      layout=c("stack", "treeview", "phylog", "radialPhylog"), 
                      ...){
    if(!is.list(pfms)){
      if(is(pfms, "pcm")) pfms <- pcm2pfm(pfms)
        plot(pfms)
        return(invisible())
    }
    layout <- match.arg(layout)
    if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) {
      pfms <- lapply(pfms, pcm2pfm)
    }
    if (any(unlist(lapply(pfms, function(.ele) 
      !inherits(.ele, c("pfm", "psam")))))) 
        stop("pfms must be a list of pfm, pcm or psam objects")
    if(all(sapply(pfms, function(.ele) is(.ele, "psam")))){
      psam <- TRUE
    }else{
      psam <- FALSE
    }
    if (length(pfms)<2)
        stop("length of pfms less than 2")
    if (length(pfms)==2){
        pfms <- DNAmotifAlignment(pfms)
        plotMotifLogoStack(pfms)
        return(invisible(list(phylog=NULL, pfms=pfms)))
    }
    
    ##calculate the distances
    dots <- list(...)
    if(!"phylog" %in% names(dots)){
      if(!psam){
        hc <- clusterMotifs(pfms)
        pfms <- pfms[hc$order]
        pfms <- DNAmotifAlignment(pfms)
        phylog <- hclust2phylog(hc)
      }
    }
    if(layout=="treeview" && !exists("hc")){
      if(!psam){
        hc <- clusterMotifs(pfms)
      }
    }
    switch(layout,
           stack = {
               plotMotifLogoStack(pfms, ...)
           },
           treeview = {
               if(!exists("hc")){
                 stop("hc is required for plotMotifLogoStackWithTree.")
               }
              plotMotifLogoStackWithTree(pfms, hc=hc, ...)
           },
           phylog = {
             if(!exists("phylog")){
               stop("phylog is required for plotMotifStackWithPhylog.")
             }
               plotMotifStackWithPhylog(phylog=phylog, pfms=pfms, ...)
           },
           radialPhylog = {
             if(!exists("phylog")){
               stop("phylog is required for plotMotifStackWithRadialPhylog.")
             }
               args <- list(phylog=phylog, pfms=pfms, ...)
               for(i in names(args)){
                   if(i %in% c("col.leaves", "col.leaves.bg", "col.bg",
                               "col.inner.label.circle",
                               "col.outer.label.circle")){
                       args[[i]] <- args[[i]][hc$order]
                   }
               }
               do.call(plotMotifStackWithRadialPhylog, args)
           },
           plotMotifLogoStack(pfms, ...)
    )
    
    return(invisible(list(phylog=phylog, pfms=pfms)))
}

