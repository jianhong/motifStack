###############################################################################
####### motifstack


#' plot a DNA sequence logo stack
#' 
#' Plot a DNA sequence logo stack
#' 
#' 
#' @param pfms a list of objects of class \link{pfm}
#' @param layout layout of the logo stack, stack, treeview or radialPhylog
#' @param reorder logical(1). Default TRUE. Set to FALSE will do alignment but
#' keep the order of the pfms. This parameter only work for stack layout.
#' @param \dots any parameters could to pass to \link{plotMotifLogoStack},
#' \link{plotMotifLogoStackWithTree}, \link{plotMotifStackWithPhylog} or
#' \link{plotMotifStackWithRadialPhylog}
#' @return return a list contains pfms and phylog
#' @export
#' @importFrom ade4 hclust2phylog
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
#'     pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#'     motifStack(pfms, "radialPhylog")
#'     
#'     ## AA motifs
#'     pcms<-importMatrix(system.file("extdata", "prot.meme", 
#'                                    package="motifStack"),
#'                    format="meme", to="pfm")
#'     motifStack(pcms[1:5])
#'     motifStack(pcms[1:5], reorder=FALSE)
#'   }
#' 
motifStack <-function(pfms, 
                      layout=c("stack", "treeview", "phylog", "radialPhylog"),
                      reorder=TRUE,
                      ...){
    if(!is.list(pfms)){
      if(is(pfms, "pcm")) pfms <- pcm2pfm(pfms)
        plot(pfms)
        return(invisible())
    }
    layout <- match.arg(layout)
    if(!pfms[[1]]@alphabet %in% c("AA")){
      if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) {
        pfms <- lapply(pfms, pcm2pfm)
      }
      if (any(unlist(lapply(pfms, function(.ele) 
        !inherits(.ele, c("pfm", "psam")))))) 
        stop("pfms must be a list of pfm, pcm or psam objects")
    }
    if(all(sapply(pfms, function(.ele) is(.ele, "psam")))){
      psam <- TRUE
    }else{
      psam <- FALSE
    }
    alphab <- vapply(pfms, function(.ele) .ele@alphabet=="RNA", 
                     FUN.VALUE = TRUE)
    alphab <- any(alphab)
    dots <- list(...)
    if("revcomp" %in% names(dots)){
      revcomp <- dots$revcomp
    }else{
      revcomp <- rep(!alphab, length(pfms))
      dots$revcomp <- revcomp[1]
    }
    if (length(pfms)<2)
        stop("length of pfms less than 2")
    if (length(pfms)==2){
      if(pfms[[1]]@alphabet %in% c("DNA", "RNA")) {
        pfms <- DNAmotifAlignment(pfms, revcomp = revcomp)
      }else{
        if(pfms[[1]]@alphabet %in% c("AA") && layout!="stack"){
          pfms <- AAmotifAlignment(pfms)
        }else{
          if(pfms[[1]]@alphabet=="others"){
            warning("alphabet is set to 'others'. ",
                    "Currently 'others' is not support for alignment.")
          }
        }
      }
      plotMotifLogoStack(pfms)
      return(invisible(list(phylog=NULL, pfms=pfms)))
    }
    dots$motifs <- pfms
    ##calculate the distances
    if(!"phylog" %in% names(dots)){
      if(!psam){
        hc <- do.call(clusterMotifs, args=dots)
        if(layout!="stack" && reorder) pfms <- pfms[hc$order]
        if(pfms[[1]]@alphabet %in% c("DNA", "RNA")) {
          pfms <- DNAmotifAlignment(pfms, revcomp = revcomp)
        }else{
          if(pfms[[1]]@alphabet %in% c("AA") && layout!="stack"){
            pfms <- AAmotifAlignment(pfms)
          }else{
            if(pfms[[1]]@alphabet=="others"){
              warning("alphabet is set to 'others'. ",
                      "Currently 'others' is not support for alignment.")
            }
          }
        }
        phylog <- hclust2phylog(hc)
      }
    }
    if(layout=="treeview" && !exists("hc")){
      if(!psam){
        dots$motifs <- pfms
        hc <- do.call(clusterMotifs, args=dots)
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

