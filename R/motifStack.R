###############################################################################
####### motifstack
motifStack <-function(pfms, 
                      layout=c("stack", "treeview", "phylog", "radialPhylog"), 
                      ...){
    if(!is.list(pfms)){
        plot(pfms)
        return(invisible())
    }
    layout <- match.arg(layout)
    if(all(sapply(pfms, class)=="pcm")) pfms <- lapply(pfms, pcm2pfm)
    if (any(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pfm"))))) 
        stop("pfms must be a list of pfm objects")
    if (length(pfms)<2)
        stop("length of pfms less than 2")
    if (length(pfms)==2){
        pfms <- DNAmotifAlignment(pfms)
        plotMotifLogoStack(pfms)
        return(invisible(list(phylog=NULL, pfms=pfms)))
    }
    pfmList2matrixList <- function(pfms){
        m <- lapply(pfms, pfm2pwm)
        names(m) <- unlist(lapply(pfms, function(.ele) .ele@name))
        m
    }
    
    ##calculate the distances
    dots <- list(...)
    if(!"phylog" %in% names(dots)){
        jaspar.scores <- readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
        d <- motifDistances(pfmList2matrixList(pfms), DBscores=jaspar.scores)
        hc <- motifHclust(d, method="average")
        pfms <- pfms[hc$order]
        pfms <- DNAmotifAlignment(pfms)
        phylog <- hclust2phylog(hc)
    }
    if(layout=="treeview" && !exists("hc")){
        jaspar.scores <- readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
        d <- motifDistances(pfmList2matrixList(pfms), DBscores=jaspar.scores)
        hc <- motifHclust(d, mothod="average")
    }
    switch(layout,
           stack = {
               plotMotifLogoStack(pfms, ...)
           },
           treeview = {
               plotMotifLogoStackWithTree(pfms, hc=hc, ...)
           },
           phylog = {
               plotMotifStackWithPhylog(phylog=phylog, pfms=pfms, ...)
           },
           radialPhylog = {
               args <- list(phylog=phylog, pfms=pfms, ...)
               for(i in names(args)){
                   if(i %in% c("col.leaves", "col.leaves.bg", "col.bg", "col.inner.label.circle", "col.outer.label.circle")){
                       args[[i]] <- args[[i]][hc$order]
                   }
               }
               do.call(plotMotifStackWithRadialPhylog, args)
           },
           plotMotifLogoStack(pfms, ...)
    )
    
    return(invisible(list(phylog=phylog, pfms=pfms)))
}

