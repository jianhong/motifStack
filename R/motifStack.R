###############################################################################
####### motifstack
motifStack <-function(pfms, 
                      layout=c("stack", "treeview", "phylog", "radialPhylog"), 
                      ...){
    if(!is.list(pfms)){
      if(is(pfms, "pcm")) pfms <- pcm2pfm(pfms)
        plot(pfms)
        return(invisible())
    }
    layout <- match.arg(layout)
    if(all(sapply(pfms, function(.ele) is(.ele, "pcm")))) pfms <- lapply(pfms, pcm2pfm)
    if (any(unlist(lapply(pfms, function(.ele) !inherits(.ele, c("pfm", "psam")))))) 
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
    pfmList2matrixList <- function(pfms){
        m <- lapply(pfms, pfm2pwm)
        names(m) <- unlist(lapply(pfms, function(.ele) .ele@name))
        m
    }
    
    ##calculate the distances
    dots <- list(...)
    if(!"phylog" %in% names(dots)){
      if(!psam){
        jaspar.scores <- readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
        d <- motifDistances(pfmList2matrixList(pfms), DBscores=jaspar.scores)
        hc <- motifHclust(d, method="average")
        pfms <- pfms[hc$order]
        pfms <- DNAmotifAlignment(pfms)
        phylog <- hclust2phylog(hc)
      }
    }
    if(layout=="treeview" && !exists("hc")){
      if(!psam){
        jaspar.scores <- readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
        d <- motifDistances(pfmList2matrixList(pfms), DBscores=jaspar.scores)
        hc <- motifHclust(d, mothod="average")
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

