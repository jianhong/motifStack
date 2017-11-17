motifPiles <- function (phylog, pfms=NULL, pfms2=NULL, 
                        r.tree=.45, col.tree=NULL,
                        cnodes=0, labels.nodes=names(phylog$nodes), clabel.nodes=0,
                        cleaves=.2, labels.leaves=names(phylog$leaves), clabel.leaves=1,
                        col.leaves=rep("black", length(labels.leaves)),
                        col.leaves.bg=NULL, col.leaves.bg.alpha=1,
                        r.pfms=NA, r.pfms2=NA, motifScale=c("logarithmic", "linear", "none"),
                        col.pfms=NULL, col.pfms.width=0.02,
                        col.pfms2=NULL, col.pfms2.width=0.02,
                        r.anno=0, col.anno=list(),
                        pfmNameSpliter=";", rcpostfix="(RC)", ic.scale=TRUE,
                        plotIndex=FALSE, IndexCol="black", IndexCex=.8,
                        groupDistance=NA, groupDistanceLineCol="red"){
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
    for(tobechecked in c("col.leaves", "col.leaves.bg", "col.tree")){
        if(checkLength(eval(as.symbol(tobechecked)))) stop(paste("the length of", tobechecked, "should be same as the length of leaves"))
        if(checkNA(eval(as.symbol(tobechecked)))) stop(paste("contain NA in", tobechecked))
    }
    if(length(r.anno)>length(col.anno)){
        if(r.anno!=0 & length(r.anno)!=1){
            stop("length of col.anno must be same as r.anno you input")
        }
    }
    if(length(col.anno)>0){
        if(class(col.anno)!="list") stop("col.anno must be a object of list")
        for(i in 1:length(col.anno)){
            if(checkLength(col.anno[[i]])) stop(paste("the length of col.anno[[", i, "]]should be same as the length of leaves"))
            if(checkNA(col.anno[[i]])) stop(paste("contain NA in col.anno[[", i, "]]"))
        }
    }
    motifScale <- match.arg(motifScale)
    leaves.names <- names(phylog$leaves)
    nodes.number <- length(phylog$nodes)
    nodes.names <- names(phylog$nodes)
    if (length(labels.leaves) != leaves.number)
        labels.leaves <- names(phylog$leaves)
    if (length(labels.nodes) != nodes.number)
        labels.nodes <- names(phylog$nodes)
    if (r.tree < 0.05 || r.tree > .5)
        stop("'r.tree': non convenient value [0.05, .5]")
    leaves.car <- gsub("[_]", " ", labels.leaves)
    nodes.car <- gsub("[_]", " ", labels.nodes)
    opar <- par(mar = par("mar"), srt = par("srt"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow=c(1,1))
    dis <- phylog$droot
    max_Dis <- max(dis)
    if(!is.na(groupDistance)){
        groupDistance <- (max_Dis - groupDistance) / max_Dis
    }
    dis <- dis/max_Dis
    dist.leaves <- dis[leaves.names]
    dist.nodes <- dis[nodes.names]
    xlim <- c(0, max_Dis/r.tree)
    asp <- c(1, 1)
    pin <- dev.size("in")
    asp <- pin[2L]/pin[1L]
    plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
                 yaxt = "n", xlim = xlim, ylim = c(-0, 1), xaxs = "i", 
                 yaxs = "i", frame.plot = FALSE)
    
    xcar <- max_Dis*1.05 + .01
    leftMar <- .01
    
    vpheight <- 1/(leaves.number+1)
    leaves.width <- lapply(leaves.car, strwidth, units="user", cex=clabel.leaves)
    string.width <- strwidth("M", units="figure", cex=IndexCex)
    r.leaves <- max(unlist(leaves.width))
    r.leaves <- r.leaves + xcar
    r.total <- (r.leaves + leftMar)/xlim[2] + col.pfms.width + string.width
    lenPFM <- c(FALSE, FALSE)
    if(!is.null(pfms)){
        if(is.na(r.pfms)){
            r.pfms <- r.total
        }
        pfmsLen <- sapply(pfms, function(.ele) length(unlist(strsplit(.ele@name, pfmNameSpliter))))
        if(motifScale=="linear"){
            vph <- vpheight * pfmsLen
        }else{
          if(motifScale=="logarithmic"){
            vph <- vpheight*(1+log2(pfmsLen))
          }else{
            vph <- rep(vpheight, length(pfmsLen))
          }
        }
        vpwidth <- mapply(function(.ele, vpht) vpht * ncol(.ele@mat) / 2 * asp, pfms, vph)
        r.total <- r.pfms + max(vpwidth) + col.pfms2.width + 2*string.width
        lenPFM[1] <- TRUE
    }
    if(!is.null(pfms2)){
        if(is.na(r.pfms2)){
            r.pfms2 <- r.total
        }
        pfmsLen2 <- sapply(pfms2, function(.ele) length(unlist(strsplit(.ele@name, pfmNameSpliter))))
        if(motifScale=="linear"){
            vph2 <- vpheight * pfmsLen2
        }else{
          if(motifScale=="logarithmic"){
            vph2 <- vpheight*(1+log2(pfmsLen2))
          }else{
            vph2 <- rep(vpheight, length(pfmsLen2))
          }
        }
        vpwidth2 <- mapply(function(.ele, vpht) vpht * ncol(.ele@mat) / 2 * asp, pfms2, vph2)
        r.total <- r.pfms2 + max(vpwidth2) + 2*string.width
        lenPFM[2] <- TRUE
    }
    
    if(sum(lenPFM)>0){
        plotIndex <- rep(plotIndex, 2)[lenPFM]
        IndexCol <- rep(IndexCol, 2)[lenPFM]
        IndexCex <- rep(IndexCex, 2)[lenPFM]
        names(plotIndex) <- names(IndexCol) <- names(IndexCex) <- c("pfm", "pfm2")[lenPFM]
    }
    
    xx <- c(dist.leaves, dist.nodes)
    
    r.total <- sum(c(r.total, r.anno), na.rm=TRUE)
    
    ratio <- .98/r.total
    if(ratio < 1){
        xcar <- xcar * ratio
        xx <- xx * ratio
        leftMar <- leftMar * ratio
        groupDistance <- groupDistance * ratio
        clabel.leaves <- clabel.leaves * ratio
        clabel.nodes <- clabel.nodes * ratio
        IndexCex <- IndexCex * ratio
        cleaves <- cleaves * ratio
        cnodes <- cleaves * ratio
        r.anno <- r.anno * ratio
        r.leaves <- r.leaves * ratio
        if(!is.null(pfms)) vpwidth <- vpwidth * ratio
        if(!is.null(pfms2)) vpwidth2 <- vpwidth2 * ratio
        if(!is.null(pfms)) r.pfms <- r.pfms * ratio
        if(!is.null(pfms2)) r.pfms2 <- r.pfms2 * ratio
        r.total <- r.total * ratio
        if(!is.null(col.pfms)) col.pfms.width <- col.pfms.width*ratio
        if(!is.null(col.pfms2)) col.pfms2.width <- col.pfms2.width*ratio
    }
    
    xx <- xx*(xlim[2]*r.tree) + leftMar
    groupDistance <- groupDistance*(xlim[2]*r.tree) + leftMar
    
    y <- (leaves.number:1)/(leaves.number + 1)
    names(y) <- leaves.names
    
    if (is.null(col.tree)){
        col.tree <- rep("black", leaves.number)
    }
    
    if (cleaves > 0) {
        for (i in 1:leaves.number) {
            points(xx[i], y[i], pch = 21, bg=1, cex = par("cex") * cleaves, col=col.tree[i])
        }
    }
    
    yn <- rep(0, nodes.number)
    names(yn) <- nodes.names
    yy <- c(y, yn)
    col.tree.lines <- col.tree
    names(col.tree.lines) <- names(phylog$leaves)
    for (i in 1:length(phylog$parts)) {
        w <- phylog$parts[[i]]
        but <- names(phylog$parts)[i]
        yy[but] <- mean(yy[w])
        colour <- floor(rowMeans(col2rgb(col.tree.lines[w])))
        col.tree.lines[but] <- rgb(colour[1], colour[2], colour[3], maxColorValue=255)
        b <- range(yy[w])
        segments(xx[but], b[1], xx[but], b[2], col=col.tree.lines[but])
        x1 <- xx[w]
        y1 <- yy[w]
        x2 <- rep(xx[but], length(w))
        segments(x1, y1, x2, y1, col=col.tree.lines[w])
    }
    if(!is.na(groupDistance)){
        if(length(groupDistanceLineCol)!=length(groupDistance)){
            groupDistanceLineCol <- rep(groupDistanceLineCol, ceiling(length(groupDistance)/length(groupDistanceLineCol)))[1:length(groupDistance)]
        }
        for(i in 1:length(groupDistance)){
            if(groupDistance[i] > 0)
                abline(v=groupDistance[i], col=groupDistanceLineCol[i], lty=2)
        }
    }
    if (cnodes > 0) {
        for (i in nodes.names) {
            points(xx[i], yy[i], pch = 21, bg="white", cex = cnodes)
        }
    }
    if (clabel.nodes > 0) {
        scatterutil.eti(xx[names(dist.nodes)], yy[names(dist.nodes)], nodes.car, 
                        clabel.nodes)
    }
    
    for (i in 1:leaves.number) {
        rect(xcar, y[i]-1/leaves.number/2, r.leaves, y[i]+1/leaves.number/2, col=highlightCol(col.leaves.bg[i], col.leaves.bg.alpha), border=NA)
        if(clabel.leaves>0) 
            text(xcar, y[i], leaves.car[i], adj = 0, cex = par("cex") * 
                     clabel.leaves, col=col.leaves[i])
        segments(xcar, y[i], xx[i], y[i], col = grey(0.7))
    }
    
    for(metaChar in c("\\","$","*","+",".","?","[","]","^","{","}","|","(",")"))
    {
        rcpostfix <- gsub(metaChar,paste("\\",metaChar,sep=""),rcpostfix,fixed=TRUE)
    }
    drawPFMcir <- function(pfms, r.x, ym, vpwid, vphei, col.pfms.wid, idx){
        assign("tmp_motifStack_symbolsCache", list(), envir=.globals)
        ##extract names
        pfmNames <- lapply(pfms, function(.ele) .ele@name)
        if(idx=="pfm2"){
          ## sort the pfms
          freq <- lengths(strsplit(unlist(pfmNames), pfmNameSpliter))
          freq.id <- order(freq, decreasing = TRUE)
          pfms <-pfms[freq.id]
          pfmNames <- pfmNames[freq.id]
          vpwid <- vpwid[freq.id]
          vphei <- vphei[freq.id]
        }
        for(i in 1:length(pfmNames)){
            pfmname <- unlist(strsplit(pfmNames[[i]], pfmNameSpliter))
            pfmname <- gsub(paste(rcpostfix,"$",sep=""),"",pfmname)
            pfmIdx <- which(makeLeaveNames(labels.leaves) %in% makeLeaveNames(pfmname))
            if(length(pfmIdx)==0) 
                pfmIdx <- which(makeLeaveNames(names(phylog$leaves)) 
                                %in% makeLeaveNames(pfmname))
            if(length(pfmIdx)>0){
                vpy <- mean(ym[pfmIdx])
                pushViewport(viewport(x=r.x, y=vpy, width=vpwid[i], height=vphei[i], just=c(0, .5)))
                plotMotifLogoA(pfms[[i]], ic.scale=ic.scale)
                popViewport()
                if(plotIndex[idx]&length(pfms)<leaves.number) {
                    grid.text(label=i, x=r.x-col.pfms.wid, 
                              y=vpy, 
                              gp=gpar(col=IndexCol[idx], cex=IndexCex[idx]), just="right")
                }
            }else{
                warning(paste("No leave named as ", paste(pfmname, collapse=", ")), sep="")
            }
        }
        rm(list="tmp_motifStack_symbolsCache", envir=.globals)
    }
    
    if(!is.null(col.pfms)){
        rect((r.pfms-col.pfms.width)*xlim[2], y-vpheight/2, r.pfms*xlim[2], y+vpheight/2, col=col.pfms, border=NA)
    }
    
    if(!is.null(pfms) && r.pfms < .95){
        drawPFMcir(pfms, r.pfms, y, vpwidth, vph, col.pfms.width, "pfm")
    }
    
    if(!is.null(col.pfms2)){
        rect((r.pfms2-col.pfms2.width)*xlim[2], y-vpheight/2, r.pfms2*xlim[2], y+vpheight/2, col=col.pfms2, border=NA)
    }
    
    if(!is.null(pfms2) && r.pfms2 < .95){
        drawPFMcir(pfms2, r.pfms2, y, vpwidth2, vph2, col.pfms2.width, "pfm2")
    }
    
    if(sum(r.anno, na.rm=TRUE)>0){
        r.anno.diff <- c(0, 0, cumsum(r.anno[-length(r.anno)]))
        r.anno <- r.total - sum(r.anno, na.rm=TRUE) + c(0, r.anno) + r.anno.diff
        for(i in 2:length(r.anno)){
            rect(r.anno[i-1]*xlim[2], y-vpheight/2, r.anno[i]*xlim[2], y+vpheight/2, col=col.anno[[i-1]], border=NA)
        }
    }
    
    return(invisible())
}