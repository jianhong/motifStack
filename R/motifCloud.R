#' plot a DNA sequence logo cloud
#' 
#' Plot a DNA sequence logo cloud
#' 
#' 
#' @param motifSig an object of class \linkS4class{motifSig}
#' @param rcpostfix postfix for reverse-complement motif names, default: (RC)
#' @param layout layout of the logo cloud, rectangles, cloud or tree
#' @param scale A vector of length 2 indicating the range of the size of the
#' sequence logo.
#' @param rot.per proportion sequence logo with 90 degree rotation. Only work
#' for "cloud" layout
#' @param draw.box draw box for each sequence logo or not
#' @param draw.freq label frequency of each signature or not
#' @param box.col color of box for each sequence logo
#' @param freq.col color of frequency label
#' @param group.col color setting for groups
#' @param groups a named vectors of motif groups
#' @param draw.legend draw group color legend or not
#' @param font font of logo
#' @param ic.scale logical If TRUE, the height of each column is proportional
#' to its information content. Otherwise, all columns have the same height.
#' @return none
#' @export
#' @importFrom grid grid.polygon gpar pushViewport viewport grid.rect grid.text
#' popViewport unit 
#' @importFrom graphics plot.window par plot.new strheight plot.default segments
#' legend 
#' @importFrom stats runif
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
#'     library(MotIV)
#'     jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), 
#'                                     "extdata", "jaspar2010_PCC_SWU.scores"))
#'     d <- MotIV::motifDistances(lapply(pfms, pfm2pwm))
#'     hc <- MotIV::motifHclust(d, method="average")
#'     library(ade4)
#'     phylog <- ade4::hclust2phylog(hc)
#'     leaves <- names(phylog$leaves)
#'     pfms <- pfms[leaves]
#'     pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#'     motifSig <- motifSignature(pfms, phylog, groupDistance=0.1)
#'     motifCloud(motifSig)
#'   }
#' 
motifCloud <- function(motifSig, rcpostfix="(RC)", 
                       layout=c("rectangles", "cloud", "tree"), 
                       scale=c(6, .5), rot.per=.1,
                       draw.box=TRUE, draw.freq=TRUE, 
                       box.col="gray", freq.col="gray",
                       group.col=NULL, groups=NULL, draw.legend=FALSE,
                       font="Helvetica-Bold", ic.scale=TRUE)
{
  if (!inherits(motifSig, "motifSig")) 
    stop("motifSig be object of motifSig. You could try\n
         ?motifSignature\nto get a motifSig.")
  signatures=signatures(motifSig)
  freq=frequence(motifSig)
  nodelist=nodelist(motifSig)
  if((!is.null(group.col)) && (!is.null(groups))){
    if(!all(names(table(groups)) %in% names(group.col)))
      stop("some value of groups has no corresponding color in group.col")
    names(groups) <- gsub(rcpostfix,"",names(groups), fixed=TRUE)
    pfmNames <-
      as.character(unlist(lapply(signatures, function(.ele) {
        unlist(strsplit(.ele@name, ";"))
        })))
    pfmNames <- gsub(rcpostfix,"",pfmNames, fixed=TRUE)
    if(!all(pfmNames %in% names(groups)))
      stop("not all pfms has a given group")
  }
  layout <- match.arg(layout, c("cloud","rectangles", "tree"), several.ok=TRUE)
  layout <- layout[1]
  normedFreq <- freq/max(freq)
  size <- (scale[1]-scale[2])*normedFreq + scale[2]
  grid.pie <- function(x, edges = 200, radius = 0.8, col = NULL, cex=1){
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
      stop("'x' values must be positive.")
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    if (is.null(col)) 
      col <- c("white", "lightblue", "mistyrose", "lightcyan", 
               "lavender", "cornsilk")
    col <- rep(col, length.out = nx)
    twopi <- 2 * pi
    t2xy <- function(t) {
      t2p <- twopi * t
      list(x = radius * cos(t2p)/2 + 0.5, y = radius * sin(t2p)/2 + 0.5)
    }
    for (i in 1L:nx) {
      n <- max(2, floor(edges * dx[i]))
      P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
      grid.polygon(c(P$x, 0.5), c(P$y, 0.5), gp=gpar(col = col[i],
                                                     fill = col[i]))
    }
  }
  plotSignature <- function(x, y, wid, ht, just, angle, sig, freq, normedFreq){
    pushViewport(viewport(x=x, y=y, 
                          width=wid, height=ht,
                          just=just, angle=angle))
    if(draw.box) grid.rect(gp=gpar(col=box.col, 
                                   lty="dashed", 
                                   fill="transparent"))
    plotMotifLogoA(sig, font=font, ic.scale=ic.scale)
    if(draw.freq) grid.text(label=freq,x=.95,y=.95,
                            gp=gpar(col=freq.col,cex=normedFreq), 
                            just=c("right","top"))
    if((!is.null(group.col)) & (!is.null(groups))) {
      sigNames <- unlist(strsplit(sig@name,";"))
      gps <- table(groups[sigNames])
      gp.col <- group.col[names(gps)]
      pushViewport(viewport(x=.1*as.numeric(ht)/as.numeric(wid), y=.9,
                            width=.2*as.numeric(ht)/as.numeric(wid), height=.2))
      grid.pie(gps, col=gp.col, cex=normedFreq)
      popViewport()
    }
    popViewport()
  }
  
  op <- par(mar=c(0,0,0,0))
  on.exit(par(op))
  
  if(layout=="rectangles"){
    plot.new()
    boxesInch <- list()
    minWid <- .Machine$integer.max
    minHt <- .Machine$integer.max
    for(i in 1:length(signatures)){
      ht <- strheight("ACGT", units="inches", cex=size[i])
      wid <- ht*ncol(signatures[[i]]@mat)/2
      boxesInch[[i]] <- new("Rect", width=wid, height=ht)
      if(ht < minHt) minHt <- ht
      if(wid < minWid) minWid <- wid
    }
    pin <- par("pin")
    bin <- new("Rect", width=pin[1], height=pin[2])
    bin.ratio <- bin@height / bin@width
    basic1.wd <- max(boxesInch[[1]]@width, boxesInch[[1]]@height)
    if(bin.ratio < 1) {
      basic1.ht <- basic1.wd * bin.ratio
    }else{
      basic1.ht <- basic1.wd
      basic1.wd <- basic1.ht/bin.ratio
    }
    basic1 <- new("Rect", width=basic1.wd, height=basic1.ht)
    freeRect <- list(basic1)
    bestNodelist <- list()
    maxX <- 0
    maxY <- 0
    for(i in 1:length(signatures)){
      if(length(freeRect)<1){
        if(interactive()) {
          warning("signature could not be fit on page. It will not be plotted.")
        }
        break
      }
      getBestNode <- function(freeRect, boxesInch, i, bin, basic1){
        bestNode <- new("Rect")
        bestX <- basic1@width
        bestY <- basic1@height
        for(j in 1:length(freeRect)){
          if(freeRect[[j]]@width >= boxesInch[[i]]@width && 
             freeRect[[j]]@height >= boxesInch[[i]]@height){
            topSideY <- freeRect[[j]]@y + boxesInch[[i]]@height
            if(topSideY < bestY || 
               (topSideY == bestY && freeRect[[j]]@x < bestX)){
              bestNode@x <- freeRect[[j]]@x
              bestNode@y <- freeRect[[j]]@y
              bestNode@width <- boxesInch[[i]]@width
              bestNode@height <- boxesInch[[i]]@height
              bestY <- topSideY
              bestX <- freeRect[[j]]@x
            }
          }
          if(freeRect[[j]]@width >= boxesInch[[i]]@height && 
             freeRect[[j]]@height >= boxesInch[[i]]@width){
            topSideY <- freeRect[[j]]@y + boxesInch[[i]]@width
            if(topSideY < bestY || 
               (topSideY == bestY && freeRect[[j]]@x < bestX)){
              bestNode@x <- freeRect[[j]]@x
              bestNode@y <- freeRect[[j]]@y
              bestNode@width <- boxesInch[[i]]@height
              bestNode@height <- boxesInch[[i]]@width
              bestY <- topSideY
              bestX <- freeRect[[j]]@x
            }
          }
        }
        if(bestNode@height==0){
          basic2.wd <- min(boxesInch[[i]]@width, boxesInch[[i]]@height)
          basic2.ht <- basic2.wd * bin.ratio
          if(bin.ratio < 1) {
            basic2.ht <- basic2.wd * bin.ratio
          }else{
            basic2.ht <- basic2.wd
            basic2.wd <- basic2.ht/bin.ratio
          }
          basic2 <- new("Rect", width=basic1@width+basic2.wd, 
                        height=basic1@height+basic2.ht)
          freeRect <- lapply(freeRect, function(.ele){
            if(.ele@x+.ele@width==basic1@width) {
              .ele@width <- .ele@width + basic2.wd
            }
            if(.ele@y+.ele@height==basic1@height) {
              .ele@height <- .ele@height + basic2.ht
            }
            .ele
          })
          freeRect <- 
            c(freeRect, new("Rect", x=basic1@width, y=0, 
                            width=basic2.wd, height=basic2@height), 
              new("Rect", x=0, y=basic1@height, 
                  width=basic2@width, height=basic2.ht))
          basic1 <- basic2
          bestNodeRes <- Recall(freeRect, boxesInch, i, bin, basic1)
          bestNode <- bestNodeRes[["bestNode"]]
          freeRect <- bestNodeRes[["freeRect"]]
          basic1 <- bestNodeRes[["basic1"]]
        }
        list(bestNode=bestNode, freeRect=freeRect, basic1=basic1)
      }
      bestNodeRes <- getBestNode(freeRect, boxesInch, i, bin, basic1)
      bestNode <- bestNodeRes[["bestNode"]]
      freeRect <- bestNodeRes[["freeRect"]]
      basic1 <- bestNodeRes[["basic1"]]
      if(bestNode@height==0){
        if(interactive()) {
          warning("signature could not be fit on page. It will not be plotted.")
        }
        break
      }
      tobeAppend <- list()
      leftNodes <- lapply(freeRect, function(freeNode, usedNode){
        if(usedNode@x >= freeNode@x + freeNode@width || 
           usedNode@x + usedNode@width <= freeNode@x ||
           usedNode@y >= freeNode@y + freeNode@height || 
           usedNode@y + usedNode@height <= freeNode@y){
          TRUE
        }else{
          if(usedNode@x < freeNode@x + freeNode@width &&
             usedNode@x + usedNode@width > freeNode@x){
            ## New node at the top side of the used node
            if(usedNode@y > freeNode@y && 
               usedNode@y <freeNode@y + freeNode@height){
              newFreeNode <- freeNode
              newFreeNode@height <- usedNode@y - newFreeNode@y
              tobeAppend <<- c(tobeAppend, newFreeNode)
            }
            ## New node at the bottom side of the used node
            if(usedNode@y + usedNode@height < freeNode@y + freeNode@height){
              newFreeNode <- freeNode
              newFreeNode@y <- usedNode@y + usedNode@height
              newFreeNode@height <- freeNode@y + freeNode@height - 
                (usedNode@y + usedNode@height)
              tobeAppend <<- c(tobeAppend, newFreeNode)
            }
          }
          if(usedNode@y < freeNode@y + freeNode@height && 
             usedNode@y + usedNode@height > freeNode@y){
            ## New node at the left side of the used node
            if(usedNode@x > freeNode@x && 
               usedNode@x <freeNode@x + freeNode@width){
              newFreeNode <- freeNode
              newFreeNode@width <- usedNode@x - newFreeNode@x
              tobeAppend <<- c(tobeAppend, newFreeNode)
            }
            ## New node at the right side of the used node
            if(usedNode@x + usedNode@width < freeNode@x + freeNode@width){
              newFreeNode <- freeNode
              newFreeNode@x <- usedNode@x + usedNode@width
              newFreeNode@width <- freeNode@x + freeNode@width -
                (usedNode@x + usedNode@width)
              tobeAppend <<- c(tobeAppend, newFreeNode)
            }
          }
          FALSE
        }
      }, bestNode)     
      freeRect <- c(freeRect[unlist(leftNodes)], tobeAppend)
      ## prune freeRect
      if(length(freeRect)>1){
        leftNodes <- rep(TRUE, length(freeRect))
        for(m in 1:(length(freeRect)-1)){
          if(leftNodes[m]){
            for(n in (m+1):length(freeRect)){
              if(leftNodes[n]){
                if(isContainedIn(freeRect[[m]], freeRect[[n]])){
                  leftNodes[m] <- FALSE
                  break;
                }
                if(isContainedIn(freeRect[[n]], freeRect[[m]])){
                  leftNodes[n] <- FALSE
                }
              }
            }
          }
        }
        freeRect <- freeRect[leftNodes]
      }
      
      ##push bestNode
      bestNodelist <- c(bestNodelist, bestNode)
      maxX <- max(unlist(lapply(freeRect, function(.ele){.ele@x+.ele@width})))
      maxY <- max(unlist(lapply(freeRect, function(.ele){.ele@y+.ele@height})))
    }
    ratio <- min(bin@width/maxX, bin@height/maxY)
    for(i in 1:length(bestNodelist)){
      bestNode <- bestNodelist[[i]]
      ## draw logo in bestNode
      if(bestNode@width<bestNode@height){
        plotSignature(x=unit(ratio*bestNode@x, "inches"), 
                      y=unit(ratio*bestNode@y+ratio*bestNode@height, "inches"),
                      wid=unit(ratio*bestNode@height, "inches"), 
                      ht=unit(ratio*bestNode@width, "inches"), 
                      just=c("left","bottom"), 
                      angle=-90,
                      sig=signatures[[i]], 
                      freq=freq[i], 
                      normedFreq=normedFreq[i])
      }else{
        plotSignature(x=unit(ratio*bestNode@x, "inches"), 
                      y=unit(ratio*bestNode@y, "inches"), 
                      wid=unit(ratio*bestNode@width, "inches"), 
                      ht=unit(ratio*bestNode@height, "inches"), 
                      just=c("left","bottom"),
                      angle=0,
                      sig=signatures[[i]], 
                      freq=freq[i], 
                      normedFreq=normedFreq[i])
      }
    }
  }else if(layout=="cloud"){
    plot.new()
    last <- 1
    overlap <- function(x1, y1, sw1, sh1){
      s <- 0
      if(length(boxes)==0) return(FALSE)
      for(i in c(last, 1:length(boxes))){
        bnds <- boxes[[i]]
        x2 <- bnds[1]
        y2 <- bnds[2]
        sw2 <- bnds[3]
        sh2 <- bnds[4]
        if(x1 < x2) overlap <- x1+sw1 > x2-s
        else overlap <- x2 +sw2 > x1-s
        if(y1 <y2) overlap <- overlap && (y1+sh1 > y2-s)
        else overlap <- overlap && (y2 +sh2 > y1-s)
        if(overlap){
          last <- i
          return(TRUE)
        }
      }
      FALSE
    }
    thetaStep <- .1
    rStep <- .05
    plot.window(c(0,1), c(0,1), asp=1)
    boxes <- list()
    
    for(i in 1:length(signatures)){
      sig <- signatures[[i]]
      rotWord <- runif(1)<rot.per
      r <- 0
      theta <- runif(1, 0, 2*pi)
      x1 <- .5
      y1 <- .5
      ht <- strheight("ACGT", cex=size[i])
      wid <- ht * ncol(sig@mat) / 2
      if(rotWord){
        tmp <- ht
        ht <- wid
        wid <- tmp
      }
      isOverlapped <- TRUE
      while(isOverlapped){
        if(!overlap(x1-.5*wid, y1-.5*ht, wid, ht) &&
           x1-.5*wid>0 && y1-.5*ht>0 &&
           x1+.5*wid<1 && y1+.5*ht<1){
          if(rotWord){
            plotSignature(x1, y1, ht, wid, "center", 90, sig, 
                          freq[i], normedFreq[i])
          }else{
            plotSignature(x1, y1, wid, ht, "center", 0, sig, 
                          freq[i], normedFreq[i])
          }
          boxes[[length(boxes)+1]] <- c(x1-.5*wid, y1-.5*ht, wid, ht)
          isOverlapped <- FALSE
        }else{
          if(r>sqrt(.5)){
            if(interactive()) {
              warning("signature could not be fit on page. 
                      It will not be plotted.")
              }
            isOverlapped <- FALSE
          }
          theta <- theta+thetaStep
          r <- r+rStep*thetaStep/(2*pi)
          x1 <- .5+r*cos(theta)
          y1 <- .5+r*sin(theta)
        }
      }
    }
  }else{
    #trim nodelist
    #get all path for each signature
    sigNames <- as.character(unlist(lapply(signatures, 
                                           function(.ele) .ele@name)))
    sigNames <- gsub(rcpostfix, "", sigNames, fixed=TRUE)
    paths <- lapply(sigNames, function(.ele, nodelist){
      pt <- c()
      getPath <- function(nodelist, nodename, pt){
        parent <- getParentNode(nodelist, nodename)
        if(!is.null(parent[1])){
          pt <- c(pt, parent[1])
          pt <- Recall(nodelist, parent[1], pt)
        }
        pt
      }
      pt <- getPath(nodelist, .ele, pt)
      pt
    }, nodelist)
    #leave the nodes in the paths
    paths <- unique(unlist(paths))
    nodelist <- nodelist[paths]
    
    last <- 1
    overlap <- function(x1, y1, sw1, sh1){
      s <- 0
      if(length(boxes)==0) return(FALSE)
      for(i in c(last, 1:length(boxes))){
        bnds <- boxes[[i]]
        x2 <- bnds[1]
        y2 <- bnds[2]
        sw2 <- bnds[3]
        sh2 <- bnds[4]
        if(x1 < x2) overlap <- x1+sw1 > x2-s
        else overlap <- x2 +sw2 > x1-s
        if(y1 <y2) overlap <- overlap && (y1+sh1 > y2-s)
        else overlap <- overlap && (y2 +sh2 > y1-s)
        if(overlap){
          last <- i
          return(TRUE)
        }
      }
      FALSE
    }
    boxes <- list()
    step <- .1
    plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
                 xaxt = "n", yaxt = "n", xlim = c(0, 1), ylim = c(0, 1), 
                 xaxs = "i", yaxs = "i", frame.plot = FALSE)
    getRdist <- function(nodelist, Rdist=c(), nodename="Root", distp=0){
      currNode <- nodelist[[nodename]]
      distl <- distp + currNode@distl
      distr <- distp + currNode@distr
      currDist <- c(distl, distr)
      names(currDist) <- c(currNode@left, currNode@right)
      Rdist <- c(Rdist, currDist)
      if(!is.null(nodelist[[currNode@left]])) {
        Rdist <- getRdist(nodelist, Rdist, currNode@left, distp=distl)
      }
      if(!is.null(nodelist[[currNode@right]])) {
        Rdist <- getRdist(nodelist, Rdist, currNode@right, distp=distr)
      }
      Rdist
    }
    Rdist <- getRdist(nodelist)
    longestRdist <- 4*max(Rdist)
    plotUnrootedTree <- 
      function(signatures, nodelist, positions, nodename="Root", 
               AXIS=0, ANGLE=2*pi, orignalX=0.5, orignalY=0.5, longestRdist=0){
      currNode <- nodelist[[nodename]]
      RdistTotal <- 0
      if(currNode@left %in% c(paths, sigNames))
        RdistTotal <- currNode@distl
      if(currNode@right %in% c(paths, sigNames)) 
        RdistTotal <- RdistTotal + currNode@distr
      #plot left node
      start <- AXIS - ANGLE/2
      alpha <- 0
      if(currNode@left %in% c(paths, sigNames)){
        alpha <- ANGLE*currNode@distl/RdistTotal
        beta <- start + alpha/2
        x <- currNode@distl/longestRdist * cos(beta) + orignalX
        y <- currNode@distl/longestRdist * sin(beta) + orignalY
        if(!is.null(nodelist[[currNode@left]])) {
          segments(orignalX, orignalY, x, y)
          positions <- plotUnrootedTree(signatures, nodelist, positions, 
                                        nodename=currNode@left, AXIS=beta,
                                        ANGLE=alpha, orignalX=x, orignalY=y, 
                                        longestRdist=longestRdist)
        }else{
          pfmlid <- getPFMid(signatures, currNode@left)
          if(length(pfmlid)>0){
            ht <- strheight("ACGT", cex=size[pfmlid])
            wid <- ht * ncol(signatures[[pfmlid]]@mat) / 2
            segments(orignalX, orignalY, x, y)
            cp <- new("Rect", x=x, y=y, width=wid, height=ht)
            positions <- c(positions, new("pos", box=cp, beta=beta, 
                                          sig=signatures[[pfmlid]], 
                                          freq=freq[[pfmlid]],
                                          norm=normedFreq[pfmlid]))
          }
        }
      }
      #plot right node
      if(currNode@right %in% c(paths, sigNames)){
        start <- start + alpha
        alpha <- ANGLE*currNode@distr/RdistTotal
        beta <- start + alpha/2
        x <- currNode@distr/longestRdist * cos(beta) + orignalX
        y <- currNode@distr/longestRdist * sin(beta) + orignalY
        if(!is.null(nodelist[[currNode@right]])) {
          segments(orignalX, orignalY, x, y)
          positions <- plotUnrootedTree(signatures, nodelist, positions,
                                        nodename=currNode@right, AXIS=beta, 
                                        ANGLE=alpha, orignalX=x, orignalY=y, 
                                        longestRdist=longestRdist)
        }else{
          pfmrid <- getPFMid(signatures, currNode@right)
          if(length(pfmrid)>0){
            ht <- strheight("ACGT", cex=size[pfmrid])
            wid <- ht * ncol(signatures[[pfmrid]]@mat) / 2
            segments(orignalX, orignalY, x, y)
            cp <- new("Rect", x=x, y=y, width=wid, height=ht)
            positions <- c(positions, new("pos", box=cp, beta=beta, 
                                          sig=signatures[[pfmrid]], 
                                          freq=freq[[pfmrid]],
                                          norm=normedFreq[pfmrid]))
          }
        }
      }
      positions
    }
    positions <- plotUnrootedTree(signatures, nodelist, 
                                  positions=list(), longestRdist=longestRdist)
    positions <- positions[order(unlist(lapply(positions, 
                                               function(.ele) .ele@freq)), 
                                 decreasing=TRUE)]
    for(i in 1:length(positions)){
      isOverlapped <- TRUE
      thisMotif <- positions[[i]]
      x <- thisMotif@box@x
      y <- thisMotif@box@y
      wid <- thisMotif@box@width
      ht <- thisMotif@box@height
      beta <- thisMotif@beta
      x1 <- x
      y1 <- y
      while(isOverlapped){
        if(!overlap(x1-.5*wid, y1-.5*ht, wid, ht)){
          segments(x,y,x1,y1,lty="dashed", col="gray")
          #plot signature
          plotSignature(x1, y1, wid, ht, "center", 0,
                        thisMotif@sig, thisMotif@freq, thisMotif@norm)
          boxes[[length(boxes)+1]] <- c(x1-.5*wid, y1-.5*ht, wid, ht)
          isOverlapped <- FALSE
        }else{
          if((x1>=1-.5*wid || x1<=.5*wid) && (y1<=.5*ht || y1>=1-.5*ht)){
            if(interactive()) {
              warning("signature could not be fit on page. 
                      It will not be plotted.")
            }
            isOverlapped <- FALSE
          }
          r <- sqrt(wid*wid + ht*ht)*step
          if(((beta %% pi) - pi/2) < 0.0001) beta <- beta + pi/90
          x1 <- x1+r*cos(beta)
          y1 <- y1+r*sin(beta)
          if(x1-.5*wid<0) x1 <- .5*wid
          if(x1+.5*wid>1) x1 <-1-.5*wid
          if(y1-.5*ht<0) y1 <- .5*ht
          if(y1+.5*ht>1) y1 <-1-.5*ht
        }
      }
    }
  }
  #plot group legend
  if((!is.null(group.col)) && (!is.null(groups)) && draw.legend){
    legend("topright", names(group.col), fill=group.col, cex=0.5)
  }
  return(invisible())
}
