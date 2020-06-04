#' get signatures from motifs
#' 
#' extract signatures from multiple motifs by distance calculated from STAMP
#' 
#' 
#' @param pfms a list of objects of class pfm
#' @param phylog an object of class phylog
#' @param groupDistance maxmal distance of motifs in the same group
#' @param rcpostfix postfix for reverse-complement motif names, default: (RC)
#' @param min.freq signatures with frequency below min.freq will not be plotted
#' @param trim minimal information content for each position of signature
#' @param families for each family, the motif number in one signature should
#' only count as 1
#' @param sort sort the signatures by frequency or not.
#' @return an Object of class \linkS4class{motifSig}
#' @export
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
#'     jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), 
#'                                    "extdata", "jaspar2010_PCC_SWU.scores"))
#'     d <- MotIV::motifDistances(lapply(pfms, pfm2pwm))
#'     hc <- MotIV::motifHclust(d, method="average")
#'     phylog <- hclust2phylog(hc)
#'     leaves <- names(phylog$leaves)
#'     pfms <- pfms[leaves]
#'     pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#'     motifSig <- motifSignature(pfms, phylog, groupDistance=0.1)
#'   }
#' 
motifSignature <- function(pfms, phylog, groupDistance, rcpostfix="(RC)", 
                           min.freq=2, trim=0.2, families=list(), sort=TRUE){
  if (!inherits(phylog, "phylog")) 
    stop("phylog should be an object of phylog of package ade4")
  leaves.number <- length(phylog$leaves)
  if (length(pfms)!=leaves.number)
    stop("length of pfms and leaves of phylog should be identical.")
  if (all(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pfm"))))){
    pfms.class="pcm"
    if(any(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pcm")))))
      stop("pfms should be a list of objects of pfm or pcm")
  }else{
    pfms.class="pfm"
    if(any(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pfm")))))
      stop("pfms should be a list of objects of pfm or pcm")
  }
  if (missing(groupDistance))
    groupDistance <- max(phylog$droot)/10
  
  getSignature <- function(pfms, pfmnames, rcpostfix, pfms.class){
    .pfmnames <- unlist(strsplit(pfmnames, ";"))
    .pfms <- lapply(.pfmnames, function(.name, pfms) {
      pfms[[getPFMid(pfms, .name, rcpostfix=rcpostfix)]]
      }, pfms)
    if(length(.pfms)>1){
      .pfms <- DNAmotifAlignment(.pfms, rcpostfix="")
      .rown <- rownames(.pfms[[1]]@mat)
      .mats <- lapply(.rown, function(i, .pfms){
        do.call(rbind,lapply(.pfms, function(.ele, i) 
          .ele@mat[i, , drop=FALSE], i))
      }, .pfms)
      if(pfms.class=="pfms"){
        .mats <- do.call(rbind, lapply(.mats, colMeans))
      }
      else{
        .mats <- do.call(rbind, lapply(.mats, colSums))
        .mats <- pcm2pfm(.mats)
      }
      rownames(.mats) <- .rown
      colnames(.mats) <- 1:ncol(.mats)
    }else{
      if(pfms.class=="pfms") .mats <- .pfms[[1]]@mat
      else .mats <- pcm2pfm(.pfms[[1]]@mat)
    }
    new("pfm", mat=.mats, name=pfmnames, alphabet=.pfms[[1]]@alphabet,
        color=.pfms[[1]]@color, background=.pfms[[1]]@background)
  }
  
  if(groupDistance > max(phylog$droot)){
    signatures <- list(getSignature(pfms, names(phylog$leaves), 
                                    rcpostfix, pfms.class=pfms.class))
    signatures <- lapply(signatures, trimMotif, trim)
    nodelist <- list(Root=new("ouNode", 
                              left=paste(names(phylog$leaves),sep=";"), 
                              parent="Root", distl=max(phylog$droot), 
                              sizel=length(pfms), sizer=0))
    return(new("motifSig", signatures=signatures, 
               freq=length(pfms), 
               nodelist=nodelist, 
               gpcol=rep("gray30", length(pfms))))
  }
  #generate the new phylog from newick tree
  tree <- phylog$tre
  droot <- phylog$droot
  tree <- gsub(";$","",tree)
  str <- unlist(strsplit(gsub("^\\((.*?)\\)([^\\)\\(;,]+)$",
                              "\\1///\\2",tree),"///",fixed=TRUE))[1]
  nodelist <- list()
  getNodelist <- function(str, nodelist=list(), pid=0){
    ge <- gregexpr("\\(([^\\(\\),]+),([^\\(\\),]+)\\)([^\\(\\),]+)", str)
    if(ge[[1]][1]!=-1){## (A,B)C
      start <- ge[[1]]
      stop <- attr(ge[[1]],"match.length")
      trec<-c()
      for(i in 1:length(start)){
        trec<- c(trec,substr(str, start[i],start[i]+stop[i]-1))
      }
      nodes <- do.call(rbind,strsplit(trec,"\\)|\\(|,"))[,2:4,drop=FALSE]
      newnodes <- apply(nodes, 1, function(.ele) 
        new("ouNode",left=.ele[1], right=.ele[2], 
            parent=.ele[3], sizel=0, sizer=0))
      names(newnodes) <- nodes[,3]
      nodelist <- c(newnodes,nodelist)
      for(i in 1:length(trec)){
        str<-gsub(trec[i],gsub(".*?\\)(.*?)$","\\1",trec[i]),str,fixed=TRUE)
      }
      Recall(str=str, nodelist=nodelist, pid=pid)
    }else{
      ge <- gregexpr("\\((([^\\(\\),]+,){2,})([^\\(\\),]+)\\)([^\\(\\),]+)", 
                     str)
      if(ge[[1]][1]!=-1){ ##(A,B,C,...)N
        start <- ge[[1]]
        stop <- attr(ge[[1]],"match.length")
        trec<-c()
        for(i in 1:length(start)){
          trec<- c(trec,substr(str, start[i],start[i]+stop[i]-1))
        }
        for(i in 1:length(trec)){
          nodes <- strsplit(trec[i], "\\)|\\(|,")[[1]][-1]
          pnodes <- nodes[length(nodes)]
          nodes <- nodes[-length(nodes)]
          curNode <- nodes[length(nodes)]
          for(j in (length(nodes)-1):2){
            parNode <- paste("MSP", pid, sep="")
            newnodes <- new("ouNode", left=nodes[i], right=curNode,
                            parent=parNode, distl=0, distr=0,
                            sizel=0, sizer=0)
            curNode <- parNode
            pid <- pid + 1
            nodelist <- c(newnodes, nodelist)
            names(nodelist)[1] <- parNode
          }
          newnodes <- new("ouNode", left=nodes[1], right=curNode, 
                          parent=pnodes, sizel=0, sizer=0)
          nodelist <- c(newnodes, nodelist)
          names(nodelist)[1] <- pnodes
          str<-gsub(trec[i],gsub(".*?\\)(.*?)$","\\1",trec[i]),str,fixed=TRUE)
        }
        Recall(str=str, nodelist=nodelist, pid=pid)
      }else{ ##to Root one
        nodes <- unlist(strsplit(str, ","))
        if(length(nodes)>2){
          curNode <- nodes[length(nodes)]
          for(i in (length(nodes)-1):2){
            parNode <- paste("MSP", pid, sep="")
            newnodes <- new("ouNode", left=nodes[i], right=curNode,
                            parent=parNode, distl=0, distr=0, sizel=0, sizer=0)
            curNode <- parNode
            pid <- pid + 1
            nodelist <- c(newnodes, nodelist)
            names(nodelist)[1] <- parNode
          }
          Root <- new("ouNode", left=nodes[1], right=curNode, 
                      parent="Root", sizel=0, sizer=0)
          nodelist <- c(Root=Root, nodelist)
          return(nodelist)
        }else{
          Root <- new("ouNode", left=nodes[1], right=nodes[2], 
                      parent="Root", sizel=0, sizer=0)
          nodelist <- c(Root=Root, nodelist)
          return(nodelist)
        }
      }
    }
  }
  nodelist <- getNodelist(str)
  leaves <- names(phylog$leaves)
  buildTree <- function(nodelist, nodename, droot, dpar=0){
    nodelist[[nodename]]@distl <<- droot[nodelist[[nodename]]@left] - dpar
    nodelist[[nodename]]@distr <<- droot[nodelist[[nodename]]@right] - dpar
    if(nodelist[[nodename]]@left %in% leaves) nodelist[[nodename]]@sizel <<- 1
    if(nodelist[[nodename]]@right %in% leaves) nodelist[[nodename]]@sizer <<- 1
    if(!is.null(nodelist[[nodelist[[nodename]]@left]])){
      buildTree(nodelist, nodelist[[nodename]]@left, droot, 
                droot[nodelist[[nodename]]@left])
    }
    if(!is.null(nodelist[[nodelist[[nodename]]@right]])){
      buildTree(nodelist, nodelist[[nodename]]@right, droot,
                droot[nodelist[[nodename]]@right])
    }
  }
  buildTree(nodelist, "Root", droot, dpar=0)
  
  mergeNodes <- function(nodelist, leaves, groupDistance){
    l <- length(nodelist)
    for(i in 1:l){
      nodename <- names(nodelist)[i]
      currNode <- nodelist[[nodename]]
      if(currNode@left %in% leaves && currNode@right %in% leaves){
        if(currNode@distl < groupDistance && currNode@distr < groupDistance){
          parNodeInfo <- getParentNode(nodelist, nodename)
          if(!is.null(parNodeInfo)[1]){
            if(parNodeInfo[2]=="left") {
              nodelist[[parNodeInfo[1]]]@sizel <- 
                nodelist[[parNodeInfo[1]]]@sizel + max(currNode@sizel, 1) +
                max(currNode@sizer, 1)
              nodelist[[parNodeInfo[1]]]@distl <- 
                nodelist[[parNodeInfo[1]]]@distl + 
                (currNode@distl + currNode@distr)/2
              nodelist[[parNodeInfo[1]]]@left <- 
                paste(currNode@left, currNode@right, sep=";")
              leaves <- c(leaves[!leaves %in% c(currNode@left,
                                                currNode@right)], 
                          paste(currNode@left, currNode@right, sep=";"))
              nodelist[[nodename]] <- new("ouNode",left="NULL", right="NULL", 
                                          parent="NULL", sizel=0, sizer=0)
            }
            else if(parNodeInfo[2]=="right") {
              nodelist[[parNodeInfo[1]]]@sizer <- 
                nodelist[[parNodeInfo[1]]]@sizer + max(currNode@sizel, 1) + 
                max(currNode@sizer, 1)
              nodelist[[parNodeInfo[1]]]@distr <- 
                nodelist[[parNodeInfo[1]]]@distr + 
                (currNode@distl + currNode@distr)/2
              nodelist[[parNodeInfo[1]]]@right <- 
                paste(currNode@left, currNode@right, sep=";")
              leaves <- c(leaves[!leaves %in% 
                                   c(currNode@left, currNode@right)], 
                          paste(currNode@left, currNode@right, sep=";"))
              nodelist[[nodename]] <- 
                new("ouNode",left="NULL", right="NULL",
                    parent="NULL", sizel=0, sizer=0)
            }
          }
        }
      }
    }
    sel <- unlist(lapply(nodelist, function(.ele) .ele@parent=="NULL"))
    nodelist <- nodelist[!sel]
    nodelist <<- nodelist
    leaves <<- leaves
    if(length(nodelist)<l) mergeNodes(nodelist, leaves, groupDistance)
  }
  
  mergeNodes(nodelist, leaves, groupDistance)
  #get signatures and frequences (count of frequences should > min.freq)
  filterFamilies <- function(pfmnames, size, families){
    if(length(families)>0){
      pfmnames <- unlist(strsplit(pfmnames, ";"))
      if(length(pfmnames)!=size) stop("pfm names contain ';'")
      for(j in 1:length(families)){
        inters <- intersect(families[[j]], pfmnames)
        if(length(inters)>1){
          pfmnames <- pfmnames[!pfmnames %in% inters[2:length(inters)]]
        }
      }
      size <- length(pfmnames)
      pfmnames <- paste(pfmnames, sep="", collapse=";")
    }
    list(name=pfmnames,size=size)
  }
  signatures <- list()
  freq <- c()
  for(i in 1:length(nodelist)){
    if(nodelist[[i]]@sizel>=min.freq){
      ftmp <- filterFamilies(nodelist[[i]]@left, nodelist[[i]]@sizel, families)
      nodelist[[i]]@left <- ftmp[["name"]]
      nodelist[[i]]@sizel <- ftmp[["size"]]
      if(ftmp[["size"]]>=min.freq){
        signatures <- 
          c(signatures, getSignature(pfms, nodelist[[i]]@left,
                                     rcpostfix, pfms.class=pfms.class))
        freq <- c(freq, nodelist[[i]]@sizel)
      }
    }
    if(nodelist[[i]]@sizer>=min.freq){
      ftmp <- filterFamilies(nodelist[[i]]@right, nodelist[[i]]@sizer, families)
      nodelist[[i]]@right <- ftmp[["name"]]
      nodelist[[i]]@sizer <- ftmp[["size"]]
      if(ftmp[["size"]]>=min.freq){
        signatures <- c(signatures, 
                        getSignature(pfms, nodelist[[i]]@right,
                                     rcpostfix, pfms.class=pfms.class))
        freq <- c(freq, nodelist[[i]]@sizer)
      }
    }
  }
  
  #trime signatures
  signatures <- lapply(signatures, trimMotif, trim)
  ord <- unlist(lapply(signatures, function(.ele){inherits(.ele,"pfm")}))
  signatures <- signatures[ord]
  freq <- freq[ord]
  #sort signatures
  if(sort){
    ord <- order(freq, decreasing=TRUE)
    signatures <- signatures[ord]
    freq <- freq[ord]
  }
  if(length(freq)==0){
    if(interactive()) warning("All frequency are smaller than min.freq.")
    return(FALSE)
  }
  
  getGpCol <- function(sig, phylog){
    pfmNames <- lapply(sig, function(.ele){unlist(strsplit(.ele@name, ";"))})
    pfmNames <- mapply(function(.ele, .name){cbind(.ele, .name)},
                       pfmNames, paste("gps", 1:length(pfmNames), sep=""),
                       SIMPLIFY=FALSE)
    pfmNames <- do.call(rbind, pfmNames)
    pfmNames <- pfmNames[match(names(phylog$leaves), pfmNames[,1]),]
    bd.color <- c("gray80","gray30")
    in.color <- rle(pfmNames[,2])
    in.color$values <- bd.color[rep(1:2, 
                                    length(in.color$lengths))[
                                      seq_along(in.color$lengths)]]
    return(inverse.rle(in.color))
  }
  gpcol <- getGpCol(signatures, phylog)
  
  return(new("motifSig", signatures=signatures, freq=freq, 
             nodelist=nodelist, gpcol=gpcol))
}
