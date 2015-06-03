reorderUPGMAtree <- function(phylog, motifs, rcpostfix = "(RC)"){
    if(!inherits(phylog, "phylog")) stop("phylog must be an object of phylog")
    if(!all(sapply(motifs, function(.ele) inherits(.ele, c("pfm", "pcm")))))
        stop("motifs must be a list of pfm or pcm object")
    lnames <- sapply(motifs, function(.ele) .ele$name)
    for(metaChar in c("\\","$","*","+",".","?","[","]","^","{","}","|","(",")"))
    {
        rcpostfix <- gsub(metaChar,paste("\\",metaChar,sep=""),rcpostfix,fixed=TRUE)
    }
    lnames <- gsub(paste(rcpostfix,"$",sep=""),"",lnames)
    if(!all(names(phylog$leaves) %in% lnames)){
        stop("not all leaves have corresponding motif")
    }
    names(motifs) <- lnames
    motifs <- motifs[names(phylog$leaves)]
    motifs <- lapply(motifs, function(.ele){
        if(class(.ele)=="pcm") .ele <- pcm2pfm(.ele)
        .ele
    })
    
    jaspar.scores <- readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
    d <- motifDistances(lapply(motifs, pfm2pwm), DBscores=jaspar.scores)
    d <- as.matrix(d)
    
    paths <- phylog$paths
    par <- sapply(paths, function(.ele) .ele[length(.ele)-1])
    
    pars <- split(names(par), factor(par, levels=unique(par)))
    
    
    buildTree <- function(pars, ln, rn){
        re <- list()
        re[[ln]] <- if(is.null(pars[[ln]])) ln else pars[[ln]]
        re[[rn]] <- if(is.null(pars[[rn]])) rn else pars[[rn]]
        if(length(unlist(pars[re[[ln]]]))>0){
            re[[ln]] <- buildTree(pars, re[[ln]][1], re[[ln]][2])
        }
        if(length(unlist(pars[re[[rn]]]))>0){
            re[[rn]] <- buildTree(pars, re[[rn]][1], re[[rn]][2])
        }
        
        return(re)
    }
    
    shortestDist <- function(mat){
        d <- dim(mat)
        id <- which(mat==min(mat))[1]
        j <- ceiling(id/d[1])
        i <- id %% d[1]
        if(i==0) i <- d[1]
        paste(i, j, sep="")
    }
    
    sortTree <- function(tree){
        child1 <- tree[[1]]
        child2 <- tree[[2]]
        if(is.list(child1)){
            child1 <- sortTree(child1)
        }
        if(is.list(child2)){
            child2 <- sortTree(child2)
        }
        child1.flat <- unlist(child1)
        len1 <- length(child1.flat)
        child2.flat <- unlist(child2)
        len2 <- length(child2.flat)
        
        dist <- d[child1.flat[c(1, len1)], child2.flat[c(1, len2)], drop=FALSE]
        neighbor <- shortestDist(dist)
        leaves <- switch(neighbor,
                         "11"=c(rev(child1.flat), child2.flat),
                         "12"=c(child2.flat, child1.flat),
                         "21"=c(child1.flat, child2.flat),
                         "22"=c(child1.flat, rev(child2.flat)),
                         c(child1.flat, child2.flat))
    }
    
    tre <- pars[["Root"]]
    if(length(tre)==0){
        return(phylog)
    }
    tre <- buildTree(pars, tre[1], tre[2])
    phylog$leaves <- phylog$leaves[sortTree(tre)]
    phylog
}