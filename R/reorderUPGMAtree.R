#' re-order UPGMA tree
#' 
#' re-order the UPGMA tree by adjacent motif distance
#' 
#' 
#' @param phylog an object of phylog
#' @param motifs a list of objects of pfm
#' @param rcpostfix the postfix for reverse complements
#' @return an object of phylog
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @examples
#' 
#'   if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'     library("MotifDb")
#'     matrix.fly <- query(MotifDb, "Dmelanogaster")
#'     motifs <- as.list(matrix.fly)
#'     motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", names(motifs), fixed=TRUE)]
#'     names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
#'                 gsub("_FBgn[0-9]+$", "", 
#'                   gsub("[^a-zA-Z0-9]","_", 
#'                      gsub("(_[0-9]+)+$", "", names(motifs)))))
#'     motifs <- motifs[unique(names(motifs))]
#'     pfms <- sample(motifs, 50)
#'     hc <- clusterMotifs(pfms)
#'     library(ade4)
#'     phylog <- ade4::hclust2phylog(hc)
#'     pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
#'                  new("pfm",mat=.ele, name=.name)})
#'     reorderUPGMAtree(phylog, pfms)
#'   }
#' 
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
    d <- matalign(motifs)
    d <- matalignOut2dist(d)
    motifs <- lapply(motifs, function(.ele){
        if(is(.ele, "pcm")) .ele <- pcm2pfm(.ele)
        .ele
    })
    
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
