getRankedUniqueMotifs <- function(phylog, attr){
    leaves=names(phylog$leaves)
    if(is.null(names(attr))){
        names(attr) <- leaves
    }else{
        attr <- attr[leaves]
    }
    droot <- phylog$droot
    droot.max <- max(droot)
    parts <- phylog$parts
    parts.dist <- droot.max - droot[names(parts)]
    parts.dist <- sort(parts.dist)
    parts.dist <- parts.dist[names(parts.dist)!="Root"]
    
    getLeaves <- function(parts, name){
        p <- parts[[name]]
        y <- c()
        for(x in p){
            if(x %in% names(parts)) y <- c(y, Recall(parts, x))
            else y <- c(y, x)
        }
        return(y)
    }
    
    getNodes <- function(parts, name){
        p <- parts[[name]]
        y <- c()
        for(x in p){
            if(x %in% names(parts)) {
                y <- c(y, x, Recall(parts, x))
            }
        }
        return(y)
    }
    
    uni.dist <- unique(parts.dist)
    
    uni <- lapply(uni.dist, function(gpDist){
        names <- names(parts.dist)[parts.dist <= gpDist]
        groups <- lapply(names, getLeaves, parts=parts)
        names(groups) <- names
        toberemoved <- lapply(names, getNodes, parts=parts)
        toberemoved <- unique(unlist(toberemoved))
        if(length(toberemoved)>0) groups <- groups[names[!names %in% toberemoved]]
        groupedLeaves <- unlist(groups)
        ungroupedLeaves <- leaves[!leaves %in% groupedLeaves]
        groups <- c(groups, as.list(ungroupedLeaves))
        if(length(groups)<1){
            return(NULL)
        }
        groups <- do.call(rbind, mapply(function(.ele, .id) cbind(.ele, .id), groups, 1:length(groups)))
        category <- attr[groups[, 1]]
        category <- split(category, groups[,2])
        category <- lapply(category, unique)
        uni.names <- names(category)[sapply(category, length)==1]
        names <- groups[groups[, 2] %in% uni.names, 1]
        as.character(names)
    })
    
    uni.rank <- table(unlist(uni))
    uni.rank <- uni.rank[order(uni.rank, decreasing=TRUE)]
    uni.rank.value <- rle(as.numeric(uni.rank))
    uni.rank.value$values <- c(1, cumsum(uni.rank.value$lengths)[-length(uni.rank.value$lengths)]+1)
    uni.rank.value <- inverse.rle(uni.rank.value)
    names(uni.rank.value) <- names(uni.rank)
    
    uni.length <- sapply(uni, length)
    uni.length <- cbind(distance=uni.dist, uniqueMotifCount=uni.length)
    uni.length <- uni.length[order(-uni.length[,2], -uni.length[,1]),]
    uni.length <- uni.length[!duplicated(uni.length[,2]),]
    
    names(uni) <- paste("dist", uni.dist, sep="_")
    uni <- uni[paste("dist", uni.length[,1], sep="_")]
    
    return(list(uni.rank=uni.rank.value, uni.length=uni.length, uni.list=uni))
}