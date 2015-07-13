checkInteger <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
}

hex2psrgb<-function(col){
    col<-col2rgb(col)
    col<-col/255
    col<-paste(col,collapse=" ")
    col
}

motifStack_private_fontsize <- 72
coloredSymbols <- function(ncha, font, color, rname, fontsize=motifStack_private_fontsize){
    symbols<-list()
    for(i in 1:ncha){
        ps<-paste("%!PS\n/",font," findfont\n",fontsize," scalefont\n",
                  hex2psrgb(color[i])," setrgbcolor\nsetfont\nnewpath\n0 0 moveto\n(",
                  rname[i],") show",sep="")
        psfilename<-tempfile()
        psfilename <- gsub("\\", "/", psfilename, fixed=TRUE)
        # step1 create a ps file
        cat(ps,file=paste(psfilename,".ps",sep=""))
        # step2 convert it by grImport::PostScriptTrace
        PostScriptTrace(paste(psfilename,".ps",sep=""), paste(psfilename,".xml",sep=""))
        # step3 read by grImport::readPicture
        symbols[[i]]<-readPicture(paste(psfilename,".xml",sep=""))
        unlink(c(paste(psfilename,".ps",sep=""), 
                 paste("capture",basename(psfilename),".ps",sep=""), 
                 paste(psfilename,".xml",sep="")))
    }
    symbols
}

addPseudolog2<-function(x, bgNoise=0.000001){
    #ifelse(x==0, -10, log2(x)) ## !important -10.
    rdirichlet <- function (n, alpha){
        x <- rgamma(n, alpha)
        sm <- x %*% rep(1, n)
        x/as.vector(sm)
    }
    if(any(x==0)){
        x <- (1-bgNoise)*x + 
            bgNoise*rdirichlet(length(x), 1)
    }
    
    log2(x)
}

## get Information Entropy from matrix
getIE<-function(x){
    addPseudolog2(nrow(x))
}

UngappedAlignment<-function(pfms, i, threshold, minimalConsensus=0, rcpostfix="(RC)", revcomp=TRUE){
    if(class(pfms[[i]])!="pfm"){
        pcms <- pfms
        pfms <- lapply(pfms, pcm2pfm)
    }else{
        pcms <- NULL
    }
    res<-getAlignedICWithoutGap(pfms[[i-1]], pfms[[i]], threshold, revcomp)
    if(!is.null(pcms)) pfms <- pcms
    if(res$max>=minimalConsensus){
        if(res$rev){
            pfms[[i]]<-matrixReverseComplement(pfms[[i]])
            pfms[[i]]@name<-paste(pfms[[i]]@name, rcpostfix, sep="")
        }
        if(res$offset>0){
            pfms[[i]]<-addBlank(pfms[[i]], res$offset, FALSE)
        }else{
            if(res[1]<0){
                pfms[1:(i-1)]<-lapply(pfms[1:(i-1)],function(.ele) addBlank(.ele, -1*res$offset, FALSE))
            }
        }
    }
    pfms
}
getAlignedICWithoutGap<-function(pfm1, pfm2, threshold, revcomp=TRUE){
    if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    offset1<-getoffsetPosByIC(pfm1, pfm2, threshold)
    if(revcomp){
        pfm3<-matrixReverseComplement(pfm2)
        offset2<-getoffsetPosByIC(pfm1, pfm3, threshold)
    }else{
        offset2<-list(k=0,max=0)
    }
    offset<-0
    rev<-FALSE
    if(offset1$max < offset2$max){
        rev<-TRUE
        max<-offset2$max
        offset<-offset2$k
    }else{
        max<-offset1$max
        offset<-offset1$k
    }
    list(offset=offset,rev=rev,max=max)    
}
getICbyBase<-function(p, pfmj){
    re<-rep(0, 5)
    re[1:4]<-pfmj * (addPseudolog2(pfmj) - addPseudolog2(p))
    re[5]<-sum(re[1:4])
    re
}
getALLRbyBase <- function(b, pfmi, pfmj){##return 2x ALLR
    sum(pfmj * (addPseudolog2(pfmi) - addPseudolog2(b))) + 
        sum(pfmi * (addPseudolog2(pfmj) - addPseudolog2(b)))
}
checkALLR <- function(...){
    getALLRbyBase(...)>0
}
SWU <- function(pattern, subject, b, 
                match=1, mismatch=-1, gap=-1000){## motif length will never longer than 1000
    if(class(pattern)!="matrix" || class(subject)!="matrix"){
        stop("pattern and subject must be numeric matrix")
    }
    m <- ncol(pattern)
    n <- ncol(subject)
    score <- matrix(0, nrow=m+1, ncol=n+1)
    point <- matrix("none", nrow=m+1, ncol=n+1)
    max_i <- 1
    max_j <- 1
    max_score <- 0
    for(i in 1:m){
        for(j in 1:n){
            diagonal_score <- 0
            left_score <- 0
            up_score <- 0
            #calculate match score
            col1 <- subject[, j]
            col2 <- pattern[, i]
            #score by KLLR
            if(checkALLR(b, col1, col2)){
                diagonal_score <- score[i, j] + match
            }else{
                diagonal_score <- score[i, j] + mismatch
            }
            
            #calculate gap scores
            up_score <- score[i, j+1] + gap
            left_score <- score[i+1, j] + gap
            
            if(diagonal_score <=0 && up_score <=0 && left_score <=0){
                score[i+1, j+1] <- 0
                point[i+1, j+1] <- "none"
                next # terminate this iteration of the loop
            }
            
            #choose best score
            if(diagonal_score >= up_score){
                if(diagonal_score >= left_score){
                    score[i+1, j+1] <- diagonal_score
                    point[i+1, j+1] <- "diagonal"
                }else{
                    score[i+1, j+1] <- left_score
                    point[i+1, j+1] <- "left"
                }
            }else{
                if(up_score > left_score){
                    score[i+1, j+1] <- up_score
                    point[i+1, j+1] <- "up"
                }else{
                    score[i+1, j+1] <- left_score
                    point[i+1, j+1] <- "left"
                }
            }
            
            #set maximum score
            if(score[i+1, j+1] > max_score){
                max_i <- i+1
                max_j <- j+1
                max_score <- score[i+1, j+1]
            }
        }
    }
    #trace-back
    align1 <- c()
    align2 <- c()
    
    i <- max_i
    j <- max_j
    while(1){
        if(point[i, j]=="none"){
            break
        }
        
        if(point[i, j]=="diagonal"){
            align1 <- c(align1, j-1)
            align2 <- c(align2, i-1)
            i <- i-1
            j <- j-1
        }else{
            if(point[i, j]=="left"){
                align1 <- c(align1, j-1)
                align2 <- c(align2, 0)
                j <- j-1
            }else{
                if(point[i, j]=="up"){
                    align1 <- c(align1, 0)
                    align2 <- c(align2, i-1)
                    i <- i-1
                }
            }
        }
    }
    align1 <- rev(align1)
    align2 <- rev(align2)
    if(length(align1)<1 | length(align2)<1){
        k <- 0
    }else{
        k=align1[1]-align2[1]
    }
    list(k=k, max=max_score)
}
ungappedScore <- function(query, subject, b){
    if(class(query)!="matrix" || class(subject)!="matrix"){
        stop("query and subject must be numeric matrix")
    }
    m <- ncol(query)
    n <- ncol(subject)
    score <- matrix(0, nrow=m+1, ncol=n+1)
    for(i in 1:m){
        for(j in 1:n){
            value <- getALLRbyBase(b, query[, i], subject[, j])
            value <- as.numeric(value > 0)
            score[i+1, j+1] <- score[i, j] + value
        }
    }
    score[1:m, 1:n] <- 0
    max_score <- max(score)
    if(max_score==0){
        return(list(k=0, max=0))
    }
    idx <- which(score==max_score, arr.ind=TRUE)
    
    if(nrow(idx)>1){
        idxn <- apply(idx, 1, min)
        idx <- idx[which(idxn==min(idxn)), , drop=FALSE]
    }
    max_r <- idx[, "row"]
    max_c <- idx[, "col"]
    k <- max_r - max_c
    list(k=as.numeric(min(k)), max=max_score)
}
getoffsetPosByIC_dev<-function(pfm1, pfm2, threshold){
    if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    colnames(pfm1$mat) <- paste("V", 1:ncol(pfm1$mat), sep="_")
    colnames(pfm2$mat) <- paste("V", 1:ncol(pfm2$mat), sep="_")
    pfm1.cp <- trimMotif(pfm1, threshold)
    pfm2.cp <- trimMotif(pfm2, threshold)
    if(class(pfm1.cp)!="pfm" | class(pfm2.cp)!="pfm"){
        return(list(k=0, max=0))
    }
    pfm1.left <- which(colnames(pfm1$mat) %in% colnames(pfm1.cp$mat))[1]-1
    pfm2.left <- which(colnames(pfm2$mat) %in% colnames(pfm2.cp$mat))[1]-1
    res <- ungappedScore(pfm1.cp$mat, pfm2.cp$mat, b=pfm1$background)
    if(res$max>0){
        res$k <- res$k - pfm2.left + pfm1.left
    }
    res
}
getoffsetPosByIC<-function(pfm1, pfm2, threshold){
    if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    score1<-rep(0, ncol(pfm1@mat))
    score2<-rep(0, ncol(pfm2@mat))
    for(i in 1:ncol(pfm1@mat)){
        J<-ifelse(ncol(pfm1@mat)+1-i>ncol(pfm2@mat),ncol(pfm2@mat),ncol(pfm1@mat)+1-i)
        for(j in 1:J){
            ic1<-getICbyBase(pfm1@background, pfm1@mat[,i+j-1])
            ic2<-getICbyBase(pfm2@background, pfm2@mat[,j])
            ic3<-getALLRbyBase(pfm1@background, pfm1@mat[,i+j-1], pfm2@mat[,j])
            if(ic1[5]>threshold & ic2[5]>threshold){
                a<-ic1[1:4]>mean(ic1[1:4])
                b<-ic2[1:4]>mean(ic2[1:4])
                if(any((a&b)) && ic3>0.25) score1[i]<-score1[i]+1
                #if(ic3>0) score1[i]<-score1[i]+1
            }
        }
    }
    for(i in 1:ncol(pfm2@mat)){
        J<-ifelse(ncol(pfm2@mat)+1-i>ncol(pfm1@mat),ncol(pfm1@mat),ncol(pfm2@mat)+1-i)
        for(j in 1:J){
            ic2<-getICbyBase(pfm2@background, pfm2@mat[,i+j-1])
            ic1<-getICbyBase(pfm1@background, pfm1@mat[,j])
            ic3<-getALLRbyBase(pfm1@background, pfm1@mat[,j], pfm2@mat[,i+j-1])
            if(ic1[5]>threshold & ic2[5]>threshold){
                a<-ic1[1:4]>mean(ic1[1:4])
                b<-ic2[1:4]>mean(ic2[1:4])
                if(any((a&b)) && ic3>0.25) score2[i]<-score2[i]+1
                #if(ic3>0) score2[i]<-score2[i]+1
            }
        }
    }
    k1<-match(max(score1),score1)
    k2<-match(max(score2),score2)
    k=ifelse(score1[k1]<score2[k2], -1*(k2-1), (k1-1))
    max=ifelse(score1[k1]<score2[k2], score2[k2], score1[k1])
    list(k=k,max=max)
}

setClass("Rect", 
         representation(x="numeric", y="numeric", width="numeric", height="numeric"), 
         prototype=prototype(x=0, y=0, width=0, height=0)
)

setGeneric("isContainedIn", function(a, b) standardGeneric("isContainedIn"))

setMethod("isContainedIn", signature(a="Rect", b="Rect"), function(a, b){
    a@x >= b@x && a@y >= b@y && a@x+a@width <= b@x+b@width && a@y+a@height <= b@y+b@height
})

setClass("pos",
         representation(box="Rect", beta="numeric", sig="pfm", freq="numeric", norm="numeric")
)

getPFMid <- function(pfms, nodename, rcpostfix="(RC)"){
    pfmNames <- as.character(unlist(lapply(pfms, function(.ele) .ele@name)))
    pfmNames <- gsub(rcpostfix,"",pfmNames, fixed=TRUE)
    which(pfmNames==nodename)
}

getParentNode <- function(nodelist, nodename){
    for(i in 1:length(nodelist)){
        currNode <- nodelist[[i]]
        if(currNode@left==nodename) return(c(currNode@parent, "left"))
        if(currNode@right==nodename) return(c(currNode@parent, "right"))
    }
    NULL
}

makeLeaveNames <- function(ch){
    gsub(".", "_", make.names(ch), fixed=TRUE)
}