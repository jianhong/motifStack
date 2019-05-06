setClass("pfm", 
        representation(mat="matrix", name="character", alphabet="character", color="character", background="numeric"),
        validity=function(object){
            re<-TRUE
            if (any(abs(1 - apply(object@mat, 2, sum)) > 0.01)) re<-"Columns of PFM must add up to 1.0"
            if (is.null(rownames(object@mat))) re<-"rownames of PFM is empty"
            if (is.null(names(object@color))) re<-"names of color is empty"
            if (length(setdif(rownames(object@mat), names(object@color)))>0) re<-"not every alphabet in pfm has a corresponding color"
            if (is.null(names(object@background))) re<-"names of background is empty"
            if (length(setdif(rownames(object@mat), names(object@background)))>0) re<-"not every alphabet in pfm has a background"
            if (!object@alphabet %in% c("DNA", "RNA", "AA", "others")) re<-"alphabet must be either DNA, RNA, AA or others"
            if (object@alphabet=="DNA" & !all(rownames(object@mat) %in% c("A","C","G","T"))) re<-"rownames of PFM for DNA motifs must be capital A, C, G and T"
            if (object@alphabet=="RNA" & !all(rownames(object@mat) %in% c("A","C","G","U"))) re<-"rownames of PFM for RNA motifs must be capital A, C, G and U"
            if (object@alphabet=="AA" & !all(rownames(object@mat) %in% LETTERS[1:26])) re<-"rownames of PFM for AA motifs must be capital letters"
            re
        }
)

setMethod("$", "pfm", function(x, name) slot(x, name))
setReplaceMethod("$", "pfm",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

setMethod("initialize","pfm",function(.Object, mat, name, alphabet, color, background){
    if(mode(mat)!="numeric") stop("mat must be a numeric matrix")
    if(is.null(rownames(mat))) stop("rownames of PFM is empty")
    rownames(mat)<-toupper(rownames(mat))
    rname<-rownames(mat)
    if(missing(alphabet)){
        rname<-rname[order(rname)]
        if(all(rname %in% c("A","C","G","T"))) alphabet<-"DNA"
        else{
            if(all(rname %in% c("A","C","G","U"))) alphabet<-"RNA"
            else{
                if(all(rname %in% LETTERS[1:26]) & nrow(mat)==20) alphabet<-"AA"
                else alphabet<-"others"
            }
        }
    }
    if(missing(color)){
      if(alphabet!="others"){
        color<-colorset(alphabet, "auto")
      }else{
        color <- rainbow(nrow(mat))
        names(color) <- rname
      }
    }
    names(color) <- toupper(names(color))
    if(missing(background)){
        background<-rep(1/nrow(mat), nrow(mat))
        names(background)<-rname
    }
    if(length(background)<nrow(mat)){
        background <- rep(background, nrow(mat))[1:nrow(mat)]
        names(background) <- rname
    }
    .Object@mat            =    mat;
    .Object@name        =    name;
    .Object@alphabet    =    alphabet;
    .Object@color        =    color;
    .Object@background    =    background;
    .Object
})

setMethod("plot", signature(x="pfm", y="ANY"), 
    function(x, y="missing", ...){
        plotMotifLogo(pfm=x@mat, 
                      motifName=x@name, 
                      p=x@background[rownames(x@mat)], 
                      colset=x@color[rownames(x@mat)], 
                      ...)
    }
)

setAs(from="pfm", to="matrix", function(from){
    from@mat
})

## get information content profile from PFM
setMethod("getIC", signature(x="pfm"), function(x, p="missing"){
    colSums(x@mat * (addPseudolog2(x@mat) - addPseudolog2(x@background)))
})
## get information content profile from matrix
setMethod("getIC", signature(x="matrix", p="numeric"), function(x, p){
    if (any(abs(1 - apply(x, 2, sum)) > 0.01)) re<-"Columns of PFM must add up to 1.0"
    colSums(x * (addPseudolog2(x) - addPseudolog2(p)))
})
setMethod("getIC", signature(x="matrix", p="matrix"), function(x, p){
    if (any(abs(1 - apply(x, 2, sum)) > 0.01)) re<-"Columns of PFM must add up to 1.0"
    colSums(x * (addPseudolog2(x) - addPseudolog2(p)))
})

setMethod("trimMotif", signature(x="pfm", t="numeric"), function(x, t){
    .sw <- c(0, 0)
    ic <- getIC(x)
    for(i in 1:ncol(x@mat)){
        if(ic[i]>=t){
            .flag <- FALSE
            if(i<ncol(x@mat)){
                if(ic[i+1]>=t) .flag<-TRUE
                else if(ic[i]>=min(2*t, 0.5)) .flag<-TRUE
            }else{
                .flag<-TRUE
            }
            if(.flag){
                if(.sw[1]==0) .sw[1] <- i
                else .sw[2] <- i
            }
        }
    }
    if(.sw[1]<.sw[2]-1) x@mat <- x@mat[, .sw[1]:.sw[2]]
    else x <- NA
    x
})

setMethod("matrixReverseComplement", "pfm", function(x){
    if(x@alphabet!="DNA") stop("alphabet of pfm must be DNA")
    mat<-x@mat
    rc<-matrix(nrow=nrow(mat),ncol=ncol(mat))
    complements<-c(4,3,2,1)
    for(i in 1:nrow(mat)){
        for(j in 1:ncol(mat)){
            rc[complements[i],(ncol(mat)+1-j)]<-mat[i,j]
        }
    }
    rownames(rc)<-rownames(mat)
    x@mat<-rc
    x
})

setMethod("addBlank", signature(x="pfm", n="numeric", b="logical"), function(x, n, b){
    if(x@alphabet!="DNA") stop("alphabet of pfm must be DNA")
    N<-matrix(rep(x@background, n), nrow=4)
    if(b){
        N<-cbind(x@mat,N)
    }else{
        N<-cbind(N,x@mat)
    }
    rownames(N)<-rownames(x@mat)
    x@mat<-N
    x
})

## for data.frame
setMethod("as.data.frame", signature(x="pfm"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@mat, ...)
})

setMethod("format", signature(x="pfm"), function(x, ...){
  paste0(x@name, "_pfm")
})
