#' Class \code{"pfm"}
#' 
#' An object of class \code{"pfm"} represents the position frequency matrix of
#' a DNA/RNA/amino-acid sequence motif. The entry stores a matrix, which in row
#' i, column j gives the frequency of observing nucleotide/or amino acid i in
#' position j of the motif.
#' 
#' 
#' @name pfm-class
#' @aliases pfm
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("pfm", mat, name, alphabet, color, background)}.
#' @keywords classes
#' @export
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- pcm2pfm(pcm)
#' motif <- new("pfm", mat=motif, name="bin_SOLEXA")
#' plot(motif)
#' 
setClass("pfm", 
        representation(mat="matrix", name="character", alphabet="character", 
                       color="character", background="numeric", tags="list",
                       markers="list"),
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
            if (any(!is(object@markers, "marker"))) re<- "markers must be a list of marker object"
            re
        }
)

#' "pfm" methods
#' 
#' methods for pfm objects.
#' 
#' 
#' @rdname pfm-class
#' @aliases $,pfm-method $<-,pfm-method
#' @docType methods
#' @param x An object of class \code{pfm}. For \code{getIC}, if parameter p is
#' followed, x should be an object of matrix.
#' @param name Slot name.
#' @param y Not use.
#' @param p p is the background frequency.
#' @param n how many spaces should be added.
#' @param b logical value to indicate where the space should be added.
#' @param t numeric value of information content threshold for trimming.
#' @param \dots Further potential arguments passed to \code{plotMotifLogo}.
#' @param row.names,optional see as.data.frame
#' @section Methods: \describe{ \item{addBlank}{\code{signature(x="pfm",
#' n="numeric", b="logical")} add space into the position frequency matrix for
#' alignment. b is a bool value, if TRUE, add space to the 3' end, else add
#' space to the 5' end. n indicates how many spaces should be added.}
#' 
#' \item{getIC}{\code{signature(x = "pfm",)} Calculate information content
#' profile for position frequency matrix. }
#' 
#' \item{getIC}{\code{signature(x = "matrix", p = "numeric")} Calculate
#' information content profile for matrix. p is the background frequency}
#' 
#' \item{matrixReverseComplement}{\code{signature(x = "pfm")} get the reverse
#' complement of position frequency matrix.}
#' 
#' \item{plot}{\code{signature(x = "pfm")} Plots the sequence logo of the
#' position frequency matrix. }
#' 
#' \item{trimMotif}{\code{signature(x = "pfm", t= "numeric")} trim motif by
#' information content. }
#' 
#' \item{$, $<-}{Get or set the slot of \code{\link{pfm-class}}}
#' \item{as.data.frame}{convert \code{\link{pfm-class}} to a data.frame}
#' \item{format}{return the name_pfm of \code{\link{pfm-class}}} }
#' @keywords classes
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- pcm2pfm(pcm)
#' motif <- new("pfm", mat=motif, name="bin_SOLEXA")
#' getIC(motif)
#' matrixReverseComplement(motif)
#' addBlank(motif, 1, FALSE)
#' addBlank(motif, 3, TRUE)
#' as(motif,"matrix")
#' as.data.frame(motif)
#' format(motif)
#' 
#' @exportMethod `$` `$<-`
setMethod("$", "pfm", function(x, name) slot(x, name))
setReplaceMethod("$", "pfm",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })
#' @importFrom grDevices rainbow
setMethod("initialize","pfm",function(.Object, mat, name, alphabet, color, background, tags, markers){
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
    .Object@mat            =    mat
    .Object@name        =    name
    .Object@alphabet    =    alphabet
    .Object@color        =    color
    .Object@background    =    background
    if(!missing(tags)) .Object@tags = tags
    if(!missing(markers)) .Object@markers = markers
    .Object
})

#' @rdname pfm-class 
#' @aliases plot,pfm,ANY-method
#' @exportMethod plot
setMethod("plot", signature(x="pfm"), 
    function(x, y="missing", ...){
        plotMotifLogo(pfm=x, ...)
    }
)

#' @name coerce
#' @rdname pfm-class
#' @aliases coerce,pfm,matrix-method
#' @exportMethod coerce
setAs(from="pfm", to="matrix", function(from){
    from@mat
})

## get information content profile from PFM
#' @rdname pfm-class
#' @exportMethod getIC
#' @aliases getIC,pfm,ANY-method getIC,matrix,numeric-method 
#' getIC,matrix,matrix-method
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

#' @rdname pfm-class
#' @exportMethod trimMotif
#' @aliases trimMotif,pfm,numeric-method
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

#' @rdname pfm-class
#' @exportMethod matrixReverseComplement
#' @aliases matrixReverseComplement,pfm-method
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

#' @rdname pfm-class
#' @exportMethod addBlank
#' @aliases addBlank,pfm,numeric,logical-method
setMethod("addBlank", signature(x="pfm", n="numeric", b="logical"), function(x, n, b){
    if(x@alphabet!="DNA") stop("alphabet of pfm must be DNA")
    N<-matrix(rep(x@background, n), nrow=4)
    if(b){
        N<-cbind(x@mat,N)
    }else{
        N<-cbind(N,x@mat)
        if(length(x@markers)>0){
          x@markers <- lapply(x@markers, function(.ele){
            .ele@start <- .ele@start+n
            .ele@stop <- .ele@stop+n
            .ele
          })
        }
    }
    rownames(N)<-rownames(x@mat)
    x@mat<-N
    x
})

## for data.frame
#' @rdname pfm-class
#' @exportMethod as.data.frame
#' @aliases as.data.frame,pfm-method
setMethod("as.data.frame", signature(x="pfm"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@mat, ...)
})

#' @rdname pfm-class
#' @exportMethod format
#' @aliases format,pfm-method
setMethod("format", signature(x="pfm"), function(x, ...){
  paste0(x@name, "_pfm")
})
