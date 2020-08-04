#' Class \code{"pcm"}
#' 
#' An object of class \code{"pcm"} represents the position count matrix of a
#' DNA/RNA/amino-acid sequence motif. The entry stores a matrix, which in row
#' i, column j gives the counts of observing nucleotide/or amino acid i in
#' position j of the motif.
#' 
#' 
#' @name pcm-class
#' @aliases pcm
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("pcm", mat, name, alphabet, color, background)}.
#' @keywords classes
#' @export
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
#' plot(motif)
#' pcm2pfm(pcm)
#' pcm2pssm(pcm)
setClass(
  "pcm", 
  representation(mat="matrix", name="character", alphabet="character", 
                 color="character", background="numeric", tags="list",
                 markers="list"),
  validity=function(object){
    re<-TRUE
    if (checkInteger(unlist(object@mat))) re<-"Columns of pcm must be integer"
    if (is.null(rownames(object@mat))) re<-"rownames of pcm is empty"
    if (is.null(names(object@color))) re<-"names of color is empty"
    if (length(setdif(rownames(object@mat), names(object@color)))>0) re<-"not every alphabet in pcm has a corresponding color"
    if (is.null(names(object@background))) re<-"names of background is empty"
    if (length(setdif(rownames(object@mat), names(object@background)))>0) re<-"not every alphabet in pcm has a background"
    if (!object@alphabet %in% c("DNA", "RNA", "AA", "others")) re<-"alphabet must be either DNA, RNA, AA or others"
    if (object@alphabet=="DNA" & !all(rownames(object@mat) %in% c("A","C","G","T"))) re<-"rownames of pcm for DNA motifs must be capital A, C, G and T"
    if (object@alphabet=="RNA" & !all(rownames(object@mat) %in% c("A","C","G","U"))) re<-"rownames of pcm for RNA motifs must be capital A, C, G and U"
    if (object@alphabet=="AA" & !all(rownames(object@mat) %in% LETTERS[1:26])) re<-"rownames of pcm for AA motifs must be capital letters"
    if (any(!is(object@markers, "marker"))) re<- "markers must be a list of marker object"
    re
  }
)


#' "pcm" methods
#' 
#' methods for pcm objects.
#' 
#' 
#' @rdname pcm-class
#' @aliases $,pcm-method $<-,pcm-method 
#' @docType methods
#' @param x An object of class \code{pcm}. For \code{getIC}, if parameter p is
#' followed, x should be an object of matrix.  For \code{pcm2pfm}, x also could
#' be an object of matrix.
#' @param y Not use.
#' @param p p is the background frequency.
#' @param n how many spaces should be added.
#' @param b logical value to indicate where the space should be added.
#' @param background a \code{"numeric"} vector. The background frequency.
#' @param t numeric value of information content threshold for trimming.
#' @param \dots Further potential arguments passed to \code{plotMotifLogo}.
#' @param row.names,optional see as.data.frame
#' @param name slot name of pcm object.
#' @section Methods: \describe{ \item{addBlank}{\code{signature(x="pcm",
#' n="numeric", b="logical")} add space into the position count matrix for
#' alignment. b is a bool value, if TRUE, add space to the 3' end, else add
#' space to the 5' end. n indicates how many spaces should be added.}
#' 
#' \item{coerce}{\code{signature(from = "pcm", to = "matrix")}: convert object
#' pcm to matrix }
#' 
#' \item{getIC}{\code{signature(x = "pcm",)} Calculate information content
#' profile for position frequency matrix. }
#' 
#' \item{matrixReverseComplement}{\code{signature(x = "pcm")} get the reverse
#' complement of position frequency matrix.}
#' 
#' \item{plot}{\code{signature(x = "pcm")} Plots the sequence logo of the
#' position count matrix. }
#' 
#' \item{trimMotif}{\code{signature(x = "pcm", t= "numeric")} trim motif by
#' information content. }
#' 
#' \item{$, $<-}{Get or set the slot of \code{\link{pcm-class}}}
#' 
#' \item{as.data.frame}{convert \code{\link{pcm-class}} to a data.frame}
#' \item{format}{return the name_pcm of \code{\link{pcm-class}}} }
#' @keywords classes
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
#' getIC(motif)
#' matrixReverseComplement(motif)
#' as(motif,"matrix")
#' pcm2pfm(motif)
#' as.data.frame(motif)
#' format(motif)
#' 
#' @exportMethod `$` `$<-`
setMethod("$", "pcm", function(x, name) slot(x, name))
setReplaceMethod("$", "pcm",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

#' @importFrom grDevices rainbow
setMethod("initialize","pcm",function(.Object, mat, name, alphabet, color, background, 
                                      tags, markers){
    if(mode(mat)!="numeric") stop("mat must be a numeric matrix")
    if(is.null(rownames(mat))) stop("rownames of pcm is empty")
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

#' @name plot
#' @rdname pcm-class
#' @aliases plot plot,pcm,ANY-method
#' @importFrom stats4 plot
#' @exportMethod plot
setMethod("plot", signature(x="pcm"), 
          function(x, y="missing", ...){
              plot(pcm2pfm(x), ...)
          }
)

#' @name as
#' @rdname pcm-class
#' @aliases coerce,pcm,matrix-method
#' @exportMethod coerce
setAs(from="pcm", to="matrix", function(from){
    from@mat
})



#' @rdname pcm-class
#' @exportMethod trimMotif
#' @aliases trimMotif trimMotif,pcm,numeric-method
setGeneric("trimMotif", function(x, t) standardGeneric("trimMotif"))
setMethod("trimMotif", signature(x="pcm", t="numeric"), function(x, t){
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

#' @rdname pcm-class
#' @exportMethod matrixReverseComplement
#' @aliases matrixReverseComplement matrixReverseComplement,pcm-method
setGeneric("matrixReverseComplement", 
           function(x) standardGeneric("matrixReverseComplement"))
setMethod("matrixReverseComplement", "pcm", function(x){
    if(x@alphabet!="DNA") stop("alphabet of pcm must be DNA")
    mat<-x@mat
    ## double check the rownames to A, C, G, T.
    if(!all(toupper(rownames(mat))==c("A", "C", "G", "T"))){
      mat <- mat[match(c("A", "C", "G", "T"), toupper(rownames(mat))), ]
    }
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

#' @rdname pcm-class
#' @exportMethod addBlank
#' @aliases addBlank addBlank,pcm,numeric,logical-method
setGeneric("addBlank", function(x, n, b) standardGeneric("addBlank"))
setMethod("addBlank", signature(x="pcm", n="numeric", b="logical"), function(x, n, b){
    #if(x@alphabet!="DNA") stop("alphabet of pfm must be DNA")
    N<-matrix(rep(rep(0,nrow(x@mat)), n), nrow=nrow(x@mat))
    rownames(N) <- rownames(x@mat)
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

#' @rdname pcm-class
#' @exportMethod getIC
#' @aliases getIC getIC,pcm,ANY-method
setGeneric("getIC", function(x, p) standardGeneric("getIC"))
setMethod("getIC", signature(x="pcm"), function(x, p="missing"){
    mat<-pcm2pfm(x@mat)
    colSums(mat * (addPseudolog2(mat) - addPseudolog2(x@background)))
})

#' @rdname pcm-class
#' @exportMethod pcm2pfm
#' @aliases pcm2pfm pcm2pfm,pcm,ANY-method pcm2pfm,matrix,ANY-method
#' pcm2pfm,matrix,numeric-method pcm2pfm,data.frame,ANY-method
#' pcm2pfm,data.frame,numeric-method 
setGeneric("pcm2pfm", function(x, background) standardGeneric("pcm2pfm"))
setMethod("pcm2pfm", signature(x="matrix", background="numeric"), function(x, background){
    s<-apply(x,2,sum)
    s[s==0] <- NA
    p<-apply(x,1,function(.ele) .ele/s)
    n<-length(s[is.na(s)])
    bck <- background/sum(background)
    if(any(is.na(s))) {
      p[is.na(p[,1]),] <- rep(bck, each=n)
    }
    t(p)
})

setMethod("pcm2pfm", signature(x="matrix"), function(x, background="missing"){
    pcm2pfm(x, rep(1/nrow(x), nrow(x)))
})

setMethod("pcm2pfm", signature(x="data.frame", background="numeric"), function(x, background){
    pcm2pfm(as.matrix(x), background)
})

setMethod("pcm2pfm", signature(x="data.frame"), function(x, background="missing"){
    pcm2pfm(as.matrix(x), rep(1/nrow(x), nrow(x)))
})

## pcm to pfm from pcm object
setMethod("pcm2pfm", signature(x="pcm"), function(x, background="missing"){
  bck <- x@background[rownames(x@mat)]
  if(any(is.na(bck))) stop("some elements does not have background value", bck)
    .mat <- pcm2pfm(x@mat, bck)
    new("pfm", mat=.mat, name=x@name, alphabet=x@alphabet, 
        color=x@color, background=x@background, tags=x@tags,
        markers=x@markers)
})

#' @rdname pcm-class
#' @exportMethod pcm2pssm
#' @aliases pcm2pssm pcm2pssm,pcm,ANY-method pcm2pssm,matrix,ANY-method
#' pcm2pssm,matrix,numeric-method pcm2pssm,data.frame,ANY-method
#' pcm2pssm,data.frame,numeric-method 
setGeneric("pcm2pssm", function(x, background) standardGeneric("pcm2pssm"))
setMethod("pcm2pssm", signature(x="matrix", background="numeric"), function(x, background){
  ## check the background order, first order, length=4; second order, length=16
  ord <- 0
  if(length(background)==nrow(x)){
    ord <- 1
  }
  if(length(background)==2^nrow(x)+nrow(x)){
    ord <- 2
  }
  if(ord==0){
    stop("The background dimention is incorrect. For example, if it is DNA motif,
         4 or 20 (16 (2nd order) + 4(1st order)) named backgrounds are expected.")
  }
  m <- pcm2pfm(x, background)
  stopifnot(all(rownames(m)==rownames(x)))
  if(length(names(background))!=length(background)){
    stop("background must be a named vector.")
  }
  if(ord==1){
    m1 <- addPseudolog2(m) - addPseudolog2(background[rownames(m)])
  }else{
    ## 2 order background
    ## 2*log2(sum(p_i*(p_ij/q_ij)))
    letters <- rownames(m)
    ## calculate 1 order of the motif
    total <- sum(x, na.rm = TRUE)
    P1 <- rowSums(x, na.rm = TRUE)
    names(P1) <- letters
    P1 <- addPseudolog2(P1)
    ## calculate 2 order of the motif
    P2 <- lapply(letters, function(.ele){
      mat <- t(t(cbind(background[rownames(m)], m[, -ncol(m)])) *
                 m[.ele, , drop=TRUE]) 
      rownames(mat) <- paste0(letters, .ele)
      mat
    })
    names(P2) <- letters
    ## 2 order matrix
    bg <- split(background[paste0(rep(letters, length(letters)),
                                  rep(letters, each=length(letters)))], 
                rep(letters, each=length(letters)))[letters]
    Pq <- mapply(FUN=function(p, b){
      addPseudolog2(p) - addPseudolog2(b)
    }, P2[letters], bg[letters], SIMPLIFY = FALSE)
    m1 <- lapply(Pq, function(.ele){## .ele is a matrix, eg 4 X len
      colSums(.ele + P1)
    })
    m1 <- do.call(rbind, m1)
    rownames(m1) <- letters
  }
  return(m1)
})

setMethod("pcm2pssm", signature(x="matrix"), function(x, background="missing"){
  background <- rep(1/nrow(x), nrow(x))
  names(background) <- rownames(x)
  pcm2pssm(x, background)
})

setMethod("pcm2pssm", signature(x="data.frame", background="numeric"), function(x, background){
  pcm2pssm(as.matrix(x), background)
})

setMethod("pcm2pssm", signature(x="data.frame"), function(x, background="missing"){
  background <- rep(1/nrow(x), nrow(x))
  names(background) <- rownames(x)
  pcm2pssm(as.matrix(x), background)
})
## pcm to pssm from pcm object
setMethod("pcm2pssm", signature(x="pcm"), function(x, background="missing"){
  .mat <- pcm2pssm(x@mat, x@background)
  new("pssm", mat=.mat, name=x@name, alphabet=x@alphabet, 
      color=x@color, background=x@background, tags=x@tags,
      markers=x@markers)
})


#' @rdname pcm-class
#' @exportMethod as.data.frame
#' @aliases as.data.frame,pcm-method 
setMethod("as.data.frame", signature(x="pcm"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@mat, ...)
})

#' @rdname pcm-class
#' @exportMethod format
#' @aliases format,pcm-method
setMethod("format", signature(x="pcm"), function(x, ...){
  paste0(x@name, "_pcm")
})
